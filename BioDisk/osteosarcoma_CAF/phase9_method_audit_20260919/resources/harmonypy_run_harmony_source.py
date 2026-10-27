def run_harmony(
    data_mat: np.ndarray,
    meta_data,
    vars_use,
    theta=None,
    lamb=None,
    sigma=0.1,
    nclust=None,
    tau=0,
    block_size=0.05,
    max_iter_harmony=10,
    max_iter_kmeans=4,
    epsilon_cluster=1e-3,
    epsilon_harmony=1e-2,
    alpha=0.2,
    batch_prop_cutoff=1e-5,
    verbose=True,
    random_state=0,
    ncores=0,
):
    """Run Harmony batch effect correction.

    Parameters
    ----------
    data_mat : np.ndarray
        PCA embedding matrix (cells x PCs or PCs x cells)
    meta_data : dict-like
        Metadata with batch variables (cells x variables).
        Accepts pandas DataFrame, dict of arrays, or any object
        supporting ``meta_data[var]`` column access and ``.shape[0]``.
    vars_use : str or list
        Column name(s) in meta_data to use for batch correction
    theta : float or list, optional
        Diversity penalty parameter(s). Default is 2 for each batch.
    lamb : float or list, optional
        Ridge regression penalty. Default is None (auto-estimation).
        Set to a positive value for fixed lambda.
    sigma : float, optional
        Kernel bandwidth for soft clustering. Default is 0.1.
    nclust : int, optional
        Number of clusters. Default is min(N/30, 100).
    tau : float, optional
        Protection against overcorrection. Default is 0.
    block_size : float, optional
        Proportion of cells to update in each block. Default is 0.05.
    max_iter_harmony : int, optional
        Maximum Harmony iterations. Default is 10.
    max_iter_kmeans : int, optional
        Maximum k-means iterations per Harmony iteration. Default is 4.
    epsilon_cluster : float, optional
        K-means convergence threshold. Default is 1e-3.
    epsilon_harmony : float, optional
        Harmony convergence threshold. Default is 1e-2.
    alpha : float, optional
        Alpha parameter for lambda estimation. Default is 0.2.
    batch_prop_cutoff : float, optional
        Minimum batch proportion in a cluster for correction. Default is 1e-5.
    verbose : bool, optional
        Print progress messages. Default is True.
    random_state : int, optional
        Random seed for reproducibility. Default is 0.
    ncores : int, optional
        Number of BLAS threads for matrix operations. Default is 0
        (use all available cores). Set to 1 for single-threaded execution.

    Returns
    -------
    Harmony
        Harmony object with corrected data in Z_corr attribute.
    """
    # Ensure data_mat is a proper numpy array (drop non-numeric columns)
    if hasattr(data_mat, 'select_dtypes'):
        data_mat = data_mat.select_dtypes(include=[np.number]).values
    elif hasattr(data_mat, 'values'):
        data_mat = data_mat.values

    if hasattr(meta_data, 'shape'):
        N = meta_data.shape[0]
    else:
        # dict-like: get length from first column
        first_key = vars_use[0] if isinstance(vars_use, list) else vars_use
        N = len(meta_data[first_key])
    if data_mat.shape[1] == N:
        pass
    elif data_mat.shape[0] == N:
        data_mat = data_mat.T
    else:
        raise ValueError(
            f"data_mat has shape {data_mat.shape}, but meta_data has {N} cells. "
            f"One dimension of data_mat must match the number of cells."
        )

    if nclust is None:
        nclust = int(min(round(N / 30.0), 100))

    if isinstance(sigma, float) and nclust > 1:
        sigma = np.repeat(sigma, nclust)

    if isinstance(vars_use, str):
        vars_use = [vars_use]

    # Build compact batch-of-cell index (n_covariates x N, int64).
    # Passed directly to C++ — no sparse Phi construction needed.
    batch_of_cell = np.empty((len(vars_use), N), dtype=np.int64)
    phi_n = np.empty(len(vars_use), dtype=int)
    offset = 0
    for c, var in enumerate(vars_use):
        uniques, codes = np.unique(np.asarray(meta_data[var]), return_inverse=True)
        n_levels = len(uniques)
        batch_of_cell[c] = codes + offset
        phi_n[c] = n_levels
        offset += n_levels

    # Theta handling - default is 2 (matches R package)
    if theta is None:
        theta = np.repeat([2] * len(phi_n), phi_n).astype(np.float32)
    elif isinstance(theta, (float, int)):
        theta = np.repeat([theta] * len(phi_n), phi_n).astype(np.float32)
    elif len(theta) == len(phi_n):
        theta = np.repeat([theta], phi_n).astype(np.float32)
    else:
        theta = np.asarray(theta, dtype=np.float32)

    assert len(theta) == np.sum(phi_n), \
        "each batch variable must have a theta"

    # Lambda handling (matches R harmony2: NULL = auto-estimation)
    lambda_estimation = False
    if lamb is None or lamb == -1:
        lambda_estimation = True
        lamb = np.zeros(1, dtype=np.float32)
    elif isinstance(lamb, (float, int)):
        lamb = np.repeat([lamb] * len(phi_n), phi_n).astype(np.float32)
        lamb = np.insert(lamb, 0, 0).astype(np.float32)
    elif len(lamb) == len(phi_n):
        lamb = np.repeat([lamb], phi_n).astype(np.float32)
        lamb = np.insert(lamb, 0, 0).astype(np.float32)
    else:
        lamb = np.asarray(lamb, dtype=np.float32)
        if len(lamb) == np.sum(phi_n):
            lamb = np.insert(lamb, 0, 0).astype(np.float32)

    # Number of items in each category
    B = int(np.sum(phi_n))
    N_b = np.bincount(batch_of_cell.ravel(), minlength=B).astype(np.float32)
    Pr_b = (N_b / N).astype(np.float32)

    if tau > 0:
        theta = theta * (1 - np.exp(-(N_b / (nclust * tau)) ** 2))

    if verbose:
        logger.info("Running Harmony")
        logger.info("  Parameters:")
        logger.info(f"    max_iter_harmony: {max_iter_harmony}")
        logger.info(f"    max_iter_kmeans: {max_iter_kmeans}")
        logger.info(f"    epsilon_cluster: {epsilon_cluster}")
        logger.info(f"    epsilon_harmony: {epsilon_harmony}")
        logger.info(f"    nclust: {nclust}")
        logger.info(f"    block_size: {block_size}")
        if lambda_estimation:
            logger.info(f"    lamb: dynamic (alpha={alpha})")
        else:
            logger.info(f"    lamb: {lamb[1:]}")
        logger.info(f"    theta: {theta}")
        logger.info(f"    sigma: {sigma[:5]}..." if len(sigma) > 5 else f"    sigma: {sigma}")
        logger.info(f"    verbose: {verbose}")
        logger.info(f"    random_state: {random_state}")
        logger.info(f"  Data: {data_mat.shape[0]} PCs × {N} cells")
        logger.info(f"  Batch variables: {vars_use}")

    # Prepare arrays for C++ backend
    data_f64 = np.ascontiguousarray(data_mat.astype(np.float64))
    batch_of_cell_c = np.ascontiguousarray(batch_of_cell)

    # Signal lambda estimation with sentinel [-1]
    if lambda_estimation:
        lamb_cpp = np.array([-1.0], dtype=np.float64)
    else:
        lamb_cpp = lamb.astype(np.float64)

    # Set BLAS thread count (Accelerate on macOS, OpenBLAS on Linux).
    # ncores=0 means use all available cores (don't set env vars).
    if ncores > 0:
        os.environ["OMP_NUM_THREADS"] = str(ncores)
        os.environ["OPENBLAS_NUM_THREADS"] = str(ncores)

    cpp_harmony = HarmonyCpp(
        data_f64,
        batch_of_cell_c,
        Pr_b.astype(np.float64),
        sigma.astype(np.float64),
        theta.astype(np.float64),
        lamb_cpp,
        float(alpha),
        max_iter_harmony,
        max_iter_kmeans,
        float(epsilon_cluster),
        float(epsilon_harmony),
        nclust,
        float(block_size),
        phi_n.tolist(),
        float(batch_prop_cutoff),
        verbose,
        random_state if random_state is not None else 0,
        logger.info,
    )
    return Harmony(cpp_harmony)
