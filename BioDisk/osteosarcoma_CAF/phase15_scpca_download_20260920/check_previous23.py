import json, zipfile
from range_zip import RemoteZip, ROOT

samples = ['SCPCS000785', 'SCPCS000786', 'SCPCS000789', 'SCPCS000795']
remote = RemoteZip('SCPCP000023', block_size=65536)
with zipfile.ZipFile(remote) as archive:
    rows = []
    for sample in samples:
        entries = [i for i in archive.infolist() if sample in i.filename and not i.is_dir() and '_bulk/' not in i.filename]
        matrices = [i for i in entries if i.filename.endswith('.h5ad')]
        row = {'sample': sample, 'files': len(entries), 'matrix_files': [i.filename for i in matrices], 'bytes': sum(i.file_size for i in entries), 'downloaded_files': sum((ROOT / 'files' / 'SCPCP000023' / i.filename).is_file() for i in entries)}
        rows.append(row)
    print(json.dumps(rows, indent=2))
    (ROOT / 'previous23_availability.json').write_text(json.dumps(rows, indent=2))
