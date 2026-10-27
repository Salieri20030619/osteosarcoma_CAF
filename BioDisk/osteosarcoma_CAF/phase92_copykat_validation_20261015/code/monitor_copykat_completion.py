from pathlib import Path
import ctypes, datetime, fcntl, hashlib, json, os, select, struct, subprocess

O = Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase92_copykat_validation_20261015')
PY = '/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python'
SAMPLES = ['GSE234187__GSM7453693'] + ['NCC__' + s for s in ['BS7','BS1','BS12','BS14','BS9','BS5','BS4','BS3','BS2']]
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(path, value):
    tmp = path.with_suffix('.monitor.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')
    tmp.replace(path)

lock = (O/'completion_monitor.lock').open('a')
fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
libc = ctypes.CDLL(None, use_errno=True)
fd = libc.inotify_init1(os.O_NONBLOCK | os.O_CLOEXEC)
if fd < 0: raise OSError(ctypes.get_errno(), 'inotify_init1')
for directory in [O, O/'resources']:
    if libc.inotify_add_watch(fd, os.fsencode(directory), 0x8 | 0x80) < 0:
        raise OSError(ctypes.get_errno(), 'inotify_add_watch')
pids = [1040168, json.loads((O/'resources/SECONDARY_LANE_RECOVERY.json').read_text(encoding='utf-8'))['pid']]
watch = {}
for pid in pids:
    try: watch[os.pidfd_open(pid)] = pid
    except ProcessLookupError: pass
observed = {}
checks = {}
def scan():
    for sm in SAMPLES:
        receipt = O/'resources'/(sm+'_COMPLETE.json')
        if not receipt.exists(): continue
        digest = hashlib.sha256(receipt.read_bytes()).hexdigest()
        if observed.get(sm) == digest: continue
        with (O/'logs'/(sm+'_validity_recheck.log')).open('w', encoding='utf-8') as log:
            result = subprocess.run([PY, str(O/'code/validate_method_controls.py'), sm], stdout=log, stderr=subprocess.STDOUT, timeout=300)
        entry = {'utc':now(), 'validator_exit_code':result.returncode, 'receipt_sha256':digest}
        report = O/'resources'/(sm+'_VALIDITY_RECHECK.json')
        if result.returncode == 0 and report.exists():
            v = json.loads(report.read_text(encoding='utf-8'))
            entry.update(copykat_reference_pass=v['CopyKAT_valid_for_binary_identity'], full_cells=v['full_cells'], patient_id=v['patient_id'])
        observed[sm] = digest
        checks[sm] = entry
        print('VALIDATION_EVENT', sm, json.dumps(entry), flush=True)
    dump(O/'MONITOR_STATUS.json', {'utc':now(), 'monitor_pid':os.getpid(), 'waiting_on_pids':list(watch.values()), 'samples':checks, 'scientific_admission_pending':True, 'new_formal_patients':0})

scan()
while watch:
    ready, _, _ = select.select([fd] + list(watch), [], [])
    relevant = False
    for w in ready:
        if w != fd:
            print('RUNNER_EXIT_EVENT', watch.pop(w), now(), flush=True)
            os.close(w)
            relevant = True
        else:
            while True:
                try: data = os.read(fd, 65536)
                except BlockingIOError: break
                offset = 0
                while offset < len(data):
                    wd, mask, cookie, length = struct.unpack_from('iIII', data, offset)
                    offset += 16
                    name = data[offset:offset+length].rstrip(b'\0').decode('utf-8')
                    offset += length
                    if name.endswith('_COMPLETE.json') or name in ['SAMPLE_PROGRESS.json','RECOVERY_LANE_PROGRESS.json','RECOVERY_LANE_FINISHED.json']:
                        relevant = True
    if relevant: scan()
scan()
accounting = {}
for name in ['SAMPLE_PROGRESS.json','RECOVERY_LANE_PROGRESS.json']:
    p = O/name
    if p.exists():
        for row in json.loads(p.read_text(encoding='utf-8'))['samples']:
            assert row['sample'] not in accounting
            accounting[row['sample']] = row
missing = sorted(set(SAMPLES) - set(accounting))
ok = not missing and all(v.get('status','').startswith(('complete_','ineligible_')) for v in accounting.values())
dump(O/'COMPUTATION_RECONCILED.json', {'utc':now(), 'computation_exit_code':0 if ok else 1, 'expected_libraries':10, 'missing_samples':missing, 'sample_accounting':accounting, 'method_validity_checks':checks, 'scientific_admission_pending':True, 'new_formal_patients':0, 'note':'Original secondary future failed on unavailable pidfd_open; recovery lane processed its samples. This record reconciles per-sample exits and completion markers, not supervisor disappearance.'})
print('ALL_LANES_RECONCILED', json.dumps({'computation_ok':ok,'missing':missing,'checks':checks}), flush=True)
os.close(fd)
