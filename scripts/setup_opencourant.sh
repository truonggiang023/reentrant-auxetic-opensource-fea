#!/usr/bin/env bash
# Install the official OpenCourant release and exercise both solver executables.
# TLS and the GitHub-published artifact digest are required, including cached files.
set -euo pipefail

root=${OPENCOURANT_INSTALL_ROOT:-/workspace/opencourant}
release=latest-20261011
digest=6cb419677ca36d73acb96218cc0cf5e4f7d50494d9c7812ba27935e3c0f8ec94
url="https://github.com/OpenCourant/OpenCourant/releases/download/$release/OpenCourant_linux64.zip"
mkdir -p "$root"
archive="$root/OpenCourant_linux64.zip"
if [[ ! -f "$archive" ]]; then
    curl --fail --location --retry 2 "$url" --output "$archive.part"
    mv "$archive.part" "$archive"
fi
printf '%s  %s\n' "$digest" "$archive" | sha256sum --check -

python - "$root" <<'PY'
import hashlib
import os
from pathlib import Path
import shlex
import sys
import zipfile

root = Path(sys.argv[1]).resolve()
runtime = root / 'runtime'
with zipfile.ZipFile(root / 'OpenCourant_linux64.zip') as archive:
    for entry in archive.infolist():
        destination = (runtime / entry.filename).resolve()
        if not destination.is_relative_to(runtime.resolve()):
            raise SystemExit('Archive entry escapes runtime directory')
        if entry.is_dir():
            destination.mkdir(parents=True, exist_ok=True)
            continue
        contents = archive.read(entry)
        if destination.exists():
            if hashlib.sha256(destination.read_bytes()).digest() != hashlib.sha256(contents).digest():
                raise SystemExit(f'Refusing to overwrite changed runtime file: {destination}')
        else:
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(contents)
        if destination.parent.name == 'exec':
            destination.chmod(0o755)

installation = runtime / 'OpenCourant'
initialization = '\n'.join([
    'export OPENCOURANT_PATH=' + shlex.quote(str(installation)),
    'export RAD_CFG_PATH="$OPENCOURANT_PATH/hm_cfg_files"',
    'export RAD_H3D_PATH="$OPENCOURANT_PATH/extlib/h3d/lib/linux64"',
    'export OMP_STACKSIZE=400m',
    'export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}',
    'export LD_LIBRARY_PATH="$OPENCOURANT_PATH/extlib/hm_reader/linux64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"',
    'export PATH="$OPENCOURANT_PATH/exec:$PATH"',
    '',
])
(root / 'env.sh').write_text(initialization)
PY

# shellcheck disable=SC1090
source "$root/env.sh"
smoke=$(mktemp -d "$root/smoke.XXXXXX")
source_base="https://raw.githubusercontent.com/OpenCourant/OpenCourant/$release/qa-tests/miniqa/SMOKE_TEST/data"
for name in TWISBEAM_0000.rad TWISBEAM_0001.rad TWISBEAMD0B ref.extract; do
    curl --fail --location --retry 2 "$source_base/$name" --output "$smoke/$name"
done
cd "$smoke"
starter_linux64_gf -i TWISBEAM_0000.rad -np 1 > starter.log 2>&1
engine_linux64_gf -i TWISBEAM_0001.rad > engine.log 2>&1

python - "$root" "$release" "$digest" <<'PY'
import hashlib
import json
from pathlib import Path
import re
import sys

root, release, digest = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
starter, engine = Path('starter.log').read_text(), Path('engine.log').read_text()
starter_errors = re.search(r'(\d+)\s+ERROR\(S\)', starter)
starter_warnings = re.findall(r'WARNING ID\s*:\s*(\d+)', starter)
if (not starter_errors or int(starter_errors.group(1)) != 0
        or any(code != '100214' for code in starter_warnings)
        or 'NORMAL TERMINATION' not in engine):
    raise SystemExit('Starter input check or engine termination failed')
if not Path('TWISBEAMT01').is_file() or not Path('TWISBEAMT01').stat().st_size:
    raise SystemExit('Time-history output is missing')
reference = {}
for line in Path('ref.extract').read_text().splitlines():
    fields = line.replace('\\t', ' ').split()
    if len(fields) >= 3:
        reference[fields[1]] = float(fields[2])
cycle = int(reference['CYCLE'])
rows = [line.split() for line in Path('TWISBEAM_0001.out').read_text().splitlines()
        if re.match(r'^\s*' + str(cycle) + r'\s', line)]
if len(rows) != 1:
    raise SystemExit('Reference cycle did not execute exactly once')
row = rows[0]
columns = {'CYCLE': 0, 'TIME': 1, 'TSTEP': 2, 'ELTID': 4, 'ERROR': 5,
           'IENERGY': 6, 'KENERGYT': 7, 'KENERGYR': 8, 'EXTWORK': 9, 'MASERR': 10}
checks = {}
for name, column in columns.items():
    actual = float(row[column].rstrip('%'))
    # Both reference and output use the same printed precision.
    if abs(actual - reference[name]) > max(1e-12, abs(reference[name]) * 1e-8):
        raise SystemExit(f'Smoke reference mismatch for {name}: {actual} != {reference[name]}')
    checks[name] = {'actual': actual, 'reference': reference[name], 'passed': True}
commit = re.search(r'CommitID:\s*([0-9a-f]{40})', starter)
if not commit:
    raise SystemExit('Missing executable build identity')
record = {
    'backend': 'OpenCourant', 'release': release,
    'archive_url': f'https://github.com/OpenCourant/OpenCourant/releases/download/{release}/OpenCourant_linux64.zip',
    'archive_sha256': digest,
    'digest_source': f'https://github.com/OpenCourant/OpenCourant/releases/expanded_assets/{release}',
    'binary_commit': commit.group(1),
    'installed_root': str(root / 'runtime' / 'OpenCourant'),
    'smoke_directory': str(Path.cwd()),
    'starter_status': 'completed_with_warning' if starter_warnings else 'completed',
    'starter_errors': 0, 'starter_warning_ids': starter_warnings,
    'starter_warning_explanation': 'Official historical smoke input has unsupported /IOFLAG fields (100214); input check reports zero errors.',
    'engine_normal_termination': True,
    'smoke_checks': checks,
    'completed_cycles': int(re.search(r'TOTAL NUMBER OF CYCLES\s*:\s*(\d+)', engine).group(1)),
    'input_sha256': {name: hashlib.sha256(Path(name).read_bytes()).hexdigest()
                     for name in ('TWISBEAM_0000.rad', 'TWISBEAM_0001.rad', 'TWISBEAMD0B', 'ref.extract')},
    'qualification': 'Runtime and official shell smoke regression only; not honeycomb fracture/contact validation.',
}
(root / 'runtime-validation.json').write_text(json.dumps(record, indent=2) + '\n')
print(f'OpenCourant {release}: starter and engine completed; {len(checks)} reference checks passed.')
print('Runtime evidence: ' + str(root / 'runtime-validation.json'))
PY
