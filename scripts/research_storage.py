#!/usr/bin/env python3
"""Named shared research storage. Exact versions, verified bytes, no physics."""
import argparse
import base64
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / 'storage-layout.json'
SHA = re.compile(r'[0-9a-f]{64}')
BLOCKED = {'.git', '.aws', '.ssh', '.venv', 'venv', 'node_modules', '__pycache__'}
SECRET = re.compile(rb'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----|(?:AKIA|ASIA)[A-Z0-9]{16}|gh[pousr]_[A-Za-z0-9]{20,}|sk-(?:proj-)?[A-Za-z0-9_-]{24,}|X-Amz-(?:Signature|Credential)=|(?:aws_secret_access_key|secret_access_key|access_token|refresh_token|password)\s*[=:]\s*["\x27]?[A-Za-z0-9_+/=-]{12,}', re.I)


def safe(name):
    if (not isinstance(name, str) or not name or name.startswith('/') or '\\' in name
            or any(ord(c) < 32 or ord(c) == 127 for c in name)
            or any(p in {'', '.', '..'} for p in name.split('/'))):
        raise ValueError('Unsafe object path')
    return name


def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def strict_json(data):
    def pairs(items):
        result = {}
        for k, v in items:
            if k in result:
                raise ValueError('Duplicate JSON key')
            result[k] = v
        return result
    return json.loads(data, object_pairs_hook=pairs,
                      parse_constant=lambda _: (_ for _ in ()).throw(ValueError('Nonfinite JSON')))


def identity(row):
    if (type(row.get('bytes')) is not int or row['bytes'] < 0
            or not isinstance(row.get('sha256'), str) or not SHA.fullmatch(row['sha256'])):
        raise ValueError('Invalid byte identity')


def inventory(root):
    if root.is_symlink() or not root.is_dir():
        raise ValueError('Selection must be a real directory')
    rows = []
    for p in sorted(root.rglob('*')):
        relative = safe(p.relative_to(root).as_posix())
        mode = p.lstat().st_mode
        if stat.S_ISLNK(mode) or not (stat.S_ISDIR(mode) or stat.S_ISREG(mode)):
            raise ValueError('Links/special files refused')
        if BLOCKED.intersection(Path(relative).parts) or p.name == '.env' or p.name.startswith('.env.'):
            raise ValueError('Credential/environment directories refused')
        if stat.S_ISREG(mode):
            before = p.stat()
            h = hashlib.sha256()
            tail = b''
            with p.open('rb') as f:
                while block := f.read(1024 * 1024):
                    if SECRET.search(tail + block):
                        raise ValueError('Potential credential marker; inspect locally')
                    tail = block[-1024:]
                    h.update(block)
            after = p.stat()
            if (before.st_size, before.st_mtime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino):
                raise ValueError('Selection changed while hashing')
            rows.append({'path': relative, 'bytes': after.st_size, 'sha256': h.hexdigest()})
    if not rows:
        raise ValueError('Empty selection')
    return rows


class Storage:
    def __init__(self, profile=None, ambient=False):
        self.config = strict_json(CONFIG.read_text())
        self.profile = None if ambient else profile or os.environ.get('AWS_PROFILE')
        if not ambient and not self.profile and not any(os.environ.get(k) for k in (
                'AWS_ACCESS_KEY_ID', 'AWS_WEB_IDENTITY_TOKEN_FILE', 'AWS_CONTAINER_CREDENTIALS_RELATIVE_URI',
                'AWS_CONTAINER_CREDENTIALS_FULL_URI')):
            self.profile = self.config['local_profile']

    def aws(self, *args):
        cmd = ['aws', *args, '--region', self.config['region'], '--output', 'json', '--no-cli-pager']
        if self.profile:
            cmd += ['--profile', self.profile]
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if p.returncode:
            raise RuntimeError('Restricted AWS operation failed')
        return strict_json(p.stdout)

    def check(self):
        i = self.aws('sts', 'get-caller-identity')
        if i['Account'] != self.config['account_id'] or i['Arn'].endswith(':root'):
            raise ValueError('Restricted configured research identity required')
        if self.aws('s3api', 'get-bucket-versioning', '--bucket', self.config['bucket']).get('Status') != 'Enabled':
            raise ValueError('Versioning required')
        policy = self.aws('s3api', 'get-bucket-lifecycle-configuration', '--bucket', self.config['bucket'])
        if any(r.get('Status') == 'Enabled' and ('Expiration' in r or 'NoncurrentVersionExpiration' in r) for r in policy.get('Rules', [])):
            raise ValueError('Completed-object expiry requires review')
        return i

    def key(self, uri):
        prefix = 's3://' + self.config['bucket'] + '/'
        if not uri.startswith(prefix):
            raise ValueError('Unconfigured bucket')
        return safe(uri[len(prefix):])

    def get(self, pin, path):
        key = self.key(pin['uri'])
        if not pin.get('version_id') or pin['version_id'] == 'null' or not SHA.fullmatch(pin['sha256']):
            raise ValueError('Exact manifest/object version and SHA256 required')
        result = self.aws('s3api', 'get-object', '--bucket', self.config['bucket'], '--key', key,
                          '--version-id', pin['version_id'], '--checksum-mode', 'ENABLED', str(path))
        if result.get('VersionId') != pin['version_id'] or digest(path) != pin['sha256']:
            raise ValueError('Exact-version bytes differ')
        if pin.get('bytes') is not None and path.stat().st_size != pin['bytes']:
            raise ValueError('Exact-version size differs')

    def put(self, key, path, sha, size):
        safe(key)
        if size >= 5 * 1024**3 or digest(path) != sha or path.stat().st_size != size:
            raise ValueError('Changed bytes or single PUT limit')
        try:
            r = self.aws('s3api', 'put-object', '--bucket', self.config['bucket'], '--key', key,
                         '--body', str(path), '--if-none-match', '*', '--checksum-sha256',
                         base64.b64encode(bytes.fromhex(sha)).decode(), '--server-side-encryption', 'AES256',
                         '--metadata', json.dumps({'sha256': sha}))
        except RuntimeError:
            r = self.aws('s3api', 'head-object', '--bucket', self.config['bucket'], '--key', key)
            if r.get('ContentLength') != size or r.get('Metadata', {}).get('sha256') != sha:
                raise ValueError('Occupied immutable key differs')
        pin = {'uri': 's3://' + self.config['bucket'] + '/' + key,
               'version_id': r.get('VersionId'), 'sha256': sha, 'bytes': size}
        with tempfile.TemporaryDirectory() as td:
            self.get(pin, Path(td) / 'readback')
        return pin


def collection(name):
    safe(name)
    if not re.fullmatch(r'[a-z0-9][a-z0-9._/-]*', name):
        raise ValueError('Use readable lowercase names')
    patterns = (r'shared/datasets/[^/]+/[^/]+', r'shared/archives/[^/]+',
                r'reproducible/experiments/[^/]+/attempts/[^/]+', r'prospector/papers/[^/]+',
                r'prospector/reviews/[^/]+', r'irreducible/evidence/[^/]+/[^/]+')
    if not any(re.fullmatch(p, name) for p in patterns):
        raise ValueError('Collection must follow storage-layout.json roots')
    return name


def push(storage, source, name, role, rights, provenance, dry):
    collection(name)
    rows = inventory(source)
    p = strict_json(provenance.read_text())
    if not isinstance(p, dict) or not p:
        raise ValueError('Explicit provenance object required')
    if SECRET.search(provenance.read_bytes()):
        raise ValueError('Potential credential marker in provenance')
    descriptor = {'schema': 'research-named-selection/v1', 'collection': name,
                  'scientific_role': role, 'rights': rights, 'provenance': p, 'files': rows}
    release = hashlib.sha256(json.dumps(descriptor, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    root = name + '/versions/' + release
    if dry:
        return {'dry_run': True, 'collection': name, 'immutable_prefix': root,
                'files': len(rows), 'bytes': sum(r['bytes'] for r in rows), 'network_requests': 0}
    storage.check()
    def upload(r):
        pin = storage.put(root + '/files/' + r['path'], source / r['path'], r['sha256'], r['bytes'])
        return dict(r, **{'uri': pin['uri'], 'version_id': pin['version_id']})
    with ThreadPoolExecutor(max_workers=2) as pool:
        versions = list(pool.map(upload, rows))
    if inventory(source) != rows:
        raise ValueError('Source changed; incomplete attempt retained')
    doc = dict(descriptor, schema='research-named-manifest/v1', release=release, files=versions)
    with tempfile.TemporaryDirectory() as td:
        marker = Path(td) / 'manifest.json'
        marker.write_text(json.dumps(doc, sort_keys=True, indent=2) + '\n')
        pin = storage.put(root + '/manifest.json', marker, digest(marker), marker.stat().st_size)
    records = ROOT / '.work/storage/receipts'
    records.mkdir(parents=True, exist_ok=True)
    receipt = {'manifest_format': doc['schema'], 'collection': name, **pin}
    record = records / (release + '.json')
    if record.exists() and strict_json(record.read_text()) != receipt:
        raise ValueError('Previous receipt differs')
    if not record.exists():
        with record.open('x') as f:
            json.dump(receipt, f, indent=2); f.write('\n')
    return receipt


def pull(storage, pin, destination):
    if destination.exists() or destination.is_symlink():
        raise ValueError('Restore destination must be new')
    storage.check()
    with tempfile.TemporaryDirectory() as td:
        marker = Path(td) / 'manifest.json'; storage.get(pin, marker)
        doc = strict_json(marker.read_text())
        if doc.get('schema') != 'research-named-manifest/v1':
            raise ValueError('Use catalog authoritative restore route for this format')
        collection(doc['collection'])
        rows = doc['files']; names = set()
        if not isinstance(rows, list) or not rows:
            raise ValueError('Missing file records')
        for row in rows:
            identity(row); name = safe(row['path'])
            if name in names or any(name.startswith(x + '/') or x.startswith(name + '/') for x in names):
                raise ValueError('Duplicate/conflicting paths')
            names.add(name)
            if storage.key(row['uri']) != doc['collection'] + '/versions/' + doc['release'] + '/files/' + name:
                raise ValueError('Object URI outside pinned collection')
        destination.mkdir(parents=True, exist_ok=False)
        for row in rows:
            path = destination / row['path']; path.parent.mkdir(parents=True, exist_ok=True)
            storage.get(row, path)
        if inventory(destination) != [{'path': r['path'], 'bytes': r['bytes'], 'sha256': r['sha256']} for r in rows]:
            raise ValueError('Restored inventory differs')
    return {'verified_files': len(rows), 'destination': str(destination)}


def catalog(storage):
    storage.check()
    key = storage.config['catalog_pointer']
    with tempfile.TemporaryDirectory() as td:
        pointer = Path(td) / 'pointer.json'
        # Discovery may advance; returned release pin is immutable and verified.
        storage.aws('s3api', 'get-object', '--bucket', storage.config['bucket'], '--key', key, str(pointer))
        p = strict_json(pointer.read_text())
        if p.get('schema') != 'research-catalog-pointer/v1':
            raise ValueError('Unknown discovery pointer')
        release = Path(td) / 'catalog.json'; storage.get(p['catalog'], release)
        doc = strict_json(release.read_text())
        if doc.get('schema') != 'research-catalog/v1':
            raise ValueError('Unknown catalog schema')
        return {'catalog_pin': p['catalog'], 'entries': doc['entries']}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--profile'); ap.add_argument('--ambient-credentials', action='store_true')
    commands = ap.add_subparsers(dest='command', required=True)
    commands.add_parser('status'); commands.add_parser('catalog')
    lookup = commands.add_parser('resolve'); lookup.add_argument('name')
    u = commands.add_parser('push'); u.add_argument('source', type=Path); u.add_argument('--collection', required=True)
    u.add_argument('--role', required=True, choices=['observations', 'calibrated-reductions', 'fitted-summaries', 'synthetic-controls', 'generated-results', 'source-review', 'mixed-evidence'])
    u.add_argument('--rights', required=True, choices=['owned', 'approved-redistribution', 'private-preservation-unreviewed'])
    u.add_argument('--provenance', type=Path, required=True); u.add_argument('--dry-run', action='store_true')
    d = commands.add_parser('pull'); d.add_argument('destination', type=Path); d.add_argument('--manifest-uri', required=True)
    d.add_argument('--manifest-sha256', required=True); d.add_argument('--manifest-version-id', required=True)
    a = ap.parse_args(); s = Storage(a.profile, a.ambient_credentials)
    if a.command == 'status':result = {'identity': s.check()['Arn'], **s.config}
    elif a.command == 'catalog':result = catalog(s)
    elif a.command == 'resolve':
        c = catalog(s); entries = [e for e in c['entries'] if e['name'] == a.name]
        if len(entries) != 1:raise ValueError('Unknown or ambiguous catalog name')
        result = {'catalog_pin': c['catalog_pin'], 'entry': entries[0]}
    elif a.command == 'push':result = push(s, a.source, a.collection, a.role, a.rights, a.provenance, a.dry_run)
    else:result = pull(s, {'uri': a.manifest_uri, 'sha256': a.manifest_sha256, 'version_id': a.manifest_version_id}, a.destination)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    try:main()
    except Exception as e:raise SystemExit('Storage refused/failed: ' + type(e).__name__ + '; no completion claim; partial attempts retained')
