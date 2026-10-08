"""Exercise custody failures and exact-version isolation without AWS credentials."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/research_storage.py'
if not SCRIPT.exists():
    SCRIPT = ROOT / 'tools/research_storage.py'
spec = importlib.util.spec_from_file_location('named_storage', SCRIPT)
s = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)


class Fake:
    def __init__(self):self.events = []; self.objects = {}; self.config = {'bucket': 'example'}
    def check(self):self.events.append('check')
    def put(self, key, path, sha, size):
        self.events.append(key); data = path.read_bytes()
        if s.digest(path) != sha or len(data) != size:raise ValueError('changed')
        pin = {'uri': 's3://example/' + key, 'sha256': sha, 'bytes': size, 'version_id': str(len(self.objects))}
        self.objects[(pin['uri'], pin['version_id'])] = data
        return pin
    def key(self, uri):return uri.removeprefix('s3://example/')
    def get(self, pin, path):
        path.write_bytes(self.objects[(pin['uri'], pin['version_id'])])
        if s.digest(path) != pin['sha256']:raise ValueError('corrupt')


class StorageTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name); self.source = self.root / 'selection'; self.source.mkdir()
        (self.source / 'receipt.json').write_text('{"status":"failed","role":"synthetic control"}\n')
        self.provenance = self.root / 'provenance.json'; self.provenance.write_text('{"original_attempt":"failed-01","reconstruction":"new attempt only"}')

    def test_dry_run_never_checks_credentials(self):
        f = Fake()
        r = s.push(f, self.source, 'reproducible/experiments/test/attempts/failed-01', 'synthetic-controls', 'owned', self.provenance, True)
        self.assertEqual(f.events, []); self.assertEqual(r['network_requests'], 0)

    def test_manifest_last_and_restore_uses_pinned_version(self):
        f = Fake()
        with patch.object(s, 'ROOT', self.root):
            r = s.push(f, self.source, 'reproducible/experiments/test/attempts/failed-01', 'synthetic-controls', 'owned', self.provenance, False)
        self.assertTrue(f.events[-1].endswith('/manifest.json'))
        f.objects[(r['uri'], 'new-version')] = b'corrupt replacement'
        doc = json.loads(f.objects[(r['uri'], r['version_id'])])
        row = doc['files'][0]; f.objects[(row['uri'], 'new-version')] = b'bad'
        s.pull(f, r, self.root / 'restored')
        self.assertEqual((self.root / 'restored/receipt.json').read_bytes(), (self.source / 'receipt.json').read_bytes())

    def test_failed_byte_readback_has_no_marker(self):
        f = Fake(); real = f.put
        def failed(key, *args):
            if '/files/' in key:raise ValueError('readback failed')
            return real(key, *args)
        f.put = failed
        with self.assertRaises(ValueError), patch.object(s, 'ROOT', self.root):
            s.push(f, self.source, 'shared/archives/test', 'mixed-evidence', 'owned', self.provenance, False)
        self.assertFalse(any(str(x).endswith('/manifest.json') for x in f.events))

    def test_credentials_and_symlinks_refused_before_network(self):
        (self.source / '.env').write_text('x=y')
        with self.assertRaises(ValueError):s.inventory(self.source)
        (self.source / '.env').unlink(); (self.source / 'link').symlink_to(self.provenance)
        with self.assertRaises(ValueError):s.inventory(self.source)
        (self.source / 'link').unlink()
        (self.source / 'secret').write_text('aws_secret_access_key=' + 'a' * 40)
        with self.assertRaises(ValueError):s.inventory(self.source)

    def test_paths_and_duplicate_json_refused(self):
        for name in ('../x', '/x', 'x//y', 'x/./y', 'x\\y', 'x\ny'):
            with self.assertRaises(ValueError):s.safe(name)
        with self.assertRaises(ValueError):s.strict_json('{"x":1,"x":2}')
        with self.assertRaises(ValueError):s.collection('reproducible/snapshots/x')
        with self.assertRaises(ValueError):s.collection('shared/datasets/x')

    def test_restore_corruption_and_existing_destination_refused(self):
        f = Fake()
        with patch.object(s, 'ROOT', self.root):
            r = s.push(f, self.source, 'shared/archives/test', 'mixed-evidence', 'owned', self.provenance, False)
        doc = json.loads(f.objects[(r['uri'], r['version_id'])]); row = doc['files'][0]
        f.objects[(row['uri'], row['version_id'])] = b'changed bytes'
        with self.assertRaises(ValueError):s.pull(f, r, self.root / 'corrupt')
        self.assertTrue((self.root / 'corrupt').exists())
        with self.assertRaises(ValueError):s.pull(f, r, self.source)


if __name__ == '__main__':unittest.main()
