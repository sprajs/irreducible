#!/usr/bin/env python3
"""Check dispatch against registered CMake tests and discovered Cargo targets."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from ci_suites import NATIVE_FAST, CLI_FAST, native_build_targets, classify_native, native_regex, cli_targets
import re

ROOT = Path(__file__).resolve().parents[1]


class SuitesTests(unittest.TestCase):
    def test_production_and_full_build_routes_are_distinct(self):
        self.assertEqual(native_build_targets('core'), ['irred_core'])
        self.assertEqual(native_build_targets('full'), [])
        self.assertTrue(native_build_targets('fast'))
        self.assertNotIn('test_ideal_acoustic', native_build_targets('fast'))
        self.assertNotIn('test_hydrogen_relics', native_build_targets('fast'))
        with self.assertRaises(ValueError):
            native_build_targets('typo')

    def test_selection_and_full_preservation_for_new_tests(self):
        inventory = {'tests': [{'name': n} for n in (*NATIVE_FAST, 'future_scientific_peer')]}
        result = classify_native(inventory)
        self.assertEqual(result['full_only'], ['future_scientific_peer'])
        for name in result['fast']:
            self.assertRegex(name, native_regex())
            self.assertIsNone(re.fullmatch(native_regex(), name + '_peer'))
        with self.assertRaises(ValueError):
            classify_native({'tests': []})
        with self.assertRaises(ValueError):
            classify_native({'tests': inventory['tests'] * 2})
        available = [*CLI_FAST, 'future_scientific_cli']
        self.assertNotIn('future_scientific_cli', cli_targets(available, 'fast'))
        self.assertIn('future_scientific_cli', cli_targets(available, 'full'))
        with self.assertRaises(ValueError):
            cli_targets([], 'fast')

    def test_selected_targets_are_actual_registered_commands(self):
        # Configure only: no compilation or science execution. Catch target/name
        # drift against CMake's actual graph, including included fragments.
        import shutil
        cmake = str(ROOT / '.build-tools/bin/cmake')
        if not Path(cmake).exists():
            cmake = shutil.which('cmake')
        if not cmake:
            self.skipTest('CMake not installed; live graph check runs with build tools')
        with tempfile.TemporaryDirectory() as directory:
            subprocess.run([cmake, '-S', str(ROOT / 'cpp'), '-B', directory,
                            '-DIRRED_TEST_CFITSIO=OFF', '-DIRRED_REFERENCE_GMP_MPFR=OFF'],
                           check=True, stdout=subprocess.DEVNULL)
            ctest = str(Path(cmake).with_name('ctest'))
            inventory = json.loads(subprocess.check_output(
                [ctest, '--test-dir', directory, '--show-only=json-v1'], text=True))
            classify_native(inventory)
            registered = {t['name']: t for t in inventory['tests']}
            # CTest retains commands even though test executables are not built.
            for name, target in NATIVE_FAST.items():
                command = registered[name].get('command', [])
                if command:
                    self.assertEqual(Path(command[0]).name, target)
            graph = json.loads((Path(directory) / 'compile_commands.json').read_text())
            commands = '\n'.join(row['command'] for row in graph)
            for target in native_build_targets('fast'):
                self.assertIn(f'{target}.dir', commands)
        declared = {p.stem for p in (ROOT / 'tests').glob('*.rs')}
        self.assertEqual(cli_targets(declared, 'fast'), sorted(CLI_FAST))


if __name__ == '__main__':
    unittest.main()
