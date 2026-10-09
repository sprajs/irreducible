#!/usr/bin/env python3
"""Check dispatch against registered CMake tests and discovered Cargo targets."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from ci_suites import NATIVE_FAST, CLI_FAST, native_build_targets, classify_native, native_regex, cli_targets
import re
import sys
import runpy
from unittest.mock import patch
import ci

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

    def test_actual_install_dispatch_keeps_full_consumers_opt_in(self):
        for suite in ('fast', 'full'):
            with tempfile.TemporaryDirectory() as scratch:
                def execute(command, **kwargs):
                    if '--install' in command:
                        prefix = Path(command[command.index('--prefix') + 1])
                        (prefix / 'include/irred').mkdir(parents=True)
                        (prefix / 'include/irred/numerics.hpp').touch()
                with patch.object(sys, 'argv', ['check_install.py', '--suite', suite]), \
                     patch('tempfile.TemporaryDirectory') as directory, \
                     patch('subprocess.run', side_effect=execute) as run, \
                     patch('builtins.print'):
                    directory.return_value.__enter__.return_value = scratch
                    runpy.run_path(str(ROOT / 'tools/check_install.py'), run_name='__main__')
                compilations = [call.args[0] for call in run.call_args_list if '-std=c++20' in call.args[0]]
                self.assertEqual(len(compilations), 1 if suite == 'fast' else 13)
                names = [next(str(x) for x in command if str(x).endswith('.cpp')) for command in compilations]
                self.assertEqual(any(n.endswith('test_installed_ideal_acoustic.cpp') for n in names), suite == 'full')
                for model in ('quintessence.cpp', 'test_dgp_growth_installed.cpp', 'test_nfw_halo_installed.cpp', 'decaying_matter.cpp', 'test_curved_flrw_sdk.cpp', 'test_hernquist_sphere_installed.cpp'):
                    self.assertEqual(any(n.endswith('/' + model) for n in names), suite == 'full')
                for command in compilations:
                    self.assertIn(str(Path(scratch) / 'prefix/include'), command)
                    self.assertIn(str(Path(scratch) / 'prefix/lib/libirred_core.a'), command)
                    self.assertNotIn(str(ROOT / 'cpp/include'), command)

    def test_actual_native_dispatch_selects_targets_and_ctest_filter(self):
        inventory = {'tests': [{'name': n} for n in (*NATIVE_FAST, 'ideal_acoustic_contract')]}
        for suite in ('fast', 'full'):
            with patch.object(sys, 'argv', ['ci.py', 'native', '--suite', suite]), \
                 patch.object(ci.shutil, 'which', return_value='/usr/bin/g++'), \
                 patch.object(ci.subprocess, 'check_output', return_value=json.dumps(inventory)), \
                 patch.object(ci, 'run') as run, patch('builtins.print'):
                ci.main()
            calls = [call.args for call in run.call_args_list]
            build = next(c for c in calls if '--build' in c)
            ctest = calls[-1]
            self.assertIn('--no-tests=error', ctest)
            if suite == 'fast':
                self.assertIn('--target', build)
                self.assertIn('test_abi', build)
                self.assertNotIn('test_ideal_acoustic', build)
                self.assertIn('-R', ctest)
                self.assertIsNone(re.fullmatch(ctest[-1], 'ideal_acoustic_contract'))
            else:
                self.assertNotIn('--target', build)
                self.assertNotIn('-R', ctest)

    def test_actual_native_dispatch_preserves_requested_profile_and_isolates_builds(self):
        inventory = {'tests': [{'name': n} for n in (*NATIVE_FAST, 'ideal_acoustic_contract')]}
        for compiler in ('gcc', 'clang'):
            for profile, build_type, suffix in (('debug', 'Debug', ''), ('release', 'Release', '-release')):
                with self.subTest(compiler=compiler, profile=profile):
                    with patch.object(sys, 'argv', ['ci.py', 'native', '--compiler', compiler,
                                                    '--profile', profile, '--suite', 'fast']), \
                         patch.object(ci.shutil, 'which', return_value=f'/usr/bin/{compiler}++'), \
                         patch.object(ci.subprocess, 'check_output', return_value=json.dumps(inventory)) as inventory_read, \
                         patch.object(ci, 'run') as run, patch('builtins.print'):
                        ci.main()
                    calls = [call.args for call in run.call_args_list]
                    configure = next(c for c in calls if '-S' in c)
                    build = next(c for c in calls if '--build' in c)
                    ctest = calls[-1]
                    directory = f'build/ci-native-{compiler}{suffix}'
                    self.assertIn(f'-DCMAKE_BUILD_TYPE={build_type}', configure)
                    self.assertEqual(configure[configure.index('-B') + 1], directory)
                    self.assertEqual(build[build.index('--build') + 1], directory)
                    self.assertEqual(ctest[ctest.index('--test-dir') + 1], directory)
                    query = inventory_read.call_args.args[0]
                    self.assertEqual(query[query.index('--test-dir') + 1], directory)
                    self.assertIn('test_abi', build)
                    self.assertNotIn('test_ideal_acoustic', build)
                    self.assertIn('-R', ctest)

    def test_actual_cli_dispatch_preserves_full_and_release(self):
        metadata = [{'name': n, 'kind': ['test']} for n in (*CLI_FAST, 'sound_horizon_cli')]
        for suite in ('fast', 'full'):
            with patch.object(sys, 'argv', ['ci.py', 'cli', '--suite', suite, '--profile', 'release']), \
                 patch.object(ci, 'targets', return_value=metadata), \
                 patch.object(ci, 'run') as run, patch('builtins.print'):
                ci.main()
            command = run.call_args.args
            self.assertIn('--release', command)
            self.assertIn('--locked', command)
            self.assertIn('--offline', command)
            self.assertEqual('sound_horizon_cli' in command, suite == 'full')
            self.assertIn('cli', command)

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
