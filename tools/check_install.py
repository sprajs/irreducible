#!/usr/bin/env python3
"""Check a fresh native install using a durable independent consumer source."""
import argparse, pathlib, subprocess, tempfile, shutil
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--profile",choices=["debug","release"],default="debug")
parser.add_argument("--suite",choices=["fast","full"],default="fast",help="fast installed ABI/ownership smoke or all scientific consumers")
arguments=parser.parse_args()
profile=arguments.profile
root=pathlib.Path(__file__).resolve().parents[1]
native=root/("build/native" if profile=="debug" else "build/native-release")
cmake=str(root/'.build-tools/bin/cmake')
with tempfile.TemporaryDirectory(prefix='irred-install-check-') as scratch:
 prefix=pathlib.Path(scratch)/'prefix'
 subprocess.run([cmake,'--install',str(native),'--prefix',str(prefix)],check=True)
 if (prefix/'include/cosmology').exists() or not (prefix/'include/irred/numerics.hpp').exists():
  raise SystemExit('incorrect installed public include root')
 sources=[root/'cpp/tests/test_installed_consumer.cpp', root/'cpp/tests/test_installed_continuous_cmb_projection.cpp', root/'cpp/examples/windowed_linear_power.cpp', root/'cpp/tests/test_effective_fluid_installed.cpp', root/'cpp/tests/test_two_deflector_forward_installed.cpp', root/'cpp/tests/test_installed_finite_opacity_source.cpp', root/'cpp/tests/test_installed_ideal_acoustic.cpp']
 if arguments.suite=='fast':
  sources=[root/'cpp/tests/test_abi.cpp']
 for source in sources:
  output=pathlib.Path(scratch)/source.stem
  subprocess.run([shutil.which('c++'),'-std=c++20','-Wall','-Wextra','-Wpedantic','-fno-fast-math','-ffp-contract=off',str(source),'-I',str(prefix/'include'),str(prefix/'lib/libirred_core.a'),'-o',str(output)],check=True)
  subprocess.run([str(output)],check=True)
 print(f'Fresh irred install and standalone {arguments.suite} consumers passed')
