#!/usr/bin/env python3
"""Check a fresh native install using a durable independent consumer source."""
import argparse, pathlib, subprocess, tempfile, shutil
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--profile",choices=["debug","release"],default="debug")
profile=parser.parse_args().profile
root=pathlib.Path(__file__).resolve().parents[1]
native=root/("build/native" if profile=="debug" else "build/native-release")
cmake=str(root/'.build-tools/bin/cmake')
with tempfile.TemporaryDirectory(prefix='irred-install-check-') as scratch:
 prefix=pathlib.Path(scratch)/'prefix'
 subprocess.run([cmake,'--install',str(native),'--prefix',str(prefix)],check=True)
 if (prefix/'include/cosmology').exists() or not (prefix/'include/irred/numerics.hpp').exists():
  raise SystemExit('incorrect installed public include root')
 for source_name in ['test_installed_consumer.cpp', 'test_installed_continuous_cmb_projection.cpp', 'test_helium_escape_history_installed.cpp']:
  output=pathlib.Path(scratch)/source_name.removesuffix('.cpp')
  subprocess.run([shutil.which('c++'),'-std=c++20','-Wall','-Wextra','-Wpedantic','-fno-fast-math','-ffp-contract=off',str(root/'cpp/tests'/source_name),'-I',str(prefix/'include'),str(prefix/'lib/libirred_core.a'),'-o',str(output)],check=True)
  subprocess.run([str(output)],check=True)
 print('Fresh irred install and standalone consumer passed')
