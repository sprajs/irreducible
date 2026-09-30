#!/usr/bin/env python3
"""Check a fresh native install using a durable independent consumer source."""
import pathlib, subprocess, tempfile, shutil
root=pathlib.Path(__file__).resolve().parents[1]
cmake=str(root/'.build-tools/bin/cmake')
with tempfile.TemporaryDirectory(prefix='irred-install-check-') as scratch:
 prefix=pathlib.Path(scratch)/'prefix'
 subprocess.run([cmake,'--install',str(root/'build/native'),'--prefix',str(prefix)],check=True)
 if (prefix/'include/cosmology').exists() or not (prefix/'include/irred/numerics.hpp').exists():
  raise SystemExit('incorrect installed public include root')
 output=pathlib.Path(scratch)/'consumer'
 subprocess.run([shutil.which('c++'),'-std=c++20','-Wall','-Wextra','-Wpedantic',str(root/'cpp/tests/test_installed_consumer.cpp'),'-I',str(prefix/'include'),str(prefix/'lib/libirred_core.a'),'-o',str(output)],check=True)
 subprocess.run([str(output)],check=True)
 print('Fresh irred install and standalone consumer passed')
