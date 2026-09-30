#!/usr/bin/env python3
"""One entry point: generate, CMake, then Cargo. No recursive build managers."""
import argparse,hashlib,json,pathlib,subprocess,os,shutil
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--profile",choices=["debug","release"],default="debug")
profile=parser.parse_args().profile
native_dir="build/native" if profile=="debug" else "build/native-release"
manifest_file="build/build-manifest.json" if profile=="debug" else "build/build-manifest-release.json"
cmake_profile="Debug" if profile=="debug" else "Release"
profile_flags="-g" if profile=="debug" else "-O3 -DNDEBUG"
root=pathlib.Path(__file__).resolve().parents[1]; os.chdir(root)
def run(*args): subprocess.run(args,check=True)
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
# Qualified bootstrap profile rejects invisible environment flag/target overrides.
for name in ['RUSTFLAGS','CARGO_ENCODED_RUSTFLAGS','CARGO_BUILD_TARGET','CARGO_BUILD_RUSTFLAGS','CXXFLAGS','CPPFLAGS','LDFLAGS','CC','CXX','RUSTC','RUSTC_WRAPPER','RUSTC_WORKSPACE_WRAPPER','CARGO_PROFILE_DEV_OPT_LEVEL','CARGO_PROFILE_DEV_PANIC']:
 if os.environ.get(name): raise SystemExit(f'unsupported build override {name}; use the declared profile')
for name in os.environ:
 if name.startswith(('CARGO_TARGET_','CARGO_PROFILE_')): raise SystemExit(f'unsupported build override {name}')
for directory in [root,*root.parents,pathlib.Path.home()/'.cargo']:
 for config in ['.cargo/config','.cargo/config.toml'] if directory!=pathlib.Path.home()/'.cargo' else ['config','config.toml']:
  if (directory/config).exists(): raise SystemExit(f'unsupported Cargo config override: {directory/config}')
run('python3','tools/generate_abi.py')
if not (root/'Cargo.lock').exists(): run('cargo','generate-lockfile','--offline')
cmake=str(root/'.build-tools/bin/cmake'); compiler=str(pathlib.Path(shutil.which('c++')).resolve())
run(cmake,'-S','cpp','-B',native_dir,f'-DCMAKE_BUILD_TYPE={cmake_profile}',f'-DCMAKE_CXX_COMPILER={compiler}','-DCMAKE_CXX_FLAGS=',f'-DCMAKE_CXX_FLAGS_{cmake_profile.upper()}={profile_flags}','-DCMAKE_EXE_LINKER_FLAGS=','-DCMAKE_STATIC_LINKER_FLAGS=','-DCOSMOLOGY_TEST_CFITSIO=OFF','-DIRRED_REFERENCE_GMP_MPFR=OFF')
paths=sorted([*root.glob('src/**/*.rs'),*root.glob('tests/**/*.rs'),*root.glob('cpp/**/*.cpp'),*root.glob('cpp/**/*.h'),*root.glob('cpp/**/*.hpp'),*root.glob('cpp/**/*.inc'),*root.glob('schema/*.json'),*root.glob('tools/*.py')]+[root/'Cargo.toml',root/'Cargo.lock',root/'build.rs',root/'cpp/CMakeLists.txt'])
manifest={'schema_version':1,'sources':{str(p.relative_to(root)):digest(p) for p in paths},'tools':{t:subprocess.check_output([t,'--version'],text=True).splitlines()[0] for t in ['rustc','cargo',compiler,cmake]},'compiler_executable_digest':digest(pathlib.Path(compiler)),'profile':profile,'flags':[f'CMAKE_BUILD_TYPE={cmake_profile}',*profile_flags.split(),'-std=c++20','-Wall','-Wextra','-Wpedantic','-fno-fast-math','-ffp-contract=off'],'native_compile_commands':json.loads((root/native_dir/'compile_commands.json').read_text()),'rust_flags':(['dev','opt-level=0','debuginfo=2','panic=unwind'] if profile=='debug' else ['release','opt-level=3','debuginfo=0','panic=unwind','lto=false','codegen-units=16','overflow-checks=false','incremental=false']),'target':subprocess.check_output(['rustc','-vV'],text=True),'linker':subprocess.check_output(['ld','--version'],text=True).splitlines()[0],'standard_library':subprocess.check_output([compiler,'-print-file-name=libstdc++.so'],text=True).strip(),'standard_library_digest':digest(pathlib.Path(subprocess.check_output([compiler,'-print-file-name=libstdc++.so'],text=True).strip()).resolve()),'tool_executable_digests':{t:digest(pathlib.Path(shutil.which(t)).resolve()) for t in ['rustc','cargo']},'jobs':4,'backend':'portable_cpu','panic':'unwind','git_head':subprocess.run(['git','rev-parse','HEAD'],capture_output=True,text=True).stdout.strip() or None}
revision=manifest.pop('git_head')
manifest['build_id']=hashlib.sha256(json.dumps(manifest,sort_keys=True,separators=(',',':')).encode()).hexdigest()
manifest['git_head']=revision
manifest['git_status']=subprocess.check_output(['git','status','--porcelain'],text=True)
(root/'build').mkdir(exist_ok=True); (root/manifest_file).write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
run(cmake,'--build',native_dir,'--parallel','4')
run('cargo','build','--locked','--offline','-j','4',*(['--release'] if profile=='release' else []))
