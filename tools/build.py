#!/usr/bin/env python3
"""One entry point: generate, CMake, then Cargo. No recursive build managers."""
import hashlib,json,pathlib,subprocess,os,shutil
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
run(cmake,'-S','cpp','-B','build/native','-DCMAKE_BUILD_TYPE=Debug',f'-DCMAKE_CXX_COMPILER={compiler}','-DCMAKE_CXX_FLAGS=','-DCMAKE_CXX_FLAGS_DEBUG=-g','-DCMAKE_EXE_LINKER_FLAGS=','-DCMAKE_STATIC_LINKER_FLAGS=','-DCOSMOLOGY_TEST_CFITSIO=OFF')
paths=sorted([*root.glob('src/**/*.rs'),*root.glob('tests/**/*.rs'),*root.glob('cpp/**/*.cpp'),*root.glob('cpp/**/*.h'),*root.glob('cpp/**/*.hpp'),*root.glob('cpp/**/*.inc'),*root.glob('schema/*.json'),*root.glob('tools/*.py')]+[root/'Cargo.toml',root/'Cargo.lock',root/'build.rs',root/'cpp/CMakeLists.txt'])
manifest={'schema_version':1,'sources':{str(p.relative_to(root)):digest(p) for p in paths},'tools':{t:subprocess.check_output([t,'--version'],text=True).splitlines()[0] for t in ['rustc','cargo',compiler,cmake]},'compiler_executable_digest':digest(pathlib.Path(compiler)),'flags':['CMAKE_BUILD_TYPE=Debug','-g','-std=c++20','-Wall','-Wextra','-Wpedantic','-fno-fast-math','-ffp-contract=off'],'native_compile_commands':json.loads((root/'build/native/compile_commands.json').read_text()),'rust_flags':['dev','opt-level=0','debuginfo=2','panic=unwind'],'target':subprocess.check_output(['rustc','-vV'],text=True),'linker':subprocess.check_output(['ld','--version'],text=True).splitlines()[0],'standard_library':subprocess.check_output([compiler,'-print-file-name=libstdc++.so'],text=True).strip(),'standard_library_digest':digest(pathlib.Path(subprocess.check_output([compiler,'-print-file-name=libstdc++.so'],text=True).strip()).resolve()),'tool_executable_digests':{t:digest(pathlib.Path(shutil.which(t)).resolve()) for t in ['rustc','cargo']},'jobs':4,'backend':'portable_cpu','panic':'unwind','git_head':subprocess.run(['git','rev-parse','HEAD'],capture_output=True,text=True).stdout.strip() or None}
revision=manifest.pop('git_head')
manifest['build_id']=hashlib.sha256(json.dumps(manifest,sort_keys=True,separators=(',',':')).encode()).hexdigest()
manifest['git_head']=revision
manifest['git_status']=subprocess.check_output(['git','status','--porcelain'],text=True)
(root/'build').mkdir(exist_ok=True); (root/'build/build-manifest.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
run(cmake,'--build','build/native','--parallel','4')
run('cargo','build','--locked','--offline','-j','4')
