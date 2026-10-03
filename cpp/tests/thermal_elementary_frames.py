#!/usr/bin/env python3
"""Conservative whole owning-fixture frame gate, with raw GCC .su receipts.

All nonentry helper frames are added, even when mutually exclusive. Encoder
methods use their literal source call DAG; helper and encoding phases are
sequential. Opaque primitive512 bytes remains an explicit profile
assumption, rather than a measurement of an external library's whole stack.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def frames(path):
    if not path.is_file():
        raise RuntimeError(f"missing actual compiler frame receipt: {path}")
    result=[]
    for line in path.read_text().splitlines():
        signature,amount,kind=line.split("\t")
        if kind!="static":
            raise RuntimeError(f"unbounded/dynamic frame: {line}")
        result.append((signature,int(amount)))
    return result


def run(binary):
    completed=subprocess.run([str(binary),"--facts"],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if completed.returncode or len(completed.stdout)>65536:
        raise RuntimeError(f"native lifetime/frame invocation refused: {completed.stdout!r}")
    records=[json.loads(line) for line in completed.stdout.splitlines()]
    if len(records)!=15 or records[-1][0]!="summary":
        raise RuntimeError("complete owning fixture receipt missing")
    return completed.stdout,records[-1]


def qualify(build,target,summary):
    core_target="test_thermal_elementary_no_elide" if target.endswith("no_elide") else "irred_core"
    helper_path=build/"CMakeFiles"/(core_target+".dir")/"src/thermal_elementary_postcheck.cpp.su"
    caller_path=build/"CMakeFiles"/(target+".dir")/"tests/test_thermal_elementary.cpp.su"
    helper=frames(helper_path)
    caller=frames(caller_path)
    entry_names=("postcheck_sqrt(","postcheck_log(","postcheck_exp(","prepare_elementary_log(",
                 "admit_elementary_candidate(","import_elementary_coordinate(","admit_elementary_output_guard(")
    entries=[size for name,size in helper if any(tag in name for tag in entry_names)]
    other=[size for name,size in helper if not any(tag in name for tag in entry_names)]
    if not entries:
        raise RuntimeError("compiled helper entry frames missing")
    request_bytes,wire_bytes=summary[-2:]
    main=[size for name,size in caller if "main(" in name]
    fixture=[size for name,size in caller if "fixture(" in name]
    if not main:
        raise RuntimeError("compiled owning caller frame missing")
    # Separate functions sum all active enclosing frames. If fixture was inlined,
    # the main frame owns both fixed objects. Never subtract an absent object.
    if fixture:
        caller_extra=max(0,max(main)-wire_bytes)+max(0,max(fixture)-request_bytes)
    else:
        caller_extra=max(0,max(main)-wire_bytes-request_bytes)
    dag={"reserve":(),"character":("reserve",),"text":("reserve","character"),
         "integer":("reserve","character"),"wide":("reserve","character"),
         "work":("character","integer"),"control":("character","integer"),
         "aggregate_control":("character","integer"),"flush":("character",),
         "result":("character","integer","wide","work","control","flush"),
         "refused":("text","integer","reserve","character","wide","work","control","flush"),
         "gate":("text","integer","character","work","control","flush"),"add_control":()}
    sizes={method:max([size for name,size in caller if f"Wire::{method}(" in name],default=0) for method in dag}
    def chain(method):
        return sizes[method]+max([chain(child) for child in dag[method]],default=0)
    encoder=max(chain(method) for method in dag)
    # Zero-allocation observer and failing assertion/I/O observers are separate
    # validation paths. They cannot qualify a successful fixture via this gate.
    if summary[4] or summary[5]:
        raise RuntimeError("allocation observation cannot qualify frame profile")
    unknown=[(name,size) for name,size in caller if "main(" not in name and "fixture(" not in name
             and "Wire::" not in name and "operator new" not in name and "operator delete" not in name
             and "__wrap_" not in name and "check(" not in name]
    active=max(max(entries)+sum(other),encoder)+caller_extra+sum(size for _,size in unknown)
    peak=request_bytes+wire_bytes+active+512+176+384
    receipt={"target":target,"helper_raw":helper,"caller_raw":caller,
             "active_helper_control_copy_frame_upper":active,"native_owned_peak_upper":peak,
             "encoder_literal_DAG_frame_upper":encoder,"conservative_unknown_caller_frames":unknown,
             "opaque_primitive512_profile_assumption":True,
             "compiler_receipt_sha256":[hashlib.sha256(p.read_bytes()).hexdigest() for p in (helper_path,caller_path)]}
    print(json.dumps(receipt))
    if active>256 or peak>4096:
        raise RuntimeError("whole active frame/4096 admission refused; raw failure retained above")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("build",type=Path)
    parser.add_argument("normal",type=Path)
    parser.add_argument("no_elide",type=Path)
    args=parser.parse_args()
    normal,normal_summary=run(args.normal)
    no_elide,no_elide_summary=run(args.no_elide)
    if normal!=no_elide:
        raise RuntimeError("no-elide outputs/counters/wire/layout differ")
    qualify(args.build,"test_thermal_elementary_fixture",normal_summary)
    qualify(args.build,"test_thermal_elementary_no_elide",no_elide_summary)


if __name__=="__main__":
    try:
        main()
    except (RuntimeError,ValueError) as error:
        print(f"elementary frame refusal: {error}",file=sys.stderr)
        sys.exit(1)
