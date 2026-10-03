#!/usr/bin/env python3
"""Independent exact-rational controls; no native libm reference calculation.

Every bounded reference invocation admits at most64 native records/65536 bytes.
The campaign controller starts separate sequential reference processes and does
no rational arithmetic. Arithmetic uses
exact Python integers/Fraction, with explicit source-operation/bit gates and a
whole-reference32MiB traced allocation gate. The provider's internal integer
temporaries are not inferred from the number of Fraction operations; complete
32-live-value/provider resource qualification remains separately withheld.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tracemalloc

MAX_BITS=1 << 20
MAX_OPERATIONS=2_000_000
reference_operations=0


class Refusal(RuntimeError):
    pass


class Exact:
    def __init__(self):
        self.operations=0
        self.maximum_bits=0

    def tick(self, count=1):
        global reference_operations
        if count>MAX_OPERATIONS-reference_operations:
            raise Refusal("reference exact-operation limit")
        reference_operations+=count
        self.operations+=count

    def bound(self, bits):
        if bits>MAX_BITS:
            raise Refusal("reference intermediate bit limit before operation")
        self.maximum_bits=max(self.maximum_bits,bits)

    def value(self, n, d=1):
        self.tick()
        self.bound(max(n.bit_length(),d.bit_length()))
        return Fraction(n,d)

    def add(self, a, b):
        self.tick()
        self.bound(max(a.numerator.bit_length()+b.denominator.bit_length(),
                       b.numerator.bit_length()+a.denominator.bit_length())+1)
        self.bound(a.denominator.bit_length()+b.denominator.bit_length())
        return a+b

    def mul(self, a, b):
        self.tick()
        self.bound(a.numerator.bit_length()+b.numerator.bit_length())
        self.bound(a.denominator.bit_length()+b.denominator.bit_length())
        return a*b

    def div(self, a, b):
        if not b:
            raise Refusal("reference division by zero")
        self.tick()
        self.bound(a.numerator.bit_length()+b.denominator.bit_length())
        self.bound(a.denominator.bit_length()+b.numerator.bit_length())
        return a/b

    def neg(self, a):
        self.tick()
        return -a

    def decode(self, text):
        if text is None:
            raise Refusal("unavailable scalar")
        if not isinstance(text,str) or len(text)!=20:
            raise Refusal("not a significant binary80 encoding")
        raw=bytes.fromhex(text)
        significand=int.from_bytes(raw[:8],"little")
        sign_exponent=int.from_bytes(raw[8:],"little")
        exponent=sign_exponent & 0x7fff
        if exponent==0 and significand==0:
            return self.value(0)  # Wire sign is preserved in the original record.
        if exponent==0 or exponent==0x7fff or not significand >> 63:
            raise Refusal("nonfinite/noncanonical/subnormal scalar")
        shift=exponent-16383-63
        self.bound(max(64+max(shift,0),1+max(-shift,0)))
        if max(64+max(shift,0),1+max(-shift,0))>16600:
            raise Refusal("decoded input bit limit")
        n=significand if not sign_exponent >> 15 else -significand
        return self.value(n << shift) if shift>=0 else self.value(n,1 << -shift)

    def decode64(self, text):
        if not isinstance(text,str) or len(text)!=16:
            raise Refusal("not a binary64 original coordinate")
        bits=int.from_bytes(bytes.fromhex(text),"little")
        exponent=(bits >> 52) & 2047
        significand=bits & ((1 << 52)-1)
        if exponent==0 and significand==0:
            return self.value(0)
        if exponent==0 or exponent==2047:
            raise Refusal("nonfinite/subnormal binary64 coordinate")
        significand |= 1 << 52
        shift=exponent-1023-52
        n=significand if not bits >> 63 else -significand
        return self.value(n << shift) if shift>=0 else self.value(n,1 << -shift)

    def log_piece(self, t, terms=256):
        # Alternating partial sums, independent of the native atanh graph.
        if not 0<=t<=self.value(1,2):
            raise Refusal("alternating reference domain")
        power=t
        total=self.value(0)
        for j in range(1,terms+1):
            term=self.div(power,self.value(j))
            total=self.add(total,term if j & 1 else self.neg(term))
            power=self.mul(power,t)
        next_term=self.div(power,self.value(terms+1))
        adjacent=self.add(total,next_term if (terms+1) & 1 else self.neg(next_term))
        return min(total,adjacent),max(total,adjacent)

    def log(self, x):
        if not self.value(1,1 << 128)<=x<=self.value(1 << 128):
            raise Refusal("reference log domain")
        m=x
        exponent=0
        two=self.value(2)
        while m<1:
            m=self.mul(m,two)
            exponent-=1
        while m>=2:
            m=self.div(m,two)
            exponent+=1
        if m==1:
            lower=upper=self.value(0)
        elif m<=self.value(3,2):
            lower,upper=self.log_piece(self.add(m,self.value(-1)))
        else:
            first=self.log_piece(self.value(1,2))
            second=self.log_piece(self.add(self.mul(m,self.value(2,3)),self.value(-1)))
            lower=self.add(first[0],second[0])
            upper=self.add(first[1],second[1])
        if exponent:
            first=self.log_piece(self.value(1,2))
            second=self.log_piece(self.value(1,3))
            l2=self.add(first[0],second[0])
            u2=self.add(first[1],second[1])
            e=self.value(exponent)
            lower=self.add(lower,self.mul(e,l2 if exponent>0 else u2))
            upper=self.add(upper,self.mul(e,u2 if exponent>0 else l2))
        return lower,upper

    def exp(self, n):
        if abs(n)>128:
            raise Refusal("reference exp domain")
        if n==0:
            return self.value(1),self.value(1)
        t=self.div(n if n>=0 else self.neg(n),self.value(128))
        term=self.value(1)
        total=self.value(1)
        for j in range(1,97):
            term=self.div(self.mul(term,t),self.value(j))
            total=self.add(total,term)
        omitted=self.div(self.mul(term,t),self.value(97))
        denominator=self.add(self.value(1),self.neg(self.div(t,self.value(98))))
        upper=self.add(total,self.div(omitted,denominator))
        lower=total
        for _ in range(7):
            lower=self.mul(lower,lower)
            upper=self.mul(upper,upper)
        if n<0:
            return self.div(self.value(1),upper),self.div(self.value(1),lower)
        return lower,upper


def native(binary, arguments):
    # Read admission precedes allocation beyond the fixed byte allowance.
    process=subprocess.Popen([str(binary),*arguments],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    output=process.stdout.read(65537)
    if len(output)>65536:
        process.kill()
        status=process.wait()
        print(json.dumps({"native_failure":{"arguments":arguments,"stage":"byte_limit","status":status,
                         "prefix_only":True,"original_output_hex":output.hex(),
                         "output_sha256":hashlib.sha256(output).hexdigest()}}),flush=True)
        raise Refusal("native wire byte limit")
    status=process.wait()
    print(json.dumps({"native_output":{"arguments":arguments,"stage":"collected","status":status,
                     "prefix_only":False,"original_output_hex":output.hex(),
                     "output_sha256":hashlib.sha256(output).hexdigest()}}),flush=True)
    if status:
        raise Refusal(f"native refused invocation status={status}; original receipt preserved")
    lines=output.splitlines()
    if len(lines)>64:
        print(json.dumps({"native_failure":{"arguments":arguments,"stage":"record_limit","status":status,
                         "records":len(lines),"output_sha256":hashlib.sha256(output).hexdigest()}}),flush=True)
        raise Refusal("native wire record limit")
    return [json.loads(line) for line in lines],hashlib.sha256(output).hexdigest()


def validate(record):
    exact=Exact()
    _,operation,status,_,xtext,ytext,bound,work,served,refusal,control=record
    if operation not in (0,1,2):
        raise Refusal("literal elementary operation identity")
    if control[3]!=128*control[0]+128*control[1]+32*control[2]+512:
        raise Refusal("integer control upper identity")
    if served[6]>served[5] or served[2]>served[0]:
        raise Refusal("overlapping subcount identity")
    if status!=1:
        if bound is not None or not refusal[0]:
            raise Refusal("native refusal without absence/causal receipt")
        return {"scope":"causal native refusal","exact_operations":0}
    if bound is None or ytext is None:
        raise Refusal("missing accepted scalar")
    x=exact.decode(xtext)
    y=exact.decode(ytext)
    lower,upper,error=(exact.decode(t) for t in bound)
    if lower>upper or error<0:
        raise Refusal("bound conventions")
    if operation==2:
        if lower<=0 or exact.mul(lower,lower)>x or exact.mul(upper,upper)<x:
            raise Refusal("independent squared image")
        ball_lower=exact.add(y,exact.neg(error))
        ball_upper=exact.add(y,error)
        if ball_lower<=0 or exact.mul(ball_lower,ball_lower)>x or exact.mul(ball_upper,ball_upper)<x:
            raise Refusal("independent squared radius")
    else:
        reference=exact.log(x) if operation==0 else exact.exp(x)
        if not lower<=reference[0]<=reference[1]<=upper:
            raise Refusal("reference interval not contained in image")
        if max(abs(exact.add(y,exact.neg(reference[0]))),
               abs(exact.add(reference[1],exact.neg(y))))>error:
            raise Refusal("reference discrepancy exceeds radius")
    return {"scope":"exact stored-argument image/radius","exact_operations":exact.operations,
            "maximum_intermediate_bits":exact.maximum_bits,
            "complete_provider_work_live32_qualified":False}


def controls(records,expected_ids):
    starts={r[1]:r for r in records if r[0]=="request"}
    finishes={r[1]:r for r in records if r[0]=="served"}
    if set(starts)!=set(expected_ids) or set(finishes)!=set(expected_ids):
        raise Refusal("complete per-request start/served lineage")
    if sum(r[0]=="request" for r in records)!=len(starts) or sum(r[0]=="served" for r in records)!=len(finishes):
        raise Refusal("duplicate request lineage")
    total=[0]*9
    scalar_total=0
    checked=[]
    for request_id in expected_ids:
        start=starts[request_id]
        events=[r for r in records if len(r)>1 and r[1]==request_id and r[0] in ("context","caller")]
        scalar_events=[r for r in records if isinstance(r[0],int) and r[0]==request_id]
        if not events and not scalar_events:
            raise Refusal("request final output/refusal missing")
        final=scalar_events[-1][8] if scalar_events else events[-1][7]
        delta=[b-a for a,b in zip(start[5],final)]
        served=finishes[request_id]
        if any(d<0 for d in delta) or delta!=served[2] or any(v>(1 << 64)-1 for v in final):
            raise Refusal("seed-preserving per-request delta")
        if served[3]!=sum(delta[i] for i in (0,1,4)):
            raise Refusal("per-request scalar delta")
        caps=start[4]
        # Malformed seeds are preserved inputs. Refusal cannot authorize a new
        # action; no cap/sum meaning is invented for an invalid original cache.
        if start[2]!=11:
            if sum(final[i] for i in (0,1,4))>caps[0] or final[5]>caps[1] or final[6]>caps[2] or final[7]>caps[3]:
                raise Refusal("served original cap")
        elif any(delta):
            raise Refusal("malformed caller ledger served an action")
        for i,d in enumerate(delta):
            if d>(1 << 64)-1-total[i]:
                raise Refusal("aggregate work overflow")
            total[i]+=d
        scalar_total+=served[3]
        for event in events:
            if event[0]=="context":
                _,_,phase,status,stage,image,call_work,prefix,refusal,control=event
                if status!=1:
                    if image is not None or not refusal[0]:
                        raise Refusal("context refusal exposes endpoints or misses receipt")
                    checked.append({"request":request_id,"phase":phase,"scope":"causal context refusal"})
                else:
                    exact=Exact()
                    lower,upper=(exact.decode(x) for x in image)
                    lo,hi=exact.log(exact.value(2))
                    if not lower<=lo<=hi<=upper:
                        raise Refusal("independent prepared log2 enclosure")
                    checked.append({"request":request_id,"phase":phase,"scope":"exact prepared log2 enclosure",
                                    "exact_operations":exact.operations,"maximum_intermediate_bits":exact.maximum_bits})
                if control[3]!=128*control[0]+128*control[1]+32*control[2]+512:
                    raise Refusal("context control upper identity")
            elif event[2]==1:
                originals=[r for r in records if r[0]=="coordinate" and r[1]==request_id]
                if len(originals)!=1:
                    raise Refusal("original binary64 input receipt")
                if event[3]:
                    exact=Exact()
                    if event[4] is None or exact.decode(event[4])!=exact.decode64(originals[0][2]):
                        raise Refusal("exact guarded original coordinate import")
                    original_sign=int.from_bytes(bytes.fromhex(originals[0][2]),"little") >> 63
                    emitted_sign=int.from_bytes(bytes.fromhex(event[4])[8:],"little") >> 15
                    if original_sign!=emitted_sign:
                        raise Refusal("original coordinate sign including exact signed zero")
                elif event[4] is not None:
                    raise Refusal("refused import output must remain unavailable")
        for event in scalar_events:
            checked.append({"request":request_id,**validate(event)})
    summaries=[r for r in records if r[0]=="campaign"]
    if summaries and (len(summaries)!=1 or summaries[0][1]!=len(expected_ids) or
                      summaries[0][2]!=len(expected_ids) or summaries[0][3]!=total or summaries[0][4]!=scalar_total):
        raise Refusal("whole campaign actual delta total")
    return checked


def source_identity(records,expected):
    sources=[r for r in records if r[0]=="control-source"]
    if len(sources)!=1 or sources[0][1]!=expected:
        raise Refusal("actual native control source identity")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("binary",type=Path)
    parser.add_argument("controls",type=Path)
    parser.add_argument("--scope",choices=("fixture","case","prefixes","environments","entry"))
    parser.add_argument("--case-index",type=int)
    parser.add_argument("--expected-source")
    args=parser.parse_args()
    if args.scope is None:
        # Campaign controller performs NO rational numerical work. Each spawned
        # bounded reference invocation owns its original <=64-record allowance.
        source=None;reports=[]
        plans=[("fixture",None),*(("case",i) for i in range(38)),
               ("prefixes",None),("environments",None),("entry",None)]
        for scope,index in plans:
            command=[sys.executable,str(Path(__file__).resolve()),str(args.binary),str(args.controls),"--scope",scope]
            if index is not None:command.extend(("--case-index",str(index)))
            if source is not None:command.extend(("--expected-source",source))
            process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            output=process.stdout.read(1024*1024+1)
            sys.stdout.buffer.write(output);sys.stdout.flush() # Preserve failures too.
            if len(output)>1024*1024:
                process.kill();status=process.wait()
                print(json.dumps({"bounded_reference_failure":{"scope":scope,"case_index":index,
                                 "stage":"receipt_byte_limit","status":status,"prefix_only":True,
                                 "original_output_hex":output.hex(),
                                 "output_sha256":hashlib.sha256(output).hexdigest()}}),flush=True)
                raise Refusal("bounded reference receipt byte limit")
            status=process.wait()
            print(json.dumps({"bounded_reference_output":{"scope":scope,"case_index":index,
                             "status":status,"prefix_only":False,"collected_bytes":len(output),
                             "output_sha256":hashlib.sha256(output).hexdigest()}}),flush=True)
            if status:
                raise Refusal("bounded reference invocation/protocol refused")
            report=json.loads(output.splitlines()[-1])
            if source is None:source=report["native_source_sha256"]
            if report["native_source_sha256"]!=source:
                raise Refusal("bounded reference source mismatch")
            reports.append(report)
        if sum(r["distinct_requests"] for r in reports)!=80:
            raise Refusal("complete80-request reference campaign")
        print(json.dumps({"scope":"complete bounded reference campaign","distinct_requests":80,
                          "original_distinct_requests_preserved":63,"reference_invocations":len(reports),
                          "native_source_sha256":source,"provider_internal_work_live32_qualified":False,
                          "whole_physics_or_reference_certificate":False,"reports":reports}))
        return
    tracemalloc.start()
    accepted=refused=0
    if args.scope=="fixture":
        records,digest=native(args.binary,["--facts"])
        source=records[0][3] if records and records[0][0]=="ieee80-le" else None
    else:
        records=[];digest=None;source=args.expected_source
    if not isinstance(source,str) or len(source)!=64:
        raise Refusal("actual original fixture source identity")
    receipts=[]
    invocations=[]
    if args.scope=="fixture":
        invocations.append({"arguments":["--facts"],"output_sha256":digest,"native_records":records})
        print(json.dumps(invocations[-1]),flush=True) # Preserve originals BEFORE validation.
    for record in records:
        if isinstance(record[0],int):
            receipts.append(validate(record))
            accepted+=record[2]==1
            refused+=record[2]!=1
    gates=[r for r in records if r[0]=="gate"]
    if args.scope=="fixture" and (accepted!=6 or refused or len(records)!=15 or len(gates)!=6 or not all(r[3]==1 for r in gates)):
        raise Refusal("original two-row fixture record identity")
    # Exact mandatory source corpus:38 separate native requests, never a reset.
    if args.scope=="case" and (args.case_index is None or not 0<=args.case_index<38):
        raise Refusal("literal case index")
    for case in (args.case_index,) if args.scope=="case" else ():
        records,case_digest=native(args.controls,["--case",str(case)])
        invocations.append({"arguments":["--case",str(case)],"output_sha256":case_digest,"native_records":records})
        print(json.dumps(invocations[-1]),flush=True)
        source_identity(records,source)
        scalar=[r for r in records if isinstance(r[0],int)]
        if len(scalar)!=1 or scalar[0][0]!=case+1:
            raise Refusal("one-case original identity")
        receipts.extend(controls(records,[case+1]))
        accepted+=scalar[0][2]==1
        refused+=scalar[0][2]!=1
        if tracemalloc.get_traced_memory()[1]>32*1024*1024:
            raise Refusal("whole traced reference payload limit")
    groups={"prefixes":("--prefixes",range(39,57)),"environments":("--environments",range(57,63)),
            "entry":("--entry-controls",range(63,80))}
    for argument,ids in (groups[args.scope],) if args.scope in groups else ():
        records,case_digest=native(args.controls,[argument])
        invocations.append({"arguments":[argument],"output_sha256":case_digest,"native_records":records})
        print(json.dumps(invocations[-1]),flush=True)
        source_identity(records,source)
        receipts.extend(controls(records,list(ids)))
        if tracemalloc.get_traced_memory()[1]>32*1024*1024:
            raise Refusal("whole traced reference payload limit")
    print(json.dumps({"scope":"independent exact stored-argument controls",
                      "accepted_scalar_records":accepted,"refused_scalar_records":refused,
                      "fixture_output_sha256":digest,"reference_source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      "peak_traced_bytes":tracemalloc.get_traced_memory()[1],
                      "provider_internal_work_live32_qualified":False,
                      "native_source_sha256":source,"bounded_scope":args.scope,
                      "distinct_requests":len(groups[args.scope][1]) if args.scope in groups else 1,
                      "aggregate_exact_operations":reference_operations,
                      "whole_physics_or_reference_certificate":False,"receipts":receipts,"invocations":invocations}))
    # Full bounded-process gate, including fixture and complete report encoding.
    # The original report is preserved before a failing whole-process disposition.
    final_peak=tracemalloc.get_traced_memory()[1]
    if final_peak>32*1024*1024:
        print(json.dumps({"reference_failure":{"stage":"whole_traced_payload","peak_traced_bytes":final_peak,
                         "maximum_traced_bytes":32*1024*1024,"bounded_scope":args.scope}}),flush=True)
        raise Refusal("whole bounded reference traced payload limit")


if __name__=="__main__":
    try:
        main()
    except (Refusal,ValueError,AssertionError) as error:
        print(f"elementary reference refused: {error}",file=sys.stderr)
        sys.exit(1)
