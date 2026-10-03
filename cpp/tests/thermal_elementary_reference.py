#!/usr/bin/env python3
"""Independent exact-rational controls; no native libm reference calculation.

Every native invocation admits at most64 records/65536 bytes. Arithmetic uses
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


class Refusal(RuntimeError):
    pass


class Exact:
    def __init__(self):
        self.operations=0
        self.maximum_bits=0

    def tick(self, count=1):
        if count>MAX_OPERATIONS-self.operations:
            raise Refusal("reference exact-operation limit")
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
        process.wait()
        raise Refusal("native wire byte limit")
    status=process.wait()
    if status:
        raise Refusal(f"native refused invocation status={status}; preserved output={output!r}")
    lines=output.splitlines()
    if len(lines)>64:
        raise Refusal("native wire record limit")
    return [json.loads(line) for line in lines],hashlib.sha256(output).hexdigest()


def validate(record):
    exact=Exact()
    _,operation,status,_,xtext,ytext,bound,work,served,refusal,control=record
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
    if control[3]!=128*control[0]+128*control[1]+32*control[2]+512:
        raise Refusal("integer control upper identity")
    if served[6]>served[5] or served[2]>served[0]:
        raise Refusal("overlapping subcount identity")
    return {"scope":"exact stored-argument image/radius","exact_operations":exact.operations,
            "maximum_intermediate_bits":exact.maximum_bits,
            "complete_provider_work_live32_qualified":False}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("binary",type=Path)
    parser.add_argument("controls",type=Path)
    args=parser.parse_args()
    tracemalloc.start()
    accepted=refused=0
    records,digest=native(args.binary,["--facts"])
    receipts=[]
    for record in records:
        if isinstance(record[0],int):
            receipts.append(validate(record))
            accepted+=record[2]==1
            refused+=record[2]!=1
    gates=[r for r in records if r[0]=="gate"]
    if accepted!=6 or refused or len(records)!=15 or len(gates)!=6 or not all(r[3]==1 for r in gates):
        raise Refusal("original two-row fixture record identity")
    # Exact mandatory source corpus:38 separate native requests, never a reset.
    for case in range(38):
        records,case_digest=native(args.controls,["--case",str(case)])
        if len(records)!=1:
            raise Refusal("one-case wire identity")
        receipt=validate(records[0])
        accepted+=records[0][2]==1
        refused+=records[0][2]!=1
        receipts.append(receipt)
        if tracemalloc.get_traced_memory()[1]>32*1024*1024:
            raise Refusal("whole traced reference payload limit")
    print(json.dumps({"scope":"independent exact stored-argument controls",
                      "accepted_records":accepted,"preserved_native_refusals":refused,
                      "fixture_output_sha256":digest,"reference_source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      "peak_traced_bytes":tracemalloc.get_traced_memory()[1],
                      "provider_internal_work_live32_qualified":False,
                      "whole_physics_or_reference_certificate":False,"receipts":receipts}))


if __name__=="__main__":
    try:
        main()
    except (Refusal,ValueError,AssertionError) as error:
        print(f"elementary reference refused: {error}",file=sys.stderr)
        sys.exit(1)
