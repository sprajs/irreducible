# Private elementary postchecks

The source contains one private C++ owner for enclosures of `sqrt`, `log` and
`exp` at exact stored `long double` arguments. It is a numerical prerequisite,
with a small synthetic consumer and permanent controls. It has no installed
header, CLI operation, C ABI or existing physical caller. This source checkpoint
has not been compiled or executed; arithmetic, counters, layouts, frames and
reference controls remain pending validation.

The first profile selects strict GNU C++20 on Linux/x86_64: little-endian IEEE
binary80, 64 significant bits and 16-byte storage. Declared-type basic operations
must round correctly, `nextafterl` must return the adjacent value, and the compiler
must preserve strict arithmetic without contraction or excess precision. Nontrapping
arithmetic, no heap allocation and the bounded opaque primitive frame are explicit
profile assumptions. Format/environment guards and finite controls support these
assumptions; they do not prove a universal compiler or library theorem. Candidate
`sqrtl`, `logl` and `expl` accuracy is not assumed.

The environment must use nearest rounding, x87 extended precision with masked
exceptions, and MXCSR nearest rounding with masked exceptions and no FTZ/DAZ.
Unsupported profiles, nonfinite values, unresolved normal margins and exhausted
budgets produce causal refusals. No jitter, fallback or tolerance relaxation is
available. A valid C++ caller owns nonoverlapping live objects; unsupported
lifetimes are not converted into numerical guarantees.

## Equations and domains

For positive-normal `x` and candidate `y`, square `y` with directed bounds and
bound the residual by `B >= |y*y-x|`. First earn `e0 >= B/y` and require `e0<y`.
Then `E >= B/(2*y-e0)` bounds `|sqrt(x)-y|`. The emitted image encloses
`[y-E,y+E]`; its endpoints must remain positive normal. The exact `x=y=1`
branch has image `[1,1]` and radius zero. A rounded square matching `x` is not
an exact-zero witness.

For `log(x)`, require `2^-128 <= x <= 2^128`. Exact power-of-two scaling gives
`x=2^e*m`, `1<=m<2`. With `z=(m-1)/(m+1)`, use
`log(m)=2*sum(z^(2j+1)/(2j+1),j=0..47)+tail`. The true `z` is at most `1/3`;
the positive omitted tail is below the exact `2^-144` upper addition. The same
graph prepares a retained `log(2)` enclosure once. Directed arithmetic owns
every series/reduction error. `log(1)` has exact image `[0,0]`; the actual
candidate discrepancy is still bounded independently.

For `exp(N)`, require a positive-normal candidate `a_hat` in the log domain.
Enclose `log(a_hat)` and earn `D >= |log(a_hat)-N|`, `0<=D<1`. Then
`E >= a_hat*D/(1-D)` bounds `|exp(N)-a_hat|`. Positive-normal directed image
endpoints are required. Exact `N=0,a_hat=1` has radius zero. Negative or exact-zero
`N` and log candidates are admitted by their stated conventions. Valid functions
may conservatively refuse when margins or positive images cannot be earned.

These are errors at **exact stored arguments**. Input/source/cast uncertainty,
physical constants, retained expansion normalization, derivatives, initial states,
matrices, history propagation and event/root errors require separate ownership.
This helper does not complete any thermal, acoustic, collisionless, predictive
controller or full-reference certificate.

## Ownership, work and wire records

Context, scratch and result are noncopyable/nonmovable and supplied in place.
The caller invokes each original candidate once after admission, retains it and
passes a const reference. Bound fields are directly emplaced after charged copies;
there are no owning floating returns or NRVO dependencies. Unavailable bounds and
uncalled candidates remain distinct from exact zero.

The nine work fields are `basic`, `next`, `scale`, `term`, `lib`, `write`, `copy`,
`guard`, `attempt`. Scalar work is `basic+next+lib`. Scaling is a basic subcount;
copies are a write subcount. Every explicit source-owned floating destination
assignment counts, including imports, constants, neighbors and final fields.
Library internals and compiler register moves are not claimed as instrumented
destinations. Their opaque primitive/profile costs remain explicit. A different
parent's destination-based scalar law needs an explicit mapping without overlap.

Each admitted action commits its fixed full increment before executing. An
unserved action preserves its exact requested increment and the served prefix.
Executed failed predicates count. Integer admission and receipt work use a finite
nonrecursive graph: per-call `128*T+128*G+32*I+512` is a conservative source-control
upper, not a hardware-instruction or elapsed-time claim. The owning fixture also
retains imports, candidate calls and caller output-gate control costs.

Records encode the ten significant IEEE80 bytes as 20 hex digits, excluding
padding. Signed-zero encodings map to mathematical zero while preserving their
sign bits. Null is unavailable; nonfinite/noncanonical encodings never become
zero. A local helper record precedes a separate charged caller gate; full consumer
availability requires that gate and successful complete-process termination.
Short/failed writes cause protocol failure. Receipt integer/byte work is separately
bounded, including refusal receipts, without floating formatting or allocation.

## Fixed consumer and permanent controls

The owning consumer uses the original synthetic dimensionless clock endpoints
`a=1e-10,.01` in that order, imported exactly from binary64. It evaluates log,
exp at the retained log candidate, and sqrt at independent exact arguments `4,16`.
It makes six candidate calls and seven helper attempts with one log context.
This pure-math fixture retains the source clock lineage; it runs no physical
background or acoustic model. Its fixed limits are 8192 scalar actions, 16384
writes, 16384 copies, 65536 semantic guards and 4096 native-owned bytes.

The original adversarial campaign retains all 63 distinct requests: the original
fixture, 38 named scalar cases, one cap-calibration preparation, 12 cap boundaries,
two overflow seeds, one inconsistent cache, one missing context, one raw alias and six environment
controls. A successor adds 17 requests: nine malformed-ledger caller entries,
three invalid imports, three exact valid imports, one unsupported import profile
and one repeated preparation. Its complete inventory is 80 requests. Failed
requests are never reset. Exact points, dyadics, near-one
cancellation, a rounded-square trap, wrong candidates, invalid domains/normality,
nonfinite inputs, exact unserved increments and no-elide ownership are included.

Every begun validation request retains its original caps, seeded ledger and available
input/candidate, status, stage, absence/refusal and served/unserved prefix.
Unassigned alias arguments remain unavailable. Checked campaign totals subtract
original seed counters before summing served work. Environment controls retain
requested and observed rounding/x87/MXCSR values. Re-preparation ends old endpoints
before refusing; no failed preparation exposes a retained interval.

The independent test uses Python `Fraction`: alternating log series after a
different range split, positive exp Taylor bounds followed by seven exact
squarings/reciprocal, and exact squared sqrt inequalities. Each bounded reference
process admits at most 64 native records/65536 bytes, with one aggregate
2-million exact-operation allowance across all of that process's objects.
Bit and exact-operation checks can refuse.
Each original native record and output digest is preserved before checking it,
including all 80 requests. Original bounded byte prefixes, output digests,
arguments, stages and actual termination statuses are retained on protocol
refusals. Every bounded process has a common final traced-memory gate, including
the fixture and report encoding. Fixture and scalar/prefix/environment/caller groups use
42 separate sequential reference processes; the largest caller group has exactly 64 records
including its source header and campaign summary. No record allowance is raised
for the full campaign. The outer controller performs no rational calculation and
retains each original bounded-process receipt before checking its completion.
The complete reference provider's internal operation/32-live-value resource gate
is explicitly withheld. Exact arithmetic controls do not silently import a full
reference pipeline or parent error/resource allocation.

The small owning fixture has its own executable, separate from the larger
validation harness. Actual GCC stack receipts feed a conservative whole active
helper/caller/encoder graph, with normal and no-elide outputs compared. The
256-byte active frame and 512-byte opaque primitive reservations, complete layout,
noheap and 4096-byte envelope must be earned without raising a failed cap.
Default source/flag mutation checks require an isolated checkout and direct
compiled binaries, with no competing build; ordinary control execution does not
authorize those mutations.

The serializer uses compiled literal text and fixed block graphs, with separate
1024-byte-write, 256-byte-read, 4096-integer-operation and 512-source-branch limits
per completed record. These count source byte assignments and conservative
integer/control graphs; they are not compiler instruction counts. Native receipts
emit the largest observed charged blocks. The final summary must itself be
dominated by those emitted extrema, or the whole protocol refuses. Ledger/cache
invariants are checked before every admission commit. Normal-or-zero binary64
imports require the selected profile before conversion. Actual runtime support
for these source bounds remains pending the five compiled contracts. Independent
source review found omitted reservation bookkeeping in the current small encoder
blocks. Wire accounting is withheld until its complete record graph is repaired;
charged extrema alone cannot establish the claimed bound.
