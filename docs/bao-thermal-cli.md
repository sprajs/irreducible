# Thermal BAO CLI consumer

The operation-local source for `bao.thermal_density` is being integrated. Its
shared ABI declarations, CLI registration and allocation/runtime qualification
are pending; this guide does not yet advertise an executable command.

The consumer targets the standalone [thermal BAO](bao-thermal.md) owner. One
strict inline synthetic observation supplies ordered queries, ratios and a full
SPD covariance. One native preparation retains that observation and Gaussian
factor; a coarse batch evaluates separately fixed thermal models. The models
share an observation law, rather than forming independent observations or a
joint parameter posterior. Every model declares physical baryon/CDM/other
massless densities, Kelvin temperatures, species state weights and supplied
drag provenance. Fixed physical densities do not imply H0 cancellation. Read
[thermal observables](thermal-observables.md) for the physical identity.

For ordered residuals r=y-mu, the normalized density uses

    log p = -1/2 [r^T C^-1 r + log det C + n log(2 pi)].

Prediction, residual and density output groups retain separate native states.
A density projection refusal preserves successful requested predictions and
residuals. Unavailable groups expose no numeric payload; retained capacity still
counts. Covariance is observation noise, not a numerical-error allocation.
No jitter, covariance inflation, dropped rows or inferred independence is used.

The original three-row fixture retains C=16 times the matrix with rows
(4,1,0), (1,5,1), (0,1,6), observed ratios (24,17,22), and H0=70/72
models with the same physical inputs. It is a synthetic control. Native requested
caps remain 16 MiB total, 8 MiB per source/evaluation phase and 200 million
combined callbacks; ratio, distance/ruler and density allocations remain those
of the native guide. Species weight g is finite positive; Kelvin conversion has
its separate representability gate.

The proposed coarse C boundary owns its result until destruction, including
redshift/model-ID pools and ordered provenance. Its allocation-free minimum
owner uses disengaged scientific optionals. Proposed transport quota status7
with a null result applies only when the configured GLOBAL native cap cannot
hold that minimum owner. When the global minimum fits, phase/count/string/array
quota failures use an owned bounded numerical refusal. Allocation failure has
its separate transport status; source preparation attempted, factor completed
and thermal batch attempted remain distinct.

Rust operation-derived allocation has a separate conditional 64 MiB envelope,
including its output/specification and common-recorder copies. The preexisting
16 MiB outer CLI request/context is separate. The decoder and private-library,
serde, growth, clone and short-refusal profiles require source and measured
allocation/lifetime evidence. An unsupported profile refuses before native
execution; public type sizes alone do not qualify private container allocation.
Every completed new-operation numerical/resource refusal owns error ID
`NUMERICAL_QUALIFICATION_REQUIRED`, so the current CLI exit mapping is6.
Structural/ordinary transport failures remain exit2. No completed numerical
result, released BAO compression, physical drag, fit, posterior or observational
qualification is claimed by this integration checkpoint.
