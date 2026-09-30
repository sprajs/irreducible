# Independent native verification

Durable verification lives in C++ tests and Rust ABI, wrapper, CLI and build-identity tests. See [testing](../../docs/testing.md) for actual commands, including the exclusive source/flag mutation test.

The exact rational Gaussian cases and handcrafted FITS fixtures here are intentional native test inputs, independent of production implementations. Their derivations and conventions belong with the headers and in [fixture provenance](../../cpp/tests/fixtures/README.md). A prepared fixture is not evidence that a future scientific consumer has been accepted.

Temporary oracle generators and bulk comparison outputs stay local. Ordinary public tests must not depend on private planning/evidence directories or an external comparison environment. Historical research assets remain local originals.
