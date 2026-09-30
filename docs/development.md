# Development workflow

Start with [AGENTS.md](../AGENTS.md), [architecture](architecture.md) and [scientific contracts](scientific-contracts.md). Use [getting started](getting-started.md) to prepare the toolchain.

## A coherent slice

1. Identify the consumer and qualified prerequisites. State the equation, semantics, domain and error budget.
2. Design an independent comparison and adversarial cases. Check reference ancestry and terms.
3. Implement the smallest coherent change, with explicit model/operation identity where appropriate.
4. Add native scientific tests and boundary/CLI tests. Run the affected tests and investigate discrepancies.
5. Update the public documentation and discovery metadata. Report remaining limits separately from implemented behavior.

`schema/abi.json` owns shared ABI definitions; regenerate through the build driver. Keep changes to generated bindings alongside their schema changes. Use coarse calls and explicit buffer ownership. C++ exceptions and Rust panics must not unwind across the ABI.

## Build discipline

`python3 tools/build.py` is the integrated entry point. It configures CMake, generates bindings/manifests and builds Rust against the native library. The supported development profile currently uses Debug, conservative floating-point flags and four jobs. Changing flags, compiler, architecture or backend requires recording identity and revalidating the affected contract.

Do not mutate source while another build is reading it. Coordinate shared CPU/memory budgets. The exclusive build-identity test belongs in an isolated checkout and must run directly as described in [testing](testing.md).

The C++ namespaces and some internal identifiers retain `cosmology`. The public name is Irreducible and the binary is `irred`; do not churn stable scientific identities to match branding.

## Publishing

Use small commits with explicit paths. Review `git diff --cached`; generated data and local plans/receipts must not enter the index. Keep `Cargo.lock`, source, tests, schemas and documentation in Git. See [repository maintenance](maintenance.md).

The integration owner coordinates shared workers and publishes validated milestones under the repository's standing authorization. Other contributors submit PRs. Do not merge old pre-cleanup history or force-push without explicit authorization.
