# Documentation

Irreducible provides physics and numerical calculations for agents working on cosmology. Its interface is non-interactive: discover a capability, supply a complete request, run it and inspect the result. Agents can extend the source when a scientific question needs a new model or reader, and use external tools for fitting or orchestration.

Read [AGENTS.md](../AGENTS.md) first, then choose a guide:

| Task | Guide |
| --- | --- |
| Development direction and missing prerequisites | [Roadmap](roadmap.md) and [gaps](gaps.md) |
| Build and run | [Getting started](getting-started.md) |
| See what is implemented | [Capabilities](capabilities.md), then `irred describe --json` |
| Make a request | [CLI contract](cli.md) |
| Interpret an output | [Run records](run-records.md) |
| Add a model, reader or calculation | [Development](development.md) and [architecture](architecture.md) |
| Check scientific assumptions | [Scientific contracts](scientific-contracts.md) |
| Run tests or add a fixture | [Testing](testing.md) and [fixture provenance](../cpp/tests/fixtures/README.md) |
| Maintain these docs | [Maintenance](maintenance.md) |

These docs travel with the source. Read the version for your commit and check the actual build's capabilities and qualification state before using a result.
