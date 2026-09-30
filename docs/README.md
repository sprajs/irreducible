# Irreducible documentation

**Agents start here.** These files are the canonical documentation shipped with the source. Read [AGENTS.md](../AGENTS.md), choose a task below, then discover the executable you will actually run. The public checkout is self-contained; private plans and development receipts are not prerequisites.

| Task | Read first | Then inspect |
| --- | --- | --- |
| Build the development executable | [Getting started](getting-started.md) | [Testing](testing.md) |
| Find out what exists | [Capabilities and limits](capabilities.md) | `target/debug/irred describe --json` |
| Run JSON requests | [CLI contract](cli.md) | [Run records](run-records.md) |
| Change an equation or model | [Scientific contracts](scientific-contracts.md) | [Development](development.md) and affected tests |
| Understand the source layout | [Architecture](architecture.md) | [ABI schema](../schema/abi.json) |
| Add or review a fixture | [Testing](testing.md) | [Fixture provenance](../cpp/tests/fixtures/README.md) |
| Submit a change | [Contributing](../CONTRIBUTING.md) | [Development](development.md) |
| Maintain public documentation | [Repository maintenance](maintenance.md) | [Wiki source](wiki/Home.md) |

## First actions

```sh
# After the dependency preparation in getting-started.md:
python3 tools/build.py
target/debug/irred describe --json
target/debug/irred run tests/fixtures/exact-add.json runs/agent-smoke
```

Inspect operation IDs, build identity, qualification state, receipt and exit code. Do not infer an operation from a roadmap or filename. If a capability is absent, report the gap or develop it with its prerequisites.

Keep four questions separate: did the run execute, is the numerical calculation qualified for this use, is the inference adequate, and what interpretation follows under the assumptions? A passing fixture answers only the question it actually tests.

This documentation travels with the code. Read the version at the commit you use and update it when a public contract changes. The [Wiki](https://github.com/sprajs/irreducible/wiki) points here so that instructions do not diverge.
