# Contributing

Help turn a scientific question into a calculation we can trust. Models, readers, numerical fixes, independent tests and clearer documentation are all welcome, including work written with an agent.

Work on a branch and open a coherent pull request. A substantial change is welcome when its parts belong together and a reader can understand and test it. Use useful intermediate commits and explain dependencies; split unrelated work rather than imposing a line-count limit. Explain the problem, what changes and how you checked it. Keep code, tests and the documentation needed to use the change together. Push work branches promptly after each meaningful checkpoint commit for backup, including before a PR or handover. The owner gives agents standing authorization to create PRs and merge ready changes after a deliberate review pass and green required CI on the latest integrated candidate; direct pushes to `main` still need a specific request. CI checks PRs and merged `main`, rather than every unpublished branch push. Agents may review their own PRs in a separate deliberate pass, fix actionable findings and merge without asking again. Green checks alone do not replace review; a separate human reviewer is not required. Do not force-push without authorization.

For physics or statistics, state the assumptions, supported domain and error budget, and bring a comparison that can challenge the implementation. Preserve original inputs and failed results. Documentation-only changes need working links and truthful examples, rather than a full scientific test campaign.

[Development](docs/development.md) explains the workflow; [testing](docs/testing.md) lists the checks. Agents should also read [AGENTS.md](AGENTS.md).
