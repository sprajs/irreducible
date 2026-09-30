# Cosmology

Cosmology is a planned agent-operated scientific engine: a Rust agent/data layer and one owned C++20 scientific library for observations, shared physics, numerical calculations and inference. CUDA and distributed execution are planned where their accuracy and performance can be demonstrated.

**Current state: planning only.** No scientific executable, supported production command or qualified numerical implementation exists yet.

- [Main plan](Plan/README.md) explains the design and links the topic owners.
- [Build order](Plan/BUILD-ORDER.md) gives the sequential starting path and comparison gates.
- [Status](Plan/STATUS.md) distinguishes planning work from implementation and scientific qualification.
- [Agent instructions](AGENTS.md) cover both developing the system and conducting analyses with it.

The project begins with mathematical and statistical foundations, then reproduces selected workflows from the earlier supernova research. That research is a source of datasets, cases and known limitations; its outputs are not automatically scientific truth or independently verified results of this project.
