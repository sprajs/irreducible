# Shared datasets and evidence

All three repositories use [storage-layout.json](../storage-layout.json) and
[research_storage.py](../scripts/research_storage.py). Reproducible owns the
canonical transport implementation; update its identical copies, configuration
and tests in Prospector and Irreducible together. This transport performs no
scientific calculation and does not change a historical source or SDK pin.

## Where material belongs

| Material | Persistent collection |
|---|---|
| Shared input declarations and eligible unique inputs | `shared/datasets/<dataset>/<version>/` |
| Complete successful, failed or refused experiment attempts | `reproducible/experiments/<experiment>/attempts/<original-attempt>/` |
| Versioned paper acquisition declarations | `prospector/papers/<versioned-id>/` |
| Original review and acquisition evidence | `prospector/reviews/<batch>/` |
| Engine comparison/build evidence | `irreducible/evidence/<component>/<attempt>/` |
| Selected historical or unfinished material | `shared/archives/<descriptive-name>/` |

Choose a meaningful dataset/release/experiment name. Keep the original attempt
identity in provenance; a readable name cannot relabel a failed attempt. Immutable
releases live below these collections in `versions/<content-hash>/`. The hash is
an internal version identity; use collection names to navigate. The shared
catalog maps readable names to exact manifests rather than making extra copies.

`shared/catalogs/research-data.json` is a small discovery pointer. `catalog`
downloads its immutable release and verifies the release's SHA256 and VersionId.
Copy the selected entry's manifest URI, SHA256, VersionId, format and restore route
into a run's provenance. Never pin a mutable `latest`/discovery pointer. Catalog
entries can reference Git, acquisition recipes, old snapshots, exact-version
approved evidence, or private preservation bundles; use their declared route.
An old UUID snapshot is a historical directory copy. Its receipt and bytes remain
intact while readable catalog entries provide navigation. Local eviction requires exact-version readback and an unchanged-file inventory; remote history is never deleted. Completed objects and old versions have no automatic expiry.

## Keep source and working copies small

Git holds source, candidate definitions, original paper reviews, configs and
concise findings. Reference exact commit/path, including changed source identities.
Keep pinned URLs/release versions, byte lengths, hashes and acquisition instructions
for reliably downloadable inputs; exclude cached copies from preservation bundles.
Keep reconstruction instructions for cheap products. Preserve expensive results,
failed/null evidence, unavailable originals and unique unfinished work. Missing
provenance needs inspection, not a guessed role or regeneration promise.

Local `data/`, `downloads/`, `results/`, `runs/`, `papers/` and `.work/` remain
ignored working storage according to each repo's ignore rules. Restore into a new
directory, verify bytes, then pass explicit local paths to the current controller
or reader. No controller downloads from a mutable catalog or silently changes its
request/model/build identity. Exact original inputs and historical receipts stay intact in Git or verified S3 custody; local copies are disposable only after the matching recovery route passes byte checks.

## S3 is the persistent data store

Git is authoritative for source, reviews, candidates, configurations and concise
findings. The configured S3 bucket is authoritative for all useful bulk data,
full evidence, results and unfinished research state needed by a later agent.
Ignored directories are temporary scratch/cache, never a cross-run handoff.
Do not rely on another chat's local path. Read the verified catalog first; restore
only needed exact versions into fresh scratch and pass explicit verified paths
to controllers. This policy does not change historical model/request/SDK pins.

Before ending or handing off meaningful work, publish the selected completed,
failed or interrupted attempt and record its manifest URI/SHA256/VersionId/format
in the concise Git account or an S3 handoff. Unique partial work belongs in a
named preservation collection. A blocked upload must be reported as unpreserved;
retain its local bytes until resolved. Never delete the only copy. New acquisition
may use the pinned upstream source, but any useful input retained between runs
must have eligible verified S3 custody. If copying restrictions prevent custody,
record that blocker and the acquisition recipe rather than claiming it persisted.
Tools, environments and cheap products can be regenerated in scratch.

`push --evict-local` publishes first, then freshly verifies the exact completion
manifest and every remote file before removing the unchanged local regular files.
`evict` does the same for an existing pin; it defaults to verification only and
needs `--apply` to remove files. Changed/extra files, symlinks, credentials and
Git-tracked material block eviction. It does not recursively delete directories.
Keep exact restore pins before eviction. Interruption journals are local scratch;
the remote completion manifest is the durable recovery authority.

A newly published attempt need not wait for scientific catalog promotion to be
recoverable: `list --collection COLLECTION` discovers its completed manifests.
This bounded listing includes only completion markers, reports pagination and
returns exact pins. Follow `next_cursor` with `--cursor` until `complete_listing`
is true. It never treats partial uploads or the mutable catalog pointer as a
scientific pin, and it does not qualify the result.

```sh
python scripts/research_storage.py list --collection shared/archives/my-handoff
python scripts/research_storage.py push .work/my-handoff \
  --collection shared/archives/my-handoff --role mixed-evidence --rights owned \
  --provenance .work/handoff-provenance.json --evict-local
python scripts/research_storage.py evict data/restored-attempt \
  --manifest-uri s3://BUCKET/COLLECTION/versions/RELEASE/manifest.json \
  --manifest-sha256 PINNED_SHA256 --manifest-version-id PINNED_VERSION_ID --apply
```

## Commands

Install AWS CLI v2 and Python 3.11+. Credentials come from the restricted local
`research` profile or the runtime's standard AWS environment/role chain. A clone
contains public settings only. Never use root or copy local credentials to a
cloud runtime. The cloud identity reads `shared/`, `prospector/`, `irreducible/` and
`reproducible/`, and writes only `reproducible/`; it cannot publish engine evidence
or the shared catalog. Use its own restricted identity, not the local key.

```sh
python scripts/research_storage.py status
python scripts/research_storage.py catalog
python scripts/research_storage.py resolve desi-bao-dr2
python scripts/research_storage.py --ambient-credentials catalog
python scripts/research_storage.py push results/my-experiment/original-attempt \
  --collection reproducible/experiments/my-experiment/attempts/original-attempt \
  --role generated-results --rights owned --provenance .work/attempt-provenance.json --dry-run
python scripts/research_storage.py pull data/restored-attempt \
  --manifest-uri s3://BUCKET/COLLECTION/versions/RELEASE/manifest.json \
  --manifest-sha256 PINNED_SHA256 --manifest-version-id PINNED_VERSION_ID
```

Remove `--dry-run` to publish the reviewed selection. The dry run makes no network
requests. Provenance is an explicit JSON object: repository/commit and dirty state,
candidate/paper/request/input identities, actual executable/build identity, seeds,
policies, dependency/axis/calibration/selection ancestry, reconstruction route and
checks/limits where relevant. Do not supply credentials. Separate observations,
calibrated reductions, fitted summaries, synthetic controls and generated results.
`mixed-evidence` requires per-file roles in provenance. Shared ancestry and unknown
cross-covariance do not establish independence. An upload is not qualification.

The uploader screens regular bytes for common credential markers, refuses
credential/environment directories, unsafe paths, links and special files, and
uses create-only writes with SHA256 admission. **Screen compressed containers
separately before selecting them**; the transport's regular-file scan is not an
archive inspector. Never select an entire mixed `.work` directory. Files must fit
the 5 GiB single PUT limit. Both upload and restore download each exact object
version and verify all byte lengths and SHA256 hashes. The completion manifest is
published last. A collision cannot overwrite an older attempt; identical partial
uploads can be resumed. Failures retain partial objects/destinations for inspection.

`--rights private-preservation-unreviewed` records owner-authorized private custody,
not a publication license. Check explicit restrictions even for private cloud
copies and retain source notices. Public URLs do not imply permission to
redistribute. For scientific evidence subject to explicit approved selection use
Reproducible's [exact-version archive guide](https://github.com/sprajs/reproducible/blob/main/docs/archive.md)
and `archive_experiment.py`; its evidence/redistribution contract remains separate.
Historical `cloud_storage.py` snapshots remain restorable by their original route.

## Catalog stewardship

An integration owner publishes immutable catalog releases under
`shared/catalogs/releases/<sha256>.json`, checks every referenced manifest/version,
then conditionally advances the discovery pointer against its previous ETag (or
create-only for the first pointer). The catalog records `manifest_format`, exact
manifest pin, rights, scientific role, dependencies and authoritative restore
route. Do not promote an unfinished upload, substitute a scientific identity, or
drop a failed result while refreshing navigation. Catalog curation is deliberate;
a new directory upload does not automatically become an accepted dataset.

The 2026-10-08 preservation audit catalogs ignored data/results, custody archives
and unfinished worktrees across all three repositories. Its detailed manifest and
restrictions live in the private shared archive; source code and scientific pins
remain in Git. The catalog exposes concise names and exact restore routes, with
unresolved source/covariance/rights items retained. Cloud executor availability
must be tested in a real runtime before claiming cloud execution is connected.
