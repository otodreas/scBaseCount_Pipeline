# GCS to R2 migration

The raw scBaseCount `h5ad` files used by this project were copied from Google Cloud Storage to Cloudflare R2 in two migrations on the Lund University bioinformatics server. The first migration created the initial mirror on May 6, 2026. The second added datasets selected by the revised metadata on July 30, 2026. The R2 mirror avoids repeated GCS downloads and does not incur R2 egress charges.

The server log at `logs/migrate_gcs_to_r2.log` was ignored by Git and remains on that server. The committed CSV manifests listed below are the durable run records.

## Historical repository snapshots

The current migration runner does not reproduce both historical selections. In particular, `pipelines/migrate_gcs_to_r2.py` now subtracts a baseline CSV from a source CSV. The May runner instead processed every row in one CSV.

At the May snapshot, the runner was located at the repository root as `migrate_gcs_to_r2.py`. It accepted `--datasets` and `--dry-run`, defaulted to `output/metadata/datasets.csv`, and selected rows by loading that CSV and iterating over it directly:

```python
datasets = pd.read_csv(args.datasets)

for n, (_, row) in enumerate(datasets.iterrows(), start=1):
```

The historical selection is visible at [lines 81-88 of the May runner](https://github.com/otodreas/scBaseCount_Pipeline/blob/d368c023d64368df65c6f66cfab42f1877d6bef9/migrate_gcs_to_r2.py#L81-L88). It had no `--baseline` argument, CSV validation, or accession subtraction. Those were added in commit [`4b92630`](https://github.com/otodreas/scBaseCount_Pipeline/commit/4b92630) on July 30 before the second migration.

Use the post-run snapshot for the migration being inspected:

| Migration | Repository snapshot | Checkout |
| --- | --- | --- |
| May 6 initial mirror | [d368c02](https://github.com/otodreas/scBaseCount_Pipeline/commit/d368c023d64368df65c6f66cfab42f1877d6bef9) | `git checkout d368c023d64368df65c6f66cfab42f1877d6bef9` |
| July 30 metadata delta | [68c0425](https://github.com/otodreas/scBaseCount_Pipeline/commit/68c04253addf9fd01ca6671771a36628b50afab9) | `git checkout 68c04253addf9fd01ca6671771a36628b50afab9` |

These snapshots were committed after their transfers. They preserve the code, input CSVs, dependency lock, and run manifests together. Manifest timestamps use the server's local clock and do not include a timezone.

## GCS access at the time

On December 18, 2025, Arc moved its documented access path to `gs://arc-institute-virtual-cell-atlas` in [PR #21](https://github.com/ArcInstitute/arc-virtual-cell-atlas/pull/21). The scBaseCount documentation available during both migrations described the bucket as Requester Pays and said:

> **Note**: The new bucket is subject to Requester Pays. Users can access up to 2TB of data per month for free before fees apply.

That wording appears in the scBaseCount README snapshots current on [May 6, 2026](https://github.com/ArcInstitute/arc-virtual-cell-atlas/blob/48603659fe2ade696986605c6e22fc7ccbca4f6a/scBaseCount/README.md#L4-L9) and [July 30, 2026](https://github.com/ArcInstitute/arc-virtual-cell-atlas/blob/cfcdd6a7709ac25a1e347f3b709fea2182c6b7b7/scBaseCount/README.md#L4-L9).

Arc's documentation was internally inconsistent. Its [Python tutorial at the time](https://github.com/ArcInstitute/arc-virtual-cell-atlas/blob/48603659fe2ade696986605c6e22fc7ccbca4f6a/scBaseCount/tutorial-py.ipynb) created a GCS filesystem without a requester billing project and showed `gsutil` commands without the `-u` billing-project option. Google documents that a strictly enforced [Requester Pays](https://cloud.google.com/storage/docs/requester-pays) bucket rejects requests that supply neither a billing project nor owner-level billing permission.

Both migrations in this repository used anonymous Google credentials without a requester billing project, and their manifests show that the downloads succeeded. Anonymous public access therefore worked without requester-side GCS charges on May 6 and July 30, 2026. This records the bucket's observed behavior, but does not establish whether Arc had enabled its Requester Pays flag or what source-side costs Arc incurred. These migrations did not claim the 2 TB Marketplace allowance through a subscribed project.

On August 28, 2026, Arc clarified the access requirements in [PR #25](https://github.com/ArcInstitute/arc-virtual-cell-atlas/pull/25) after a user downloaded 1.44 TiB through an unsubscribed billing project and was charged standard GCS rates. The [current Arc instructions](https://github.com/ArcInstitute/arc-virtual-cell-atlas/blob/main/README.md#accessing-the-data) require users to subscribe a billing-enabled project through Marketplace and pass that exact project with each request. That clarification concerns requests explicitly billed to a project and does not explain why historical anonymous access succeeded.

The historical anonymous transfer worked, but `scripts/storage/gcs.py` does not implement Arc's current documented access procedure. Do not assume it will keep working if the bucket configuration changes.

## May 6 initial mirror

The May runner processed all 801 rows in its source CSV without a baseline. It uploaded 798 objects and skipped three objects that already existed in R2 with matching GCS MD5 metadata.

Tracked inputs, code, and outputs in snapshot `d368c02`:

- Source: [output/metadata/datasets.csv](https://github.com/otodreas/scBaseCount_Pipeline/blob/d368c023d64368df65c6f66cfab42f1877d6bef9/output/metadata/datasets.csv), containing 801 rows in this snapshot.
- Runner: [migrate_gcs_to_r2.py](https://github.com/otodreas/scBaseCount_Pipeline/blob/d368c023d64368df65c6f66cfab42f1877d6bef9/migrate_gcs_to_r2.py).
- Storage code: [scripts/gcs/client.py](https://github.com/otodreas/scBaseCount_Pipeline/blob/d368c023d64368df65c6f66cfab42f1877d6bef9/scripts/gcs/client.py) and [scripts/r2/client.py](https://github.com/otodreas/scBaseCount_Pipeline/blob/d368c023d64368df65c6f66cfab42f1877d6bef9/scripts/r2/client.py).
- Environment definition: [.env.example](https://github.com/otodreas/scBaseCount_Pipeline/blob/d368c023d64368df65c6f66cfab42f1877d6bef9/.env.example), [pyproject.toml](https://github.com/otodreas/scBaseCount_Pipeline/blob/d368c023d64368df65c6f66cfab42f1877d6bef9/pyproject.toml), and [uv.lock](https://github.com/otodreas/scBaseCount_Pipeline/blob/d368c023d64368df65c6f66cfab42f1877d6bef9/uv.lock).
- Completed manifest: [output/migration/20260506_165832/run.csv](20260506_165832/run.csv).
- Preflight inputs: [tests/datasets_test.csv](https://github.com/otodreas/scBaseCount_Pipeline/blob/d368c023d64368df65c6f66cfab42f1877d6bef9/tests/datasets_test.csv), [tests/datasets_test2.csv](https://github.com/otodreas/scBaseCount_Pipeline/blob/d368c023d64368df65c6f66cfab42f1877d6bef9/tests/datasets_test2.csv), and [tests/datasets_test3.csv](https://github.com/otodreas/scBaseCount_Pipeline/blob/d368c023d64368df65c6f66cfab42f1877d6bef9/tests/datasets_test3.csv).

The tracked preflight manifests record the small tests that preceded the completed transfer:

| Run | Accession | Recorded outcome |
| --- | --- | --- |
| [`output/migration/20260506_154825/run.csv`](20260506_154825/run.csv) | `SRX12708356` | skipped |
| [`output/migration/20260506_155355/run.csv`](20260506_155355/run.csv) | `SRX12708356` | skipped |
| [`output/migration/20260506_155521/run.csv`](20260506_155521/run.csv) | `SRX22996378` | failed because `data/` was not writable |
| [`output/migration/20260506_160248/run.csv`](20260506_160248/run.csv) | `SRX22996378` | uploaded after the local permission problem was resolved |
| [`output/migration/20260506_161016/run.csv`](20260506_161016/run.csv) | `SRX12708356` | skipped |
| [`output/migration/20260506_163944/run.csv`](20260506_163944/run.csv) | `SRX12708356` | skipped |
| [`output/migration/20260506_164701/run.csv`](20260506_164701/run.csv) | `SRX17412822` | dry run |
| [`output/migration/20260506_164720/run.csv`](20260506_164720/run.csv) | `SRX17412822` | uploaded |

## July 30 metadata delta

The July runner compared the 1,816-row source against the 772-row baseline by accession and selected 1,048 rows that were not represented by the baseline. The completed attempt uploaded 1,029 objects and skipped 19 objects already present in R2 with matching MD5 metadata.

Tracked inputs, code, and outputs in snapshot `68c0425`:

- Source: [output/metadata/datasets_v2.csv](https://github.com/otodreas/scBaseCount_Pipeline/blob/68c04253addf9fd01ca6671771a36628b50afab9/output/metadata/datasets_v2.csv), containing 1,816 rows.
- Baseline: [output/metadata/datasets.csv](https://github.com/otodreas/scBaseCount_Pipeline/blob/68c04253addf9fd01ca6671771a36628b50afab9/output/metadata/datasets.csv), containing 772 rows.
- Runner: [pipelines/migrate_gcs_to_r2.py](https://github.com/otodreas/scBaseCount_Pipeline/blob/68c04253addf9fd01ca6671771a36628b50afab9/pipelines/migrate_gcs_to_r2.py).
- Storage code: [scripts/storage/gcs.py](https://github.com/otodreas/scBaseCount_Pipeline/blob/68c04253addf9fd01ca6671771a36628b50afab9/scripts/storage/gcs.py), [scripts/storage/r2.py](https://github.com/otodreas/scBaseCount_Pipeline/blob/68c04253addf9fd01ca6671771a36628b50afab9/scripts/storage/r2.py), and [scripts/storage/transfer.py](https://github.com/otodreas/scBaseCount_Pipeline/blob/68c04253addf9fd01ca6671771a36628b50afab9/scripts/storage/transfer.py).
- Environment definition: [.env.example](https://github.com/otodreas/scBaseCount_Pipeline/blob/68c04253addf9fd01ca6671771a36628b50afab9/.env.example), [pyproject.toml](https://github.com/otodreas/scBaseCount_Pipeline/blob/68c04253addf9fd01ca6671771a36628b50afab9/pyproject.toml), and [uv.lock](https://github.com/otodreas/scBaseCount_Pipeline/blob/68c04253addf9fd01ca6671771a36628b50afab9/uv.lock).

| Attempt | Recorded outcome |
| --- | --- |
| [`output/migration/20260730_145313/run.csv`](20260730_145313/run.csv) | All 1,048 rows failed before any GCS request because the client could not resolve a project identifier.* |
| [`output/migration/20260730_150109/run.csv`](20260730_150109/run.csv) | 1,029 uploaded and 19 skipped. |

*Commit [a35144f](https://github.com/otodreas/scBaseCount_Pipeline/commit/a35144f) supplied a project identifier while retaining anonymous credentials


## Files outside Git

- `.env` contained local credentials and was intentionally ignored.
- `data/` held temporary downloads and was ignored. Files downloaded by the runner were deleted after each upload attempt.
- `logs/migrate_gcs_to_r2.log` was ignored and remains on the server where the transfers ran.
- The mirrored objects reside in R2. Their object keys and source MD5 values are recorded in the completed manifests.
