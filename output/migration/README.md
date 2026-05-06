# Migration

On 2026-05-06, I migrated all of the datasets we will use (see `[datasets.csv](output/metadata/datasets.csv)`) from [scBaseCount](https://github.com/ArcInstitute/arc-virtual-cell-atlas/blob/main/scBaseCount/README.md) on Google Cloud to Cloudflare R2. There is a 2TB/month limit for free downloads on the Google Cloud hosted scBaseCount dataset, but we have no egress fees on our R2 storage.

The transfer was done on the Lund University bioinformatics server, so that's where the log `migrate_gcs_to_r2.log` can be found.

The final run in `[migration](output/migration)` is the one to look at. The ones before were tests, transferring over a few small `h5ad` files.