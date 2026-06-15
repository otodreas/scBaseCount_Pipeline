"""Generate zero-shot scGPT cell embeddings from a raw h5ad and write a .npz artifact."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import scanpy as sc
from scgpt.tasks import embed_data


def run_embed(
    h5ad_path: Path,
    model_dir: Path,
    out_path: Path,
    gene_col: str = "gene_symbols",
    max_length: int = 1200,
    batch_size: int = 64,
    device: str = "cpu",
) -> Path:
    """Load raw h5ad, compute scGPT embeddings, and write obs_names + X to out_path."""
    adata = sc.read_h5ad(h5ad_path)
    adata.obs_names_make_unique()

    embed_data(
        adata,
        model_dir=model_dir,
        gene_col=gene_col,
        max_length=max_length,
        batch_size=batch_size,
        device=device,
        use_fast_transformer=False,
        return_new_adata=False,
    )

    if "X_scGPT" not in adata.obsm:
        raise RuntimeError("embed_data did not populate adata.obsm['X_scGPT']")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        out_path,
        obs_names=np.asarray(adata.obs_names, dtype=str),
        X=np.asarray(adata.obsm["X_scGPT"], dtype=np.float32),
    )
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate scGPT cell embeddings from a raw h5ad.")
    parser.add_argument("--h5ad", type=Path, required=True, help="Path to raw h5ad file")
    parser.add_argument("--model-dir", type=Path, required=True, help="scGPT checkpoint directory")
    parser.add_argument("--out", type=Path, required=True, help="Output .npz path (obs_names + X)")
    parser.add_argument("--gene-col", default="gene_symbols", help="adata.var column with gene symbols")
    parser.add_argument("--max-length", type=int, default=1200)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--device", default="cpu", help="cpu or cuda")
    args = parser.parse_args()

    out = run_embed(
        h5ad_path=args.h5ad,
        model_dir=args.model_dir,
        out_path=args.out,
        gene_col=args.gene_col,
        max_length=args.max_length,
        batch_size=args.batch_size,
        device=args.device,
    )
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
