import scanpy as sc
from cytetype import CyteType
from dotenv import load_dotenv
from shared.repo import REPO_ROOT

load_dotenv()

INPUT = REPO_ROOT / "output/atlas/v2/post/production/atlas_v2_post.h5ad"

OUTPUT = REPO_ROOT / "output/atlas/v2/post/production/cytetype/atlas_v2_post_cytetype.h5ad"

group_key = "leiden_atlas"  # leiden_uncorrected

atlas = sc.read(INPUT, backed="r")

annotator = CyteType(
    adata=atlas,
    group_key=group_key,
    # rank_key
)

atlas = annotator.run(
    study_context="mouse",
)

INPUT.unlink()
atlas.write_h5ad(OUTPUT, compression="gzip")
