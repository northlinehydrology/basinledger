"""The Pleasant Lake model (Wisconsin Central Sands), a published MODFLOW 6 model Basin Ledger didn't write, into the
cache: the second way in's first real test model (docs/design.md, "Extending an existing MODFLOW 6 model").

    uv run python examples/central-sands/fetch.py       # 16 MB download, 212 MB unzipped, into the cache

Source: Fienen, M.N., Haserodt, M.J., and Leaf, A.T., 2021, MODFLOW models used to simulate groundwater flow in the
Wisconsin Central Sands Study Area, 2012-2018: U.S. Geological Survey data release, doi:10.5066/P9BVFSGJ.
"""

from __future__ import annotations

import hashlib
import sys
import urllib.request
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from basinledger.config import cache_dir  # noqa: E402

ITEM = "https://www.sciencebase.gov/catalog/file/get/60623131d34e227131f99776?f=__disk__"
FILES = {
    "PLEASANT_LAKE_MODEL.zip": ("b4%2F6c%2F44%2Fb46c446356354f2ba79dbf75765063e21fea4815",
                                "bd65b78e8904226356bbfe8e2b2fec9bbca7b4ab4ec2289b8645c28c64337a07"),
}


def folder() -> Path:
    return cache_dir() / "central-sands"


def workspace() -> Path:
    """The simulation as published: the parent model and its LGR inset."""
    return folder() / "pleasant" / "final_model_files"


def fetch() -> Path:
    out = folder()
    out.mkdir(parents=True, exist_ok=True)
    for name, (key, sha) in FILES.items():
        path = out / name
        if not path.exists():
            req = urllib.request.Request(ITEM + key, headers={"User-Agent": "basinledger-fetch"})
            with urllib.request.urlopen(req) as r, open(path, "wb") as f:
                f.write(r.read())
        if hashlib.sha256(path.read_bytes()).hexdigest() != sha:
            raise ValueError(f"{path}: not the published file (sha256 differs); delete it and fetch again")
    if not workspace().is_dir():
        with zipfile.ZipFile(out / "PLEASANT_LAKE_MODEL.zip") as z:
            z.extractall(out / "pleasant")
    return workspace()


if __name__ == "__main__":
    print(fetch())
