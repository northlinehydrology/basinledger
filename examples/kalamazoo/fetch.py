"""The Kalamazoo River basin model (Michigan, HUC 04050003): the steady MODFLOW-NWT inset model cut from USGS's Lake
Michigan Basin regional model, the second way in's Michigan test model once converted to MODFLOW 6 (docs/design.md,
"Extending an existing MODFLOW 6 model").

ScienceBase serves the archive's model.zip (1.04 GB, all five inset models) only after a captcha, so it is
downloaded by hand from https://doi.org/10.5066/F76D5R5V (Attached Files, model.zip) into the cache folder
(~/.cache/basinledger/kalamazoo, or $BASINLEDGER_CACHE/kalamazoo) or ~/Downloads. This script checks it against the
published MD5, extracts only the Kalamazoo model (``nwt/``) and converts it to MODFLOW 6 with mf5to6 (``mf6/``):

    uv run python examples/kalamazoo/fetch.py [path/to/model.zip]

Source: Haserodt, M.J., Kauffman, L.J., and Feinstein, D.T., 2018, MODFLOW-NWT inset models from the regional Lake
Michigan Basin model in support of groundwater age calculations for glacial aquifers: U.S. Geological Survey data
release, doi:10.5066/F76D5R5V; report Feinstein and others, 2018, SIR 2018-5038, doi:10.3133/sir20185038.
"""

from __future__ import annotations

import hashlib
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from basinledger.config import cache_dir  # noqa: E402

MD5 = "928ab6bec42adf5a14922bdba15c87f7"  # ScienceBase's checksum for model.zip


def folder() -> Path:
    return cache_dir() / "kalamazoo"


def nwt_folder() -> Path:
    """The Kalamazoo inset model as published (MODFLOW-NWT)."""
    return folder() / "nwt"


def find_zip(given: str | None = None) -> Path:
    for p in ([Path(given).expanduser()] if given else []) + [folder() / "model.zip", Path.home() / "Downloads" / "model.zip"]:
        if p.exists():
            return p
    raise FileNotFoundError("model.zip not found: download it from https://doi.org/10.5066/F76D5R5V (Attached "
                            f"Files) into {folder()} or ~/Downloads")


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def extract(path: Path) -> Path:
    """The Kalamazoo model's own folder only (``kala/``, not the seeded-well example), into ``nwt_folder()``."""
    if md5(path) != MD5:
        raise ValueError(f"{path}: not the published model.zip (MD5 differs)")
    out = nwt_folder()
    with zipfile.ZipFile(path) as z:
        members = [m for m in z.namelist() if "kala/" in m.lower().replace("\\", "/") and not m.endswith("/")]
        if not members:
            raise ValueError(f"{path}: no kala/ folder in it")
        for m in members:
            rel = m.replace("\\", "/").split("kala/", 1)[1] if "kala/" in m else m.replace("\\", "/").lower().split("kala/", 1)[1]
            dest = (out / rel).resolve()
            if out.resolve() not in dest.parents:
                raise ValueError(f"{m}: outside the folder")
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(z.read(m))
    return out


def mf6_folder() -> Path:
    """The Kalamazoo model converted to MODFLOW 6 by mf5to6."""
    return folder() / "mf6"


def convert(nwt: Path | None = None) -> Path:
    """MODFLOW-NWT to MODFLOW 6 with USGS's mf5to6 (``uv run get-modflow :flopy --repo modflow6 --subset mf5to6``),
    into ``mf6_folder()``. The converter reads a copy of the name file with Unix line endings, without MNWI (MNW2's
    output reporting, which mf5to6 doesn't convert) and the three files only MNWI writes, and with the other output
    files marked REPLACE. The model's input is untouched. MNW2's wells become MAW wells."""
    import shutil
    import subprocess

    from basinledger.coupling.modflow6 import find_library
    nwt = nwt or nwt_folder()
    exe = find_library().parent / "mf5to6"
    if not exe.exists():
        exe = Path(shutil.which("mf5to6") or "")
    if not exe.is_file():
        raise FileNotFoundError("mf5to6 not found: uv run get-modflow :flopy --repo modflow6 --subset mf5to6")
    lines = (nwt / "KALA4x-MPfix2.nam").read_text().replace("\r", "").splitlines()
    skip = ("MNWI ",)
    drop = ("mnw2fix2.wel", ".byn", ".qsm")
    keep = [ln + (" REPLACE" if ln.startswith("DATA") else "") for ln in lines
            if ln.strip() and not ln.startswith(skip) and not any(d in ln for d in drop)]
    (nwt / "convert.nam").write_text("\n".join(keep) + "\n")
    before = {p.name for p in nwt.iterdir()}
    r = subprocess.run([str(exe), "convert.nam", "kala"], cwd=nwt, capture_output=True, text=True)
    if "Conversion successful" not in r.stdout:
        raise RuntimeError(f"mf5to6 failed:\n{r.stdout[-2000:]}")
    out = mf6_folder()
    out.mkdir(parents=True, exist_ok=True)
    for p in nwt.iterdir():
        if p.name not in before:
            shutil.move(str(p), out / p.name)
    return out


if __name__ == "__main__":
    nwt = nwt_folder() if (nwt_folder() / "KALA4x-MPfix2.nam").exists() else \
        extract(find_zip(sys.argv[1] if len(sys.argv) > 1 else None))
    print(convert(nwt))
