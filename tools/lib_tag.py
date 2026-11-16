"""The wheel platform tag for MODFLOW 6's shared library (libmf6) from USGS's release, read from the library itself,
or an error if it couldn't load on a machine of that tag as it is:

    python tools/lib_tag.py ~/.local/share/flopy/bin/libmf6.dylib      # macosx_14_0_arm64, win_amd64, manylinux_2_35_x86_64...

- Windows: a 64-bit x86 DLL importing only Windows' own DLLs (read from its import table).
- Linux: an x86-64 or ARM library needing only glibc's own libraries; tagged manylinux with the newest glibc symbol
  version it uses.
- macOS: arm64 or x86_64, linked only to libSystem, tagged with the oldest macOS it was built for.
"""

from __future__ import annotations

import re
import struct
import subprocess
import sys
from pathlib import Path

# DLLs every Windows has (System32), besides the api-ms-win-* sets of the Universal C runtime
WINDOWS_DLLS = {"kernel32.dll", "kernelbase.dll", "ntdll.dll", "msvcrt.dll", "ucrtbase.dll", "user32.dll", "gdi32.dll",
                "advapi32.dll", "shell32.dll", "shlwapi.dll", "ole32.dll", "oleaut32.dll", "comdlg32.dll", "ws2_32.dll",
                "bcrypt.dll", "crypt32.dll", "secur32.dll", "rpcrt4.dll", "version.dll", "psapi.dll", "dbghelp.dll",
                "winmm.dll", "iphlpapi.dll", "shcore.dll", "userenv.dll", "imagehlp.dll", "dnsapi.dll",
                "netapi32.dll", "wininet.dll", "setupapi.dll", "cfgmgr32.dll"}
GLIBC = re.compile(r"^(libc|libm|libdl|libpthread|librt|libutil|ld-linux[\w.-]*)\.so(\.\d+)*$")


def pe_imports(data: bytes, pe: int) -> set[str]:
    """The DLLs a PE32+ file imports, from its import table."""
    nsect, = struct.unpack_from("<H", data, pe + 6)
    opt_size, = struct.unpack_from("<H", data, pe + 20)
    opt = pe + 24
    if struct.unpack_from("<H", data, opt)[0] != 0x20B:
        raise SystemExit("not a PE32+ (64-bit) file")
    rva, size = struct.unpack_from("<II", data, opt + 112 + 8)          # data directory 1: imports
    sections = [struct.unpack_from("<IIII", data, opt + opt_size + 40 * i + 8) for i in range(nsect)]

    def off(r):
        for vsize, vaddr, rawsize, rawptr in sections:
            if vaddr <= r < vaddr + max(vsize, rawsize):
                return rawptr + r - vaddr
        raise SystemExit(f"address {r:#x} in no section")

    names, d = set(), off(rva) if size else None
    while d is not None:
        name_rva = struct.unpack_from("<I", data, d + 12)[0]
        if name_rva == 0:
            break
        o = off(name_rva)
        names.add(data[o:data.index(b"\0", o)].decode().lower())
        d += 20
    return names


def windows(data: bytes) -> str:
    pe = struct.unpack_from("<I", data, 0x3C)[0]
    if data[pe:pe + 4] != b"PE\0\0":
        raise SystemExit("not a PE file")
    machine = struct.unpack_from("<H", data, pe + 4)[0]
    if machine != 0x8664:
        raise SystemExit(f"PE machine {machine:#x}, not x86-64")
    other = sorted(n for n in pe_imports(data, pe) if not (n in WINDOWS_DLLS or n.startswith("api-ms-win-")))
    if other:
        raise SystemExit(f"needs DLLs a plain Windows doesn't have: {other}")
    return "win_amd64"


def linux(lib: Path, data: bytes) -> str:
    machine = struct.unpack_from("<H", data, 18)[0]
    arch = {62: "x86_64", 183: "aarch64"}.get(machine)
    if arch is None:
        raise SystemExit(f"ELF machine {machine}")
    dyn = subprocess.run(["readelf", "-d", str(lib)], check=True, capture_output=True, text=True).stdout
    needed = re.findall(r"\(NEEDED\)\s+Shared library: \[([^\]]+)\]", dyn)
    other = [n for n in needed if not GLIBC.match(n)]
    if other:
        raise SystemExit(f"needs libraries beyond glibc's: {other}")
    vers = subprocess.run(["readelf", "-V", str(lib)], check=True, capture_output=True, text=True).stdout
    newest = max(((int(a), int(b)) for a, b in re.findall(r"GLIBC_(\d+)\.(\d+)", vers)), default=(2, 17))
    newest = max(newest, (2, 17))
    return f"manylinux_{newest[0]}_{newest[1]}_{arch}"


def macos(lib: Path) -> str:
    own = subprocess.run(["otool", "-D", str(lib)], check=True, capture_output=True, text=True).stdout.splitlines()[1:]
    libs = subprocess.run(["otool", "-L", str(lib)], check=True, capture_output=True, text=True).stdout.splitlines()[1:]
    names = [line.split()[0] for line in libs if line.strip()]
    other = [n for n in names if n not in [o.strip() for o in own] and not n.startswith("/usr/lib/libSystem")]
    if other:
        raise SystemExit(f"linked to libraries a plain Mac doesn't have: {other}")
    load = subprocess.run(["otool", "-l", str(lib)], check=True, capture_output=True, text=True).stdout
    m = re.search(r"cmd LC_BUILD_VERSION.*?minos (\d+)\.(\d+)", load, re.S)
    if not m:
        raise SystemExit("no LC_BUILD_VERSION: can't tell the oldest macOS it runs on")
    arch = subprocess.run(["lipo", "-archs", str(lib)], check=True, capture_output=True, text=True).stdout.split()
    if len(arch) != 1 or arch[0] not in ("arm64", "x86_64"):
        raise SystemExit(f"architectures {arch}")
    major, minor = int(m[1]), int(m[2])
    return f"macosx_{major}_{minor if major < 11 else 0}_{arch[0]}"


def tag(lib: Path) -> str:
    data = lib.read_bytes()
    if data[:2] == b"MZ":
        return windows(data)
    if data[:4] == b"\x7fELF":
        return linux(lib, data)
    if data[:4] in (b"\xcf\xfa\xed\xfe", b"\xca\xfe\xba\xbe"):
        return macos(lib)
    raise SystemExit(f"{lib}: not a library this knows")


if __name__ == "__main__":
    print(tag(Path(sys.argv[1])))
