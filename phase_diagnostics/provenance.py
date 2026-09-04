"""Small, deterministic run provenance records."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shlex
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping


SAFE_ENV_PREFIXES = ("HK_", "OMP_", "MKL_", "CUDA_", "CONDA_", "PYTHONHASHSEED")


def _source_files(repo: Path) -> list[Path]:
    """Return the small source surface that determines diagnostic outputs."""
    paths: set[Path] = set()
    for pattern in (
        "phase_diagnostics/*.py",
        "scripts/run_phase*.py",
        "scripts/run_phase*.sh",
        "scripts/generate_phase*.py",
        "tests/test_phase_diagnostics.py",
    ):
        paths.update(path for path in repo.glob(pattern) if path.is_file())
    for relative in (
        "Makefile",
        "blind.c",
        "hickit.h",
        "run_blind_p9016_minimal.c",
        "test_blind_estep_scores.c",
        "eval/evaluate_p9016_baseline.py",
        "pytest.ini",
    ):
        path = repo / relative
        if path.is_file():
            paths.add(path)
    return sorted(paths)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _git(repo: Path, *args: str) -> str:
    return subprocess.check_output(
        ("git", *args), cwd=repo, text=True, stderr=subprocess.DEVNULL
    ).rstrip("\r\n")


def collect_provenance(
    repo: Path, command: Iterable[str], seed: int, config: Mapping[str, object],
    inputs: Iterable[Path] = (),
) -> dict[str, object]:
    status = _git(repo, "status", "--porcelain")
    source_files = _source_files(repo)
    git_diff = subprocess.check_output(
        ("git", "diff", "--binary", "HEAD", "--", *[str(path.relative_to(repo)) for path in source_files]),
        cwd=repo,
    )
    safe_env = {
        key: value for key, value in sorted(os.environ.items())
        if any(key == prefix or key.startswith(prefix) for prefix in SAFE_ENV_PREFIXES)
    }
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python": sys.version.replace("\n", " "),
        "conda_default_env": os.environ.get("CONDA_DEFAULT_ENV", ""),
        "git_commit": _git(repo, "rev-parse", "HEAD"),
        "git_dirty": bool(status),
        "git_status": status.splitlines(),
        "git_diff_sha256": hashlib.sha256(git_diff).hexdigest(),
        "source_sha256": {
            str(path.relative_to(repo)): file_sha256(path) for path in source_files
        },
        "command": shlex.join(list(command)),
        "seed": int(seed),
        "environment": safe_env,
        "configuration": dict(config),
        "input_sha256": {str(path): file_sha256(path) for path in inputs if path.is_file()},
    }


def write_provenance(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
