"""Shared helpers: config loading, paths, label-level mapping.

Both the Python tool scripts and the metrics script import from here so the
granularity levels and ground-truth conventions stay identical everywhere.
"""
from __future__ import annotations

import os
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPO_ROOT / "config" / "config.yaml"

LEVELS = ("coarse", "medium", "fine")


def load_config(path: str | os.PathLike | None = None) -> dict:
    with open(path or CONFIG_PATH) as fh:
        return yaml.safe_load(fh)


def repo_path(*parts: str) -> Path:
    return REPO_ROOT.joinpath(*parts)


def map_to_level(raw_labels, level: str, mapping: dict) -> list:
    """Collapse raw author labels to one granularity level.

    `mapping` is a dict raw_label -> level_label, built per dataset from the
    author annotation. Unmapped labels become None so metrics can mask them
    (an unmapped truth label is not scored, never counted as wrong).
    """
    if level not in LEVELS:
        raise ValueError(f"level must be one of {LEVELS}, got {level!r}")
    return [mapping.get(lbl) for lbl in raw_labels]


def prediction_file(cfg: dict, dataset: str, tool: str) -> Path:
    """Canonical output path every tool writes to: cell_id, predicted_label."""
    return repo_path(cfg["paths"]["predictions"]) / f"{dataset}__{tool}.csv"
