from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class JsonlSnapshot:
    """One accession-keyed row from a contexts JSONL file."""

    accession: str
    raw_line: str
    record: dict[str, Any]


@dataclass
class CompareResult:
    """Summary of differences between two contexts JSONL snapshots."""

    deleted: list[str] = field(default_factory=list)
    added: list[str] = field(default_factory=list)
    identical: list[str] = field(default_factory=list)
    changed: list[str] = field(default_factory=list)
    changed_field_counts: Counter[str] = field(default_factory=Counter)


def _parse_source(spec: str) -> tuple[str, Path | None]:
    """Return `(git_ref, path)` when spec is a git ref, else `(None, file_path)`."""
    if ":" in spec and not Path(spec).exists():
        ref, _, path = spec.partition(":")
        if ref and path:
            return ref, Path(path)
    path = Path(spec)
    if path.exists():
        return None, path
    return spec, Path("output/context/contexts.jsonl")


def _load_jsonl_text(text: str) -> dict[str, JsonlSnapshot]:
    """Parse JSONL text into accession-keyed snapshots."""
    rows: dict[str, JsonlSnapshot] = {}
    for line_no, line in enumerate(text.splitlines(), start=1):
        raw_line = line.rstrip("\n")
        if not raw_line:
            continue
        record = json.loads(raw_line)
        accession = record.get("accession")
        if not accession:
            raise ValueError(f"Line {line_no} is missing accession")
        if accession in rows:
            raise ValueError(f"Duplicate accession {accession!r}")
        rows[accession] = JsonlSnapshot(accession=accession, raw_line=raw_line, record=record)
    return rows


def load_jsonl_source(spec: str) -> dict[str, JsonlSnapshot]:
    """Load contexts JSONL from a file path or git ref such as `HEAD:output/context/contexts.jsonl`."""
    git_ref, path = _parse_source(spec)
    if git_ref is not None:
        repo_path = path.as_posix() if path is not None else "output/context/contexts.jsonl"
        text = subprocess.check_output(["git", "show", f"{git_ref}:{repo_path}"], text=True)
        return _load_jsonl_text(text)
    if path is None or not path.exists():
        raise FileNotFoundError(f"Could not resolve source {spec!r}")
    return _load_jsonl_text(path.read_text())


def _normalize_record(record: dict[str, Any]) -> dict[str, Any]:
    """Return a copy with list order normalized for semantic comparison."""
    normalized = json.loads(json.dumps(record))
    run_accessions = normalized.get("runAccessions")
    if isinstance(run_accessions, list):
        normalized["runAccessions"] = sorted(run_accessions)
    return normalized


def _diff_paths(old: Any, new: Any, prefix: str = "") -> list[str]:
    """Return dotted paths where two JSON values differ."""
    if type(old) is not type(new):
        return [prefix or "<root>"]
    if isinstance(old, dict):
        diffs: list[str] = []
        for key in sorted(set(old) | set(new)):
            child_prefix = f"{prefix}.{key}" if prefix else key
            if key not in old or key not in new:
                diffs.append(child_prefix)
            else:
                diffs.extend(_diff_paths(old[key], new[key], child_prefix))
        return diffs
    if isinstance(old, list):
        return [] if old == new else [prefix or "<root>"]
    return [] if old == new else [prefix or "<root>"]


def compare_snapshots(
    old_rows: dict[str, JsonlSnapshot],
    new_rows: dict[str, JsonlSnapshot],
    *,
    semantic: bool = False,
) -> CompareResult:
    """Compare two accession-keyed snapshots and classify rows as deleted, added, identical, or changed."""
    result = CompareResult()
    old_keys = set(old_rows)
    new_keys = set(new_rows)

    result.deleted = sorted(old_keys - new_keys)
    result.added = sorted(new_keys - old_keys)

    for accession in sorted(old_keys & new_keys):
        old_row = old_rows[accession]
        new_row = new_rows[accession]
        if not semantic:
            if old_row.raw_line == new_row.raw_line:
                result.identical.append(accession)
            else:
                result.changed.append(accession)
                top_level = {path.split(".")[0] for path in _diff_paths(old_row.record, new_row.record)}
                for field_name in top_level:
                    result.changed_field_counts[field_name] += 1
            continue

        if _normalize_record(old_row.record) == _normalize_record(new_row.record):
            result.identical.append(accession)
        else:
            result.changed.append(accession)
            top_level = {
                path.split(".")[0]
                for path in _diff_paths(
                    _normalize_record(old_row.record),
                    _normalize_record(new_row.record),
                )
            }
            for field_name in top_level:
                result.changed_field_counts[field_name] += 1

    return result


def _print_result(result: CompareResult, *, list_changed: bool, verbose: bool) -> None:
    """Print a human-readable comparison summary."""
    print(f"deleted:   {len(result.deleted)}")
    print(f"added:     {len(result.added)}")
    print(f"identical: {len(result.identical)}")
    print(f"changed:   {len(result.changed)}")

    if result.deleted:
        print("\nDeleted accessions:")
        for accession in result.deleted:
            print(f"  {accession}")

    if result.added:
        print("\nAdded accessions:")
        for accession in result.added:
            print(f"  {accession}")

    if verbose and result.changed_field_counts:
        print("\nTop-level fields changed among modified rows:")
        for field_name, count in result.changed_field_counts.most_common():
            print(f"  {field_name}: {count}")

    if list_changed and result.changed:
        print("\nChanged accessions:")
        for accession in result.changed:
            print(f"  {accession}")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare two contexts JSONL files or a git ref against the working tree file.",
    )
    parser.add_argument(
        "old",
        help="Baseline source: file path, git ref (HEAD), or git object path (HEAD:output/context/contexts.jsonl).",
    )
    parser.add_argument(
        "new",
        nargs="?",
        default="output/context/contexts.jsonl",
        help="New source to compare against. Defaults to output/context/contexts.jsonl.",
    )
    parser.add_argument(
        "--semantic",
        action="store_true",
        help="Treat runAccessions order differences as identical.",
    )
    parser.add_argument(
        "--list-changed",
        action="store_true",
        help="Print changed accession IDs.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print top-level field counts for changed rows.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Compare two contexts JSONL sources and exit non-zero unless only deletions occurred."""
    args = _build_parser().parse_args(argv)
    old_rows = load_jsonl_source(args.old)
    new_rows = load_jsonl_source(args.new)
    result = compare_snapshots(old_rows, new_rows, semantic=args.semantic)
    _print_result(result, list_changed=args.list_changed, verbose=args.verbose)

    only_deletions = not result.added and not result.changed
    if only_deletions:
        print("\nOK: only deletions")
        return 0

    print("\nNOT OK: file differs beyond deletions")
    return 1


if __name__ == "__main__":
    sys.exit(main())
