#!/usr/bin/env python3
"""
collapse_bvbrc_genome_by_accession.py

Purpose
-------
Process a tab-delimited BV-BRC genome metadata table and collapse rows for
non-empty genome.genbank_accessions values when multiple genome.genome_id
values are present for the same accession.

For each non-empty value in the "genome.genbank_accessions" column:
    - If all rows for that accession have the same genome.genome_id, keep all
      rows unchanged.
    - If more than one distinct genome.genome_id is present, keep only the
      single row with the most recent genome.date_modified value.

For rows with an empty genome.genbank_accessions value:
    - Pass the row through unchanged.

A new column, "genome.collapsed_ct", is added:
    - For non-empty genome.genbank_accessions groups, this is the total number
      of original rows for that accession.
    - For rows with empty genome.genbank_accessions, this is set to 1.

Input
-----
A TSV file with a header row.

Output
------
A TSV file with the original columns plus "genome.collapsed_ct".

Notes
-----
- This script assumes the column names are:
      genome.genbank_accessions
      genome.genome_id
      genome.date_modified
- The date comparison uses Python datetime parsing. Several common ISO-like
  formats are supported.
- If multiple rows tie for the newest genome.date_modified within a collapsed
  accession group, the first such row in input order is kept.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from datetime import datetime, timezone
from typing import Dict, List, Tuple


# ---------------------------------------------------------------------------
# Column names used in the input file.
# Change these here if the upstream schema changes.
# ---------------------------------------------------------------------------
ACCESSION_COL = "genome.genbank_accessions"
GENOME_ID_COL = "genome.genome_id"
DATE_MODIFIED_COL = "genome.date_modified"
COLLAPSED_COL = "genome.collapsed_ct"


# ---------------------------------------------------------------------------
# Supported datetime parsing logic.
#
# We try a few common formats first. If the string ends with "Z", we convert it
# to a timezone-aware UTC timestamp. If only a date is provided, it is treated
# as midnight on that date.
#
# If none of the supported formats match, we raise a ValueError so that the
# caller gets a clear message about the offending row.
# ---------------------------------------------------------------------------
def parse_date(date_str: str) -> datetime:
    """
    Parse a genome.date_modified string into a datetime object.

    Supported examples:
        2024-01-31
        2024-01-31 14:22:05
        2024-01-31T14:22:05
        2024-01-31T14:22:05Z
        2024-01-31T14:22:05+00:00

    Parameters
    ----------
    date_str : str
        The input date string from genome.date_modified.

    Returns
    -------
    datetime
        Parsed datetime object.

    Raises
    ------
    ValueError
        If the date string cannot be parsed.
    """
    value = date_str.strip()
    if not value:
        raise ValueError("Empty genome.date_modified value")

    # Handle trailing Z as UTC in a way datetime.fromisoformat accepts.
    if value.endswith("Z"):
        value = value[:-1] + "+00:00"

    # Try Python's ISO parser first.
    try:
        return datetime.fromisoformat(value)
    except ValueError:
        pass

    # Fall back to a few common explicit formats.
    formats = [
        "%Y-%m-%d",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%d %H:%M",
        "%Y/%m/%d",
        "%Y/%m/%d %H:%M:%S",
        "%m/%d/%Y",
        "%m/%d/%Y %H:%M:%S",
    ]
    for fmt in formats:
        try:
            return datetime.strptime(value, fmt)
        except ValueError:
            continue

    raise ValueError(f"Could not parse genome.date_modified value: {date_str!r}")


# ---------------------------------------------------------------------------
# Normalize datetimes so that naive and timezone-aware timestamps can be
# compared safely. Naive datetimes are treated as UTC for comparison purposes.
# ---------------------------------------------------------------------------
def normalize_datetime(dt: datetime) -> datetime:
    """
    Convert a datetime into a comparable form.

    If the datetime is naive (no timezone info), treat it as UTC.
    If it is timezone-aware, convert it to UTC.

    Parameters
    ----------
    dt : datetime
        Input datetime.

    Returns
    -------
    datetime
        Timezone-aware UTC datetime.
    """
    if dt.tzinfo is None:
        return dt.replace(tzinfo=timezone.utc)
    return dt.astimezone(timezone.utc)


# ---------------------------------------------------------------------------
# Read all rows from the TSV file while preserving header order and row order.
# ---------------------------------------------------------------------------
def read_tsv(path: str) -> Tuple[List[str], List[Dict[str, str]]]:
    """
    Read a TSV file into memory.

    Parameters
    ----------
    path : str
        Path to the input TSV file.

    Returns
    -------
    tuple[list[str], list[dict[str, str]]]
        Header fieldnames and list of row dictionaries.
    """
    with open(path, "r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input TSV appears to be missing a header row")
        rows = list(reader)
        return reader.fieldnames, rows


# ---------------------------------------------------------------------------
# Validate that the required columns are present before doing any processing.
# ---------------------------------------------------------------------------
def validate_columns(fieldnames: List[str]) -> None:
    """
    Ensure the required columns exist in the input file.

    Parameters
    ----------
    fieldnames : list[str]
        Header names from the input TSV.

    Raises
    ------
    ValueError
        If one or more required columns are missing.
    """
    required = {ACCESSION_COL, GENOME_ID_COL, DATE_MODIFIED_COL}
    missing = sorted(required - set(fieldnames))
    if missing:
        raise ValueError(
            "Input TSV is missing required column(s): " + ", ".join(missing)
        )


# ---------------------------------------------------------------------------
# Apply the collapsing rules described in the header comment.
#
# We preserve input order for:
#   - rows with empty genome.genbank_accessions
#   - groups that do not need collapsing
#
# For groups that do need collapsing, we keep exactly one row: the newest by
# genome.date_modified. If there is a tie, the first row among tied rows in
# input order is kept.
# ---------------------------------------------------------------------------
def collapse_rows(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    """
    Collapse rows by non-empty genome.genbank_accessions where needed.

    Parameters
    ----------
    rows : list[dict[str, str]]
        Input rows from the TSV.

    Returns
    -------
    list[dict[str, str]]
        Output rows after collapsing and annotating genome.collapsed_ct.
    """
    # Group row indices by non-empty accession while preserving original row
    # order globally.
    accession_to_indices: Dict[str, List[int]] = defaultdict(list)
    for idx, row in enumerate(rows):
        accession = row[ACCESSION_COL].strip()
        if accession:
            accession_to_indices[accession].append(idx)

    # Determine which row indices to keep, and what genome.collapsed_ct value to
    # assign to each kept row.
    keep_indices = set()
    collapsed_counts: Dict[int, int] = {}

    # First handle rows with empty accession values: pass through unchanged and
    # set genome.collapsed_ct = 1.
    for idx, row in enumerate(rows):
        accession = row[ACCESSION_COL].strip()
        if not accession:
            keep_indices.add(idx)
            collapsed_counts[idx] = 1

    # Now process each non-empty accession group.
    for accession, indices in accession_to_indices.items():
        group_rows = [rows[i] for i in indices]
        genome_ids = {r[GENOME_ID_COL].strip() for r in group_rows}
        original_count = len(indices)

        # If the accession group maps to only one genome ID, keep all rows.
        if len(genome_ids) <= 1:
            for idx in indices:
                keep_indices.add(idx)
                collapsed_counts[idx] = original_count
            continue

        # Otherwise, keep only the row with the newest genome.date_modified.
        best_idx = None
        best_dt = None

        for idx in indices:
            row = rows[idx]
            parsed_dt = normalize_datetime(parse_date(row[DATE_MODIFIED_COL]))

            if best_dt is None or parsed_dt > best_dt:
                best_dt = parsed_dt
                best_idx = idx

        assert best_idx is not None, "Internal error: failed to select a best row"
        keep_indices.add(best_idx)
        collapsed_counts[best_idx] = original_count

    # Build the final output list in original input order.
    output_rows: List[Dict[str, str]] = []
    for idx, row in enumerate(rows):
        if idx in keep_indices:
            new_row = dict(row)
            new_row[COLLAPSED_COL] = str(collapsed_counts[idx])
            output_rows.append(new_row)

    return output_rows


# ---------------------------------------------------------------------------
# Write the processed rows back out as TSV.
# ---------------------------------------------------------------------------
def write_tsv(path: str, fieldnames: List[str], rows: List[Dict[str, str]]) -> None:
    """
    Write rows to a TSV file, appending genome.collapsed_ct to the header if needed.

    Parameters
    ----------
    path : str
        Output TSV path.
    fieldnames : list[str]
        Original header names from the input file.
    rows : list[dict[str, str]]
        Processed output rows.
    """
    output_fields = list(fieldnames)
    if COLLAPSED_COL not in output_fields:
        output_fields.append(COLLAPSED_COL)

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=output_fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# Command-line interface.
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments containing input and output paths.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Collapse BV-BRC genome metadata rows by non-empty "
            "'genome.genbank_accessions' when multiple genome IDs are present, "
            "keeping only the row with the newest genome.date_modified."
        )
    )
    parser.add_argument(
        "input_tsv",
        help="Input TSV file (for example: bv-brc-cache/BVBRC_genome.txt)",
    )
    parser.add_argument(
        "output_tsv",
        help="Output TSV file to write",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main program entry point.
# ---------------------------------------------------------------------------
def main() -> int:
    """
    Run the TSV collapsing workflow.

    Returns
    -------
    int
        Exit status code (0 for success, non-zero for failure).
    """
    args = parse_args()

    try:
        fieldnames, rows = read_tsv(args.input_tsv)
        validate_columns(fieldnames)
        output_rows = collapse_rows(rows)
        write_tsv(args.output_tsv, fieldnames, output_rows)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
