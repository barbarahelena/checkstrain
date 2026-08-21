#!/usr/bin/env python3
"""Filter a Bakta annotation table down to the osmoadaptation gene panel that
osmotool profiles, joining in osmotool's per-family counts for the same genome.

osmotool's *.gene_counts.tsv always lists every family in its panel (count 0
for families with no hit -- see osmotool.quantifier.count_hits), so it is
used here as the authoritative "genes of interest" list rather than
hardcoding a copy of the panel that could drift out of sync with the
osmo_refdb release actually used upstream.

Nextflow's `template` mechanism rewrites backslash escapes (e.g. \\n, \\t) in
this file before Python sees it, since the whole file is interpolated as a
Groovy string -- so avoid embedding \\n in single-line string literals; use
print()/multi-line strings with real newlines instead.
"""

import csv
import platform
import re

# Interpolated by Nextflow from the process' inputs/task context. `prefix` is
# computed here rather than via a `def prefix = ...` in the process' script
# block: template files are rendered against the process' own binding (inputs,
# task, workflow, params), which does NOT include locally-scoped `def`
# variables declared in the script block -- those only exist in that block's
# Groovy closure, not the template's binding, and referencing one from here
# raises "No such property" at render time.
bakta_tsv       = "$bakta_tsv"
gene_counts_tsv = "$gene_counts_tsv"
prefix          = "${task.ext.prefix ?: meta.id}"

# Bakta sometimes suffixes a gene symbol with "_2", "_3", ... when a family
# has multiple paralogous hits in the same genome (e.g. "opuAA_2"); strip
# that before falling back to a second lookup against the osmotool panel.
PARALOG_SUFFIX_RE = re.compile(r"_[0-9]+\$")


def load_osmotool_panel(path):
    """Return {lowercase family name: {family, raw_count, norm_col, norm_value}}."""
    panel = {}
    with open(path, newline="") as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            header = line.strip().split("\t")
            break
        if header is None:
            return panel
        norm_col = header[2] if len(header) > 2 else "norm"
        reader = csv.DictReader(fh, fieldnames=header, delimiter="\t")
        for row in reader:
            gene = row["gene"]
            panel[gene.lower()] = {
                "family": gene,
                "raw_count": row.get("raw_count", ""),
                "norm_col": norm_col,
                "norm_value": row.get(norm_col, ""),
            }
    return panel


def bakta_rows(path):
    """Yield (fieldnames, DictReader) for every annotated feature in a Bakta TSV.

    Bakta's *.tsv starts with a few '#'-prefixed comment lines (tool/version,
    command line, ...) before the real, also '#'-prefixed, column header
    ("#Sequence Id\tType\t...\tDbXrefs"); everything after that is data.
    """
    with open(path, newline="") as fh:
        lines = fh.readlines()
    header_idx = next(
        i for i, line in enumerate(lines) if line.startswith("#Sequence Id")
    )
    fieldnames = lines[header_idx].lstrip("#").strip().split("\t")
    reader = csv.DictReader(lines[header_idx + 1 :], fieldnames=fieldnames, delimiter="\t")
    return fieldnames, reader


panel = load_osmotool_panel(gene_counts_tsv)
fieldnames, reader = bakta_rows(bakta_tsv)

out_fields = fieldnames + [
    "osmotool_family",
    "osmotool_raw_count",
    "osmotool_norm_col",
    "osmotool_norm_value",
]

n_hits = 0
with open(f"{prefix}.genes_of_interest.tsv", "w", newline="") as out_f:
    writer = csv.DictWriter(out_f, fieldnames=out_fields, delimiter="\t")
    writer.writeheader()
    for row in reader:
        gene = (row.get("Gene") or "").strip()
        if not gene:
            continue
        key = gene.lower()
        hit = panel.get(key) or panel.get(PARALOG_SUFFIX_RE.sub("", key))
        if hit is None:
            continue
        n_hits += 1
        writer.writerow(
            {
                **row,
                "osmotool_family": hit["family"],
                "osmotool_raw_count": hit["raw_count"],
                "osmotool_norm_col": hit["norm_col"],
                "osmotool_norm_value": hit["norm_value"],
            }
        )

print(f"{prefix}: {n_hits} Bakta-annotated gene(s) matched the osmotool panel ({len(panel)} families)")

with open("versions.yml", "w") as fh:
    print('"$task.process":', file=fh)
    print(f"    python: {platform.python_version()}", file=fh)
