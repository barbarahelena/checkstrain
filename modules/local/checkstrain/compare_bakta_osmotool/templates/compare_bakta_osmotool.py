#!/usr/bin/env python3
"""Cross-check Bakta's osmoadaptation gene calls against osmotool's own
coordinate-level hits for the same genome, to surface genes where the two
independent gene-calling methods (Bakta's annotation pipeline vs. osmotool's
own Prodigal + DIAMOND search) disagree.

Matching is done by gene name, not genomic coordinates: Bakta renumbers
contigs during annotation (e.g. "contig_5"), while osmotool's
gene_coordinates.tsv keeps the original input FASTA's contig names (e.g.
"JPJG01000046.1") -- the two tools' contig identifiers don't correspond, so
coordinate overlap isn't a usable join key here.

Bakta's rows are grouped by their already-matched `osmotool_family` column
(see select_genes_of_interest.py) rather than the raw Bakta `Gene` symbol,
since that column is already normalised against this same panel (paralog
suffixes like "_2" stripped) -- reusing it avoids re-deriving that logic.

Nextflow's `template` mechanism rewrites backslash escapes (e.g. \\n, \\t) in
this file before Python sees it -- avoid embedding \\n in single-line string
literals; use print()/multi-line strings with real newlines instead.
"""

import csv
import platform
from collections import defaultdict

# Interpolated by Nextflow from the process' inputs/task context (see
# select_genes_of_interest.py for why `prefix` is computed here rather than
# via a `def prefix = ...` in the process' script block).
genes_of_interest_tsv = "$genes_of_interest_tsv"
gene_coordinates_tsv  = "$gene_coordinates_tsv"
prefix                = "${task.ext.prefix ?: meta.id}"


def bakta_loci_by_gene(path):
    by_gene = defaultdict(list)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = (row.get("osmotool_family") or "").strip()
            if not gene:
                continue
            locus = "{}:{}-{}({})/{}".format(
                row.get("Sequence Id", "-"),
                row.get("Start", "-"),
                row.get("Stop", "-"),
                row.get("Strand", "-"),
                row.get("Locus Tag", "-"),
            )
            by_gene[gene].append(locus)
    return by_gene


def osmotool_loci_by_gene(path):
    by_gene = defaultdict(list)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = (row.get("gene") or "").strip()
            if not gene:
                continue
            locus = "{}:{}-{}({})/{}".format(
                row.get("contig", "-"),
                row.get("start", "-"),
                row.get("end", "-"),
                row.get("strand", "-"),
                row.get("protein_id", "-"),
            )
            by_gene[gene].append(locus)
    return by_gene


bakta_by_gene = bakta_loci_by_gene(genes_of_interest_tsv)
osmotool_by_gene = osmotool_loci_by_gene(gene_coordinates_tsv)

all_genes = sorted(set(bakta_by_gene) | set(osmotool_by_gene))

out_fields = ["gene", "status", "bakta_hits", "bakta_loci", "osmotool_hits", "osmotool_loci"]

n_concordant = 0
n_bakta_only = 0
n_osmotool_only = 0

with open(f"{prefix}.bakta_osmotool_comparison.tsv", "w", newline="") as out_f:
    writer = csv.DictWriter(out_f, fieldnames=out_fields, delimiter="\t")
    writer.writeheader()
    for gene in all_genes:
        b = bakta_by_gene.get(gene, [])
        o = osmotool_by_gene.get(gene, [])
        if b and o:
            status = "concordant"
            n_concordant += 1
        elif b:
            status = "bakta_only"
            n_bakta_only += 1
        else:
            status = "osmotool_only"
            n_osmotool_only += 1
        writer.writerow(
            {
                "gene": gene,
                "status": status,
                "bakta_hits": len(b),
                "bakta_loci": ";".join(b) or "-",
                "osmotool_hits": len(o),
                "osmotool_loci": ";".join(o) or "-",
            }
        )

print(
    f"{prefix}: {n_concordant} concordant, {n_bakta_only} bakta-only, "
    f"{n_osmotool_only} osmotool-only ({len(all_genes)} total genes)"
)

with open("versions.yml", "w") as fh:
    print('"$task.process":', file=fh)
    print(f"    python: {platform.python_version()}", file=fh)
