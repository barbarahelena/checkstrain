# checkstrain

Aim: to compare the output of bakta and osmotool for a gene panel of interest.
For each bacterial strain in the input:

1. **Downloads** the genome assembly from NCBI ([`datasets`](https://www.ncbi.nlm.nih.gov/datasets/) CLI), given a name + NCBI assembly accession.
2. **Annotates** the genome with [Bakta](https://github.com/oschwengers/bakta).
3. **Selects the genes of interest**: cross-references Bakta's annotation against the osmoadaptation gene panel that [osmotool](../osmotool) profiles (run on the same genome via `osmotool annotate`), producing a filtered table of just those genes with Bakta's annotation fields plus osmotool's per-family counts.
4. **Compares the two independent gene calls**: Bakta's filtered gene calls vs. osmotool's own coordinate-level hits (`gene_coordinates.tsv`), matched by gene name, to flag genes each method found that the other didn't.

## Usage

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    -profile docker
```

`samplesheet.csv`:

```csv
name,accession
Escherichia coli K-12 MG1655,GCF_000005845.2
Bacillus subtilis 168,GCF_000009045.1
```

- `name` — free-text strain/species name (used to derive the sample id).
- `accession` — the NCBI **assembly** accession (`GCF_...`/`GCA_...`), not a taxon or BioSample ID.

Smoke-test with a small built-in samplesheet:

```bash
nextflow run main.nf -profile test,docker --outdir results
```

## Key params

| Param               | Default    | Description                                                                 |
|---------------------|------------|-------------------------------------------------------------------------------|
| `--input`            | `null`     | Samplesheet CSV (`name,accession`). Required.                                 |
| `--outdir`           | `null`     | Output directory. Required.                                                   |
| `--bakta_database`   | `null`     | Path to a pre-built Bakta database; downloaded automatically if not set.      |
| `--osmo_db`          | `null`     | Path to a pre-built `osmo_refdb` release; downloaded automatically if not set.|
| `--osmo_db_release`  | `latest`   | `osmo_refdb` release to fetch when `--osmo_db` is not set.                    |

## Output

Per sample, under `outdir/`:

- `bakta_bakta/<id>.tsv`, `.gff3`, `.faa`, ... — full Bakta annotation.
- `osmotool_annotate/<id>.gene_counts.tsv`, `.systems.tsv`, `.gene_coordinates.tsv` — full osmotool profile of the same genome.
- `checkstrain_select_genes_of_interest/<id>.genes_of_interest.tsv` — Bakta's annotation rows restricted to genes on osmotool's panel, with `osmotool_family`, `osmotool_raw_count`, `osmotool_norm_col`, `osmotool_norm_value` appended.
- `checkstrain_compare_bakta_osmotool/<id>.bakta_osmotool_comparison.tsv` — one row per gene name seen by either method, with a `status` of `concordant` (both found it), `bakta_only`, or `osmotool_only`, plus each side's loci (`bakta_loci`/`osmotool_loci`, semicolon-joined if a gene has multiple hits on one side). Matched by gene name, not coordinates — Bakta renumbers contigs during annotation (`contig_5`, ...) while osmotool's `gene_coordinates.tsv` keeps the input FASTA's original contig names, so the two aren't coordinate-comparable.
