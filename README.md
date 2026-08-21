# checkstrain

Small Nextflow (DSL2) pipeline that, for each bacterial strain in the input:

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

## Running on Snellius

Same pattern as this project's other pipelines: a `snellius` profile (SLURM executor,
partition `genoa`, Singularity) is included in `nextflow.config`. A pre-built Bakta
database already exists at `/projects/prjs1784/db/bakta_db`, so pass it with
`--bakta_database` to skip the ~60GB download on every run.

```bash
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --job-name=checkstrain
#SBATCH --output=checkstrain_%j.out

eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run /projects/prjs1784/pipelines/checkstrain/main.nf \
    -profile snellius \
    --input /projects/prjs1784/salt/data/samplesheet.csv \
    --outdir /projects/prjs1784/salt/results_checkstrain \
    --bakta_database /projects/prjs1784/db/bakta_db \
    -resume
```

`--osmo_db` is left unset here, so `OSMOTOOL_DOWNLOAD_DB` fetches `osmo_refdb` (latest
release) once per run automatically; pass `--osmo_db /path/to/refdb/<release>` once you
have a copy on disk to skip that download too, the same way `--bakta_database` does.

## Output

Per sample, under `outdir/`:

- `bakta_bakta/<id>.tsv`, `.gff3`, `.faa`, ... — full Bakta annotation.
- `osmotool_annotate/<id>.gene_counts.tsv`, `.systems.tsv`, `.gene_coordinates.tsv` — full osmotool profile of the same genome.
- `checkstrain_select_genes_of_interest/<id>.genes_of_interest.tsv` — Bakta's annotation rows restricted to genes on osmotool's panel, with `osmotool_family`, `osmotool_raw_count`, `osmotool_norm_col`, `osmotool_norm_value` appended.
- `checkstrain_compare_bakta_osmotool/<id>.bakta_osmotool_comparison.tsv` — one row per gene name seen by either method, with a `status` of `concordant` (both found it), `bakta_only`, or `osmotool_only`, plus each side's loci (`bakta_loci`/`osmotool_loci`, semicolon-joined if a gene has multiple hits on one side). Matched by gene name, not coordinates — Bakta renumbers contigs during annotation (`contig_5`, ...) while osmotool's `gene_coordinates.tsv` keeps the input FASTA's original contig names, so the two aren't coordinate-comparable.

## Why annotate twice (Bakta + osmotool)?

osmotool's `annotate` mode calls its own ORFs (Prodigal) and only reports the osmoadaptation gene panel — it doesn't carry Bakta's richer per-gene fields (locus tag, product description, DbXrefs, EC/COG, ...). Running Bakta gives the full annotation; running osmotool on the same FASTA gives the authoritative panel of gene family names it profiles (its `*.gene_counts.tsv` always lists all families, with count 0 where absent — see `osmotool.quantifier.KNOWN_FAMILIES`). The selection step joins the two by gene symbol so the final table has Bakta's annotation detail but only for genes osmotool actually tracks.
