# FusionCheck

**FusionCheck** is a Perl tool designed to detect and annotate gene fusion events from DNA/RNA sequencing data. It analyzes reads from specified regions in a BAM file, performs local realignment using BLAT to identify potential fusion breakpoints, annotates them with gene information from refGene and repetitive elements from RepeatMasker, and outputs structured fusion event reports.

## Features

- **Flexible Input Regions**: Accepts breakpoint coordinates directly via command line (e.g., `chr20:43961871`) or batch input via a BED file.
- **Accurate Fusion Detection**: Extracts reads containing splice (S) or intron (N) CIGAR operations and realigns them with BLAT to confirm fusion structures; supports both single-end and paired-end data.
- **Gene Annotation**: Utilizes tabix for fast queries of the refGene database to associate breakpoints with gene names and transcriptional directions.
- **Repeat Element Filtering**: Automatically annotates repetitive elements near breakpoints (rmsk) to help exclude false positives.
- **Multi-level Outputs**: Generates detailed alignment information per read, summary of fusion events, lists of discordant read pairs spanning breakpoints, and an optional filtering module (`--filter`) for high-confidence calls.
- **Performance Monitoring**: Automatically logs runtime and memory usage, facilitating resource evaluation.

## Dependencies

Ensure the following software is installed and executable:

- **Perl** (≥5.10) and modules: `Cwd`, `File::Spec`, `Getopt::Long`, `feature`, `List::Util`, `Time::HiRes` (usually core modules)
- **samtools** (≥1.9) – for reading BAM files
- **tabix** (≥1.9) – for indexed queries of compressed refGene and rmsk files
- **BLAT** – for sequence alignment (must support `.2bit` genome files)
- Reference genome: **hg19** in 2bit format (e.g., `hg19.2bit`)
- Annotation files (must be bgzip-compressed and tabix-indexed):
  - refGene.txt.gz and .tbi
  - rmsk.txt.gz and .tbi

## Installation

```bash
git clone https://github.com/Jwei-BI/FusionCheck.git
cd FusionCheck
chmod +x fusionCheck_regionTarget.pl
```

You can either modify the default paths at the top of the script (e.g., `$blat`, `$hg19`, `$refgene`) or specify them via command-line options each run.

## Usage

```bash
perl fusionCheck_regionTarget.pl --bam <bam_file> [options]
```

### Required Arguments

- `--bam` : Path to the input BAM file (must be sorted and indexed).
- **Region specification** (choose one):
  - `--regions "chr1:12345 chr2:67890"` : Space-separated list of breakpoint coordinates.
  - `--bed <bed_file>` : BED file describing two breakpoints per line (format described below).

### Main Options

| Option | Description | Default |
|--------|-------------|---------|
| `--extension <int>` | Bases to extend around each breakpoint when extracting reads | 100 |
| `--min_match <int>` | Minimum number of matching bases required for a fusion event | 20 |
| `--prefix <str>` | Prefix for output files | fusionCheck |
| `--out_dir <dir>` | Output directory | current directory |
| `--fast` | Enable fast mode: directly use hg19.2bit region parameters to restrict BLAT search range | off |
| `--filter` | Enable filtering module to generate `.fusion.filter.xls` | off |
| `--help` | Display help message | off |

### External Tool Path Options

| Option | Description | Default Example |
|--------|-------------|-----------------|
| `--blat <path>` | Path to BLAT executable | /path/to/blat |
| `--samtools <path>` | Path to samtools executable | /software/samtools/1.9/samtools |
| `--tabix <path>` | Path to tabix executable | /software/htslib/1.9/tabix |
| `--hg19 <path>` | Path to hg19.2bit file | /ref/hg19.2bit |
| `--refgene <path>` | Path to compressed refGene.txt.gz | /ref/refGene.txt.sort.gz |
| `--rmsk <path>` | Path to compressed rmsk.txt.gz | /ref/hg19.rmsk.sort.gz |

### Input File Format

#### BED File Format (for `--bed`)

Each line must contain 10 tab-separated fields describing a candidate fusion pair:

| Column | Description |
|--------|-------------|
| 1 | Chromosome of breakpoint 1 |
| 2 | Start position of breakpoint 1 (0-based) |
| 3 | End position of breakpoint 1 |
| 4 | Gene name for breakpoint 1 (optional) |
| 5 | Strand of breakpoint 1 (+/-) |
| 6 | Chromosome of breakpoint 2 |
| 7 | Start position of breakpoint 2 (0-based) |
| 8 | End position of breakpoint 2 |
| 9 | Gene name for breakpoint 2 (optional) |
| 10 | Strand of breakpoint 2 (+/-) |

The program automatically extends each breakpoint by `--extension` bases in both directions to extract reads.

### Examples

#### Basic usage (single breakpoint pair)

```bash
perl fusionCheck_regionTarget.pl \
    --bam tumor.bam \
    --regions "chr20:43961871 chr6:117650683" \
    --extension 100 \
    --min_match 50 \
    --blat /usr/local/bin/blat \
    --hg19 /data/reference/hg19.2bit \
    --refgene /data/annot/refGene.txt.gz \
    --rmsk /data/annot/rmsk.txt.gz \
    --tabix /usr/local/bin/tabix \
    --samtools /usr/local/bin/samtools \
    --prefix sample1 \
    --out_dir ./results
```

#### Using a BED file and enabling filtering

```bash
perl fusionCheck_regionTarget.pl \
    --bam tumor.bam \
    --bed fusions.bed \
    --extension 200 \
    --filter \
    --prefix sample1_filtered
```

#### Fast mode (restrict BLAT search range)

```bash
perl fusionCheck_regionTarget.pl \
    --bam tumor.bam \
    --regions "chr17:57697534" \
    --fast \
    --prefix fast_mode
```

## Output Files

After execution, the following files are created in `--out_dir` (or current directory) with the specified `<prefix>`:

| File | Description |
|------|-------------|
| `<prefix>.read.scan.xls` | Detailed information for each read supporting a fusion, including breakpoint coordinates, genes, and alignment score. |
| `<prefix>.fusion.xls` | Summary per fusion event, listing breakpoints, number of supporting reads, number of spanning read pairs, and average score. |
| `<prefix>.spanning.xls` | List of discordant read pairs that span both breakpoints for each fusion (useful for validation). |
| `<prefix>.fusion.all.xls` | Complete information for all fusion events, including repeat element annotations. |
| `<prefix>.fusion.filter.xls` | (Only if `--filter` is used) High-confidence fusion events after quality filtering. |
| `tmp/` | Temporary directory with intermediate results for each region; can be manually removed after the run. |

### Output Field Descriptions

- **read_name**: Read name, possibly with `/1` or `/2` suffix indicating paired-end.
- **fusion_name**: Format `Gene1-direction1::Gene2-direction2`, where direction indicates the breakpoint relative to the gene’s transcriptional orientation (forward/reversed).
- **breakpoint1/2**: Chromosomal position of the breakpoint, format `chr:pos`.
- **rmsk1/2**: Repeat element information near the breakpoint, format `chr:start-stop,repeat_name,class`; multiple entries separated by semicolons.
- **reads_count**: Number of unique reads supporting the fusion event.
- **spanning_count**: Number of discordant read pairs that span both breakpoints (provides additional support).
- **score_mean**: Average alignment score (matching bases / read length) across supporting reads.

## How It Works

1. **Region Extraction**: Based on user-specified breakpoints, reads overlapping the breakpoints plus flanking regions (extended by `--extension`) are extracted from the BAM file.
2. **Candidate Read Filtering**: Reads containing `N` (intron) or `S` (soft-clip) in their CIGAR strings are retained, as they may indicate fusions or structural variants.
3. **BLAT Realignment**: Candidate read sequences are aligned to the hg19 reference genome using BLAT; output is in PSL format.
4. **Fusion Structure Parsing**:
   - If a single read aligns to two distinct regions (block=2), it is directly considered a potential fusion.
   - For complex alignments (multiple reads each aligning to different regions), combinatorial analysis is performed to identify the most plausible fusion pair.
5. **Breakpoint Gene Annotation**: Tabix queries against the refGene database identify genes overlapping the breakpoints, prioritizing transcripts with complete CDS.
6. **Repeat Annotation**: Breakpoint regions are queried against the rmsk database to mark repetitive elements.
7. **Integration and Filtering**: All events are aggregated, support metrics are calculated, and optional filtering removes low-confidence events or those located within repetitive regions.

## Important Notes

- The BAM file must be sorted and indexed (`.bai`).
- refGene and rmsk files must be compressed with `bgzip` and indexed with `tabix`.
- Default BLAT parameters are optimized for short reads; adjust `$blat_args` in the script if needed.
- Fast mode (`--fast`) restricts BLAT alignment to the specified regions only, greatly speeding up the process but may miss distant fusion partners.

## Author

Junwei Ma
Email: 2254170457@qq.com  
Date: 2024-05-30

## License

This project is licensed under the GPL-3.0 License. See the [LICENSE](LICENSE) file for details.