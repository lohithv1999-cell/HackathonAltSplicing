# HackathonAltSplicing

![image](https://github.com/user-attachments/assets/3f94eb52-b412-4ec2-96fc-67f5e3fb34d3)

## Quick Start

This is a fork of [haessar/HackathonAltSplicing](https://github.com/haessar/HackathonAltSplicing),
extending Workflow 2 with a quality-control framework, configurable library
strandedness and annotation handling, and a split-half stability check. These changes
are intended for merge upstream.

Clone the repository and change directory:

```
git clone https://github.com/lohithv1999-cell/HackathonAltSplicing.git
cd HackathonAltSplicing
```

### Isoform detection pipeline

The Snakemake workflow in isoform_detection/ (Workflow 2 in the figure above) assigns every read in an aligned BAM to an assembled transcript, tagging it with an integer IC value, and reports a set of quality-control metrics alongside. It runs on bulk RNA-seq or given a barcode to cell-type mapping, on single-cell data split by cell type. Slurm bash scripts for running on an HPC system are in the scripts directory.

#### 1. Environment

Install conda if you do not have it (Miniforge), then:
```
conda create -n altsplice -c conda-forge -c bioconda --file isoform_detection/requirements.txt --yes
conda activate altsplice
```

kallisto is pinned separately and installed by Snakemake itself, from isoform_detection/env/kallisto.yml, provided you pass --use-conda when you run. You do not need to install it yourself.

Check the tools are on the path before going further:

##### On Windows

This will not work on Windows directly. The pipeline depends on stringtie, samtools, gffread, gffcompare and kallisto, a
ll of which come from bioconda, and bioconda builds for Linux and macOS only. 
Installing from requirements.txt in a Windows conda prompt fails to solve the environment rather than reporting the platform as the reason, 
which makes it easy to misread as a dependency conflict.

Use the Windows Subsystem for Linux instead. In PowerShell as administrator:

```
powershell
wsl --install
```

That installs WSL2 with Ubuntu. Restart, set a username and password when prompted, then install Miniforge inside the Linux environment and follow the instructions above unchanged:

```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

If you use VS Code, install the WSL extension and open the repository with `code .` from inside the WSL terminal, or use **Connect to WSL** from the command palette. That gives you the editor on Windows with the terminal, Python interpreter and conda environment all running inside Linux, which is usually the most comfortable way to work.

One thing worth getting right from the start: keep the repository and your data inside the Linux filesystem, under ~/, rather than on a Windows drive reached through /mnt/c/. Cross-filesystem access in WSL is slow, and this pipeline reads and writes large BAM files throughout, so the difference is substantial rather than marginal.

Alternatively, if you have access to a Linux server or HPC system, run it there.

#### 2. Inputs

You need three things, none of which are in this repository:

| Input	Notes | Notes |
|---|---|
| **Aligned BAM(s)** |Coordinate sorted. For 10x data, the `possorted_genome_bam.bam` from Cell Ranger. |
| **Reference genome FASTA** | Transcript sequences are extracted from this; required even when no annotation is supplied. |
| **Reference annotation GFF3** |	Optional |

Supplying an annotation lets StringTie use known transcript models to inform assembly. 
Without one, transcripts are reconstructed de novo from read coverage alone, which the pipeline supports 
and which still recovers a comparable transcript set on the data tested here. 
An annotation is worth supplying where a good one exists, particularly for cell types with few reads, 
since it provides structure the reads alone may not support.

For single-cell data you also need a CSV mapping cell barcodes to cell types. 
Examples are in `analysis/` - cc_barcode_Mira_1.csv and similar - and the format is two columns:

```
cell_barcode,Cluster
AAACCCAAGAATCCCT-1,Neuron 1
AAACCCACATAGGTTC-1,Stem C
```

The cell types in this CSV determine what the pipeline treats as a sample: 
it reads the unique values and produces one output per cell type. 
If no CSV is given the pipeline runs in bulk mode, treating the BAM as a single sample.

#### 3. Configuration

Everything is set in isoform_detection/config.yaml. 
The committed copy contains paths from the machine it was last run on, 
so it will not work unmodified. At minimum change:

```
yaml
run_identity: "my_run"          # names the output directory; use a fresh name per run

mixed_bams:                      # your BAM(s)
  - "/path/to/possorted_genome_bam.bam"

metadata_csv:                    # barcode-to-cell-type CSV; leave empty [] for bulk mode
  - "/path/to/barcodes.csv"

genome: "/path/to/genome.fa"         # required
gff: "/path/to/annotation.gff3"  # optional; use "" to assemble de novo
```

Three settings are worth understanding rather than accepting:

library_type = "fr", "rf", or "" if unstranded. This matters more than it looks. 
StringTie takes the strand of a spliced alignment from the XS:A: tag, 
and discards spliced reads that lack one unless the library type is declared. 
Cell Ranger does not write XS, so leaving this empty on 10x data silently throws away every junction-spanning read 
and produces an assembly with no spliced transcripts at all. If unsure which of fr and rf applies, 
run both on one sample and compare each against your reference annotation with gffcompare; 
the wrong orientation matches almost nothing.

read_mode = "single" or "paired", fragment_length, fragment_sd - kallisto cannot estimate fragment length from single-end data, 
so these are supplied directly. They feed the model that produces the assignment probabilities.

zw_threshold - the probability above which a read is assigned, default 0.8. 
This is a real choice rather than a formality: retention falls by around nineteen percentage points between cutoffs of 0.5 and 0.95 on our data.

#### 4. Run

From inside isoform_detection/:

```
cd isoform_detection
snakemake -n --cores 8            # dry run: check the job list looks right
snakemake --cores 8 --use-conda   # run
```

--use-conda is required, since kallisto is installed into its own environment. Raise --cores to whatever the machine has.

A full cohort takes hours, so give it time before assuming something has stalled. 
If a run is interrupted, Snakemake leaves a lock behind; clear it and continue with:

```
snakemake --unlock
snakemake --cores 8 --use-conda --rerun-incomplete
```

#### 5. Outputs

Everything is written to ../results/<run_identity>/, that is, alongside the repository rather than inside isoform_detection/:

```
results/<run_identity>/
├── fastq/                        extracted reads
├── kallisto/                     index, pseudoalignments, abundance estimates
├── qc/
│   ├── completeness/             per cell type
│   ├── ic0/                      why unassigned reads failed
│   ├── sensitivity/              retention across ZW cutoffs
│   ├── splithalf/                stability, when enabled
│   │   ├── bams/                 the split halves
│   │   ├── gffcmp/               gffcompare output per cell type (.tracking, .stats, .loci, .annotated.gtf)
│   │   ├── unstable_loci/        loci whose counts disagree most between halves
│   │   ├── <cell>_stability.tsv  per cell type
│   │   └── ALL_stability.tsv     the collated table
│   ├── QC_Dashboard.html         all metrics on one page; open in a browser
│   └── summary_statistics.tsv    headline figures for the run
├── stringtie/                    assembled transcripts per cell type (GTF, FASTA)
└── wp3_add_ic_tag/               the tagged BAMs - the main output
```

Start with QC_Dashboard.html. An exported chart can be traced back to the run that produced it.

#### 6. Optional QC

Two checks are controlled from the config:

```
yaml
run_ic0_profiling: true    # on by default; adds a pass over each BAM
run_splithalf: false       # off by default; runs assembly and pseudoalignment twice
```

IC:i:0 profiling classifies each unassigned read as either scored but below threshold or never pseudo-aligned. 
It costs one additional pass over each BAM, which is modest against the rest of the run, so it is left on.

Split-half stability is off because it is expensive. It splits every cell type in two and runs assembly and pseudoalignment on each half. 
Treat it as a pre-flight check rather than something to run routinely; 
enable it, read the report, and use it to decide whether to commit to a full tagged run. 
It requires an assembly containing spliced transcripts; 
on an assembly with none it reports that it cannot be assessed rather than returning a misleading zero.

### JBrowse 2 with plugin and custom tracks

For running JBrowse 2 with the plugin developed for Workflow 1, 
ensure Node.js is installed (e.g. `sudo apt install nodejs npm`) and follow the guide in `jbrowse-plugin-bedfeaturecoloring/README.md`. 
To serve a custom `config.json` for loading in tracks for analysis, in a new terminal:

```
cd analysis/jbrowse_env/
npx serve . --cors -p 3001
```

and navigate to `http://localhost:3000/?config=http://localhost:3001/config.json` in a web browser.

## Contributing

### 🔖 Issue Labelling

I've created issues that reflect the work packages that I outlined in the introductory slides (for some of them it made sense to split them in two). 
Work package #1 is signified by `(WP1)`, etc. I've used labels to designate which "theme" they belong to, 
whether they have a significant MARS component, and whether they are `Coding` or `Research` heavy. 
When determining what you want to work on, you can filter the labels that appeal to you or are a suitable match to your skillset 
(e.g. if you have no interest in using MARS, you could filter out the `MARS` label).

### 🔀 Branching Strategy

As there will be lots of us working feverishly in the same repo, I recommend we use **feature branches** for all development work rather than `main`.  

- Create a new branch from `main` for each new feature or fix (perhaps to work on a GitHub issue):  
  ```bash
  git checkout -b your-feature-name
  ```
- Once you're happy with changes, open a pull request to merge into `main` and someone can review them.

- This will help prevent conflicts, support collaboration, and maintain a clean commit history.