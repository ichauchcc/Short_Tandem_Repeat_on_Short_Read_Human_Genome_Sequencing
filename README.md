# Short Tandem Repeat Analysis on Short-Read Human Genome Sequencing

**Repository:** `Short_Tandem_Repeat_on_Short_Read_Human_Genome_Sequencing`  
**Platform:** NCI Gadi (PBS Pro)  
**Primary tools:** ExpansionHunter (EH) - ExpansionHunter-v5.0.0-linux_x86_64; ExpansionHunter Denovo (EHdn) - ExpansionHunterDenovo-v0.9.0; samtools - 1.23.1, Python3  
**Focus:** Human germline short tandem repeat (STR) expansions, with emphasis on inherited cardiac disease and sudden cardiac death (SCD).

> **Research use only.** STR calls are candidates, not molecular diagnoses. Confirm clinically relevant expansions with an orthogonal assay and interpret them in the context of phenotype, inheritance, assay limitations and gene–disease validity.

## Contents

1. [Overview](#1-overview)
2. [Biological mechanisms](#2-biological-mechanisms)
3. [Inheritance and interpretation](#3-inheritance-and-interpretation)
4. [Cardiac disease and sudden cardiac death](#4-cardiac-disease-and-sudden-cardiac-death)
5. [Detection strategies](#5-detection-strategies)
6. [Requirements and Gadi environment](#6-requirements-and-gadi-environment)
7. [Directory structure](#7-directory-structure)
8. [Step 1 — input quality control](#8-step-1--input-quality-control)
9. [Step 2 — known-locus genotyping with ExpansionHunter](#9-step-2--known-locus-genotyping-with-expansionhunter)
10. [Step 3 — genome-wide discovery with EHdn](#10-step-3--genome-wide-discovery-with-ehdn)
11. [Step 4 — prioritisation and validation](#11-step-4--prioritisation-and-validation)
12. [Reporting and reproducibility](#12-reporting-and-reproducibility)
13. [Limitations](#13-limitations)
14. [References and resources](#14-references-and-resources)

## 1. Overview

STRs are tandem arrays of short nucleotide motifs, often 1–6 bp long; disease-associated tandem repeat expansions can also involve longer repeat units. Expanded alleles may be missed by conventional SNV/indel pipelines because reads map ambiguously and expansions can exceed read or fragment length.

This workflow performs **complementary analyses** of aligned, short-read whole-genome sequencing (WGS) data:

| Method | Objective | Input | Main output | Key limitation |
|---|---|---|---|---|
| **ExpansionHunter** | Genotype *specified* STR loci | Indexed BAM/CRAM, matched reference FASTA, variant catalog | Per-locus repeat estimates, VCF, JSON, supporting alignments | Very long/complex expansions and repeat interruptions may be poorly resolved |
| **ExpansionHunter Denovo** | Detect *candidate* expansions without a predefined locus list | Short-read BAM/CRAM; preferably a comparable cohort | Per-sample STR profiles and cohort outlier/case-control statistics | Does not provide definitive allele sizes or precise locus boundaries |

**Important input distinction:** despite the historical description “WGS FASTA files”, EH and EHdn operate on **aligned sequencing reads (BAM/CRAM)**, not a consensus FASTA. The reference FASTA is a separate required input for EH. Start after the upstream alignment pipeline has produced `sample.final.bam` and its index.

### Research questions

- Are established disease-associated repeat loci expanded in an individual?
- Are previously uncharacterised repeat motifs/loci unusually abundant in an SCD cohort?
- Do putative expansions fit an established cardiac phenotype, mechanism and inheritance pattern?
- Can candidates be confirmed by repeat-primed PCR, long-read sequencing or another validated assay?

## 2. Biological mechanisms

Repeat expansions cause disease through several distinct mechanisms; **repeat size alone does not establish pathogenicity**.

| Mechanism | Molecular effect | Illustrative loci | Clinical interpretation |
|---|---|---|---|
| **Toxic RNA gain of function** | Expanded RNA sequesters RNA-binding proteins, causing abnormal splicing | `DMPK` CTG (DM1), `CNBP` CCTG (DM2) | Cardiac conduction disease and arrhythmia are particularly relevant in myotonic dystrophy |
| **Transcriptional silencing / loss of function** | Expansion-associated chromatin changes reduce gene expression | `FXN` intronic GAA (Friedreich ataxia), `FMR1` full mutation CGG (fragile X syndrome) | For `FXN`, reduced frataxin causes mitochondrial dysfunction and cardiomyopathy; mechanism and inheritance are locus-specific |
| **Toxic protein gain of function** | Expanded coding repeat produces an abnormal protein (often polyglutamine) | `HTT` CAG and several spinocerebellar ataxia genes | Primarily neurological phenotypes; cardiac association requires locus-specific evidence |
| **Repeat-associated non-AUG (RAN) translation / RNA toxicity** | Expanded repeats produce abnormal RNA and/or repeat-derived peptides | `C9orf72` GGGGCC, `FMR1` premutation-associated conditions | Mechanisms vary; not automatically cardiac disease genes |
| **Altered gene regulation or other locus-specific effects** | Changes in transcript abundance, splicing or chromatin | Some noncoding repeat disorders | Requires locus-specific experimental and clinical support |

Additional modifiers include repeat interruptions, somatic mosaicism, tissue-specific expansion, age-dependent penetrance and anticipation. Blood-derived repeat lengths may not represent cardiac tissue. For `FXN`, the predominant pathogenic mechanism is **reduced frataxin expression**, not a primary toxic-RNA mechanism; most affected individuals have biallelic pathogenic variants (commonly two expanded GAA alleles, or an expansion in trans with another pathogenic `FXN` allele).

## 3. Inheritance and interpretation

| Pattern | Examples | What must be established |
|---|---|---|
| **Autosomal dominant (AD)** | `DMPK` (DM1), `CNBP` (DM2), `HTT` | One disease-range expanded allele may be sufficient; account for age-dependent penetrance and anticipation |
| **Autosomal recessive (AR)** | `FXN` (Friedreich ataxia) | **Biallelic** pathogenic variation in trans; one expanded allele usually indicates carrier status rather than a diagnosis |
| **X-linked** | `FMR1`, `AFF2` | Interpret sex chromosome complement, repeat range, methylation and mosaicism; hemizygosity and heterozygosity are not interchangeable |
| **Complex / locus-specific** | `RFC1` repeat-associated ataxia | Motif identity and phase matter; two repeat alleles are not necessarily pathogenic unless they have the relevant configuration |

**Do not infer inheritance from the gene name or chromosome alone.** Disease mechanism, pathogenic repeat motif, size range, penetrance, parental origin and sex chromosome status must be interpreted using current locus-specific clinical resources. Male/female labels in job scripts are metadata; they do not replace genotype-based sex-chromosome assessment. For X-linked loci, verify the ploidy and output conventions of the chosen software version and catalog.

### Repeat-size categories

Categories such as normal, intermediate, premutation and full mutation are **locus-specific**, and sometimes depend on repeat interruptions and methylation. Avoid universal thresholds. For example, intermediate-range `DMPK` alleles should **not** be labelled disease-causing solely because they are more frequent in cases than controls.

## 4. Cardiac disease and sudden cardiac death

| Locus / condition | Inheritance | Established cardiac involvement | Relevance to SCD | Interpretation caution |
|---|---|---|---|---|
| **`DMPK` CTG — myotonic dystrophy type 1** | AD | Progressive PR/QRS prolongation, AV block, atrial and ventricular arrhythmias; sometimes cardiomyopathy | **Directly relevant:** conduction disease and arrhythmias can cause sudden death; presentations without obvious neuromuscular signs have been reported | Repeat length alone is an unreliable individual predictor of SCD; assess ECG, family history and confirm expansion |
| **`CNBP` CCTG — myotonic dystrophy type 2** | AD | Conduction abnormalities and arrhythmias; cardiac involvement generally less prominent than DM1 | Relevant differential diagnosis, with less direct SCD evidence than DM1 | Complex motif and repeat structure may require specialist testing |
| **`FXN` GAA — Friedreich ataxia** | AR | Hypertrophic cardiomyopathy, myocardial fibrosis, arrhythmias and heart failure | Cardiac complications can contribute to premature death, including sudden death | **One heterozygous expansion is generally insufficient**; evaluate second allele and neurological phenotype |
| **Other neurological repeat disorders** | Varies | Some have reported cardiac features | Consider only with appropriate phenotype and evidence | Do not presume that every repeat-associated disorder explains unexplained SCD |

**Molecular autopsy prioritisation:** For unexplained SCD, especially in a structurally normal heart, prioritise `DMPK` when conduction-system disease, arrhythmia, myotonia, cataracts, family history or unexplained sudden deaths are present. For cardiomyopathy (especially hypertrophic or mitochondrial-like pathology), consider `FXN` if systemic features or biallelic variants support Friedreich ataxia. Evaluate `CNBP` when the clinical presentation is compatible with DM2. Repeat expansions should be considered **alongside**, not instead of, conventional cardiomyopathy/channelopathy variants, structural variants, mitochondrial variants and non-genetic causes.

A 2026 *European Heart Journal* report described `DMPK` 3′UTR expansions in families with unexplained SCD without overt neuromuscular or conduction-system features; this supports including repeat analysis in selected molecular-autopsy workflows, but does not imply that every rare or intermediate `DMPK` allele is causal.

## 5. Detection strategies

### Evidence types in short-read WGS

- **Spanning reads:** contain the complete repeat and both flanks; informative for shorter alleles.
- **Flanking reads:** anchor to unique sequence and extend into the repeat.
- **In-repeat reads (IRRs):** composed predominantly of repeat sequence; support long expansions but are hard to place uniquely.
- **Paired-end information:** a uniquely mapped mate can help localise a repeat-derived read.
- **Read-depth / cohort outlier signal:** can identify expanded motifs but does not independently prove allele length or locus.

**EH** uses a targeted sequence-graph approach to estimate repeat alleles at specified loci. **EHdn** detects repeat-rich reads genome-wide and compares STR profiles across samples. They answer different questions; EHdn is **not** simply a more comprehensive EH genotyper.

### Suggested validation

| Question | Assay |
|---|---|
| Is the expected motif expanded? | Repeat-primed PCR (RP-PCR/TP-PCR), where validated for the locus |
| What is the allele size? | Flanking PCR/capillary electrophoresis for amplifiable alleles; Southern blot or long-read sequencing for larger expansions as appropriate |
| What are the motif interruptions and haplotypes? | Targeted or whole-genome long-read sequencing |
| Is silencing present? | Methylation or expression studies where mechanistically appropriate (e.g. `FMR1`, `FXN`) |
| Does the candidate segregate with disease? | Family genotyping with an assay suitable for the expansion |

## 6. Requirements and Gadi environment

### Gadi infrastructure

NCI Gadi uses **PBS Pro** and environment modules. The following commands are intended for Gadi, not an Artemis-specific filesystem. Software module names and versions **must be checked on the live system**: availability of preinstalled EH/EHdn modules has **not been verified**.

```bash
module avail 2>&1 | grep -Ei 'samtools|python|expansion|bcftools|apptainer|singularity'
module spider samtools 2>&1 || true
command -v qsub
command -v samtools || true
command -v ExpansionHunter || true
command -v ExpansionHunterDenovo || true
```

Load the modules actually present, e.g. `module load samtools/<AVAILABLE_VERSION>` and `module load python3/<AVAILABLE_VERSION>`. Do **not** paste `<AVAILABLE_VERSION>` literally. If EH or EHdn is unavailable, install a versioned binary or approved container in your project software directory, following the upstream documentation and NCI policies. Pin the exact executable version, catalog commit and container digest if used.

**Required:** indexed coordinate-sorted BAM (or supported CRAM), matching indexed GRCh38 reference FASTA, versioned EH variant catalog, reference-compatible chromosome naming, access to Gadi project storage, and appropriate human-genomics data permissions. PCR-free WGS is preferred for EH. EHdn benefits from a sufficiently large, technically matched comparator cohort.

### Configure project paths

Create `config.sh` (do not commit patient-specific identifiers or protected paths):

```bash
#!/usr/bin/env bash
export PROJECT='YOUR_NCI_PROJECT'
export ROOT="/g/data/${PROJECT}/STR_analysis"
export BAM_DIR="/scratch/${PROJECT}/Fastq-to-VCF/Final_bam"
export REF='/g/data/YOUR_PROJECT/reference/GRCh38.fa'
export CATALOG='/g/data/YOUR_PROJECT/catalogs/eh_grch38_catalog.json'
export EH_BIN='/g/data/YOUR_PROJECT/software/ExpansionHunter'
export EHDN_BIN='/g/data/YOUR_PROJECT/software/ExpansionHunterDenovo'
mkdir -p "${ROOT}"/{eh,ehdn,logs,metadata}
```

Replace placeholders with actual accessible paths. `/scratch` is temporary working storage subject to NCI retention/cleanup policies; preserve validated results and metadata in approved persistent project storage. `GRCh38.fa` must match the exact alignment reference, not merely the assembly name.

## 7. Directory structure

```text
Short_Tandem_Repeat_on_Short_Read_Human_Genome_Sequencing/
├── README.md
├── config.sh                  # local, gitignored
├── jobs/
│   ├── expansionhunter.pbs
│   └── ehdn_profile.pbs
├── catalogs/                  # provenance/links; versioned catalog if licensing permits
├── metadata/                  # de-identified cohort manifest
└── results/                   # do not commit identifiable patient outputs
```

Recommended `.gitignore`:

```gitignore
config.sh
*.bam
*.bai
*.cram
*.crai
*.fastq.gz
results/
logs/
metadata/*patient*
```

## 8. Step 1 — input quality control

```bash
source ./config.sh
sample='XXX'
bam="${BAM_DIR}/${sample}.final.bam"

[[ -s "$bam" ]] || { echo "Missing BAM: $bam" >&2; exit 1; }
[[ -s "$REF" && -s "${REF}.fai" ]] || { echo 'Missing indexed reference' >&2; exit 1; }
samtools quickcheck -v "$bam" || exit 1
[[ -s "${bam}.bai" || -s "${bam%.bam}.bai" ]] || samtools index "$bam"
samtools view -H "$bam" | grep '^@SQ' | head
samtools faidx "$REF" | head -n 2
```

Check reference sequence names/lengths against BAM `@SQ`, coverage, read length, library insert size, duplicate rate, mapping quality, PCR-free status, sex chromosome coverage and sample identity. Do not assume `sample.final.bam` is suitable solely because it exists.

## 9. Step 2 — known-locus genotyping with ExpansionHunter

### Create `jobs/expansionhunter.pbs`

```bash
#!/bin/bash
#PBS -N EH
#PBS -l ncpus=4
#PBS -l mem=16GB
#PBS -l walltime=08:00:00
#PBS -l storage=gdata/YOUR_NCI_PROJECT+scratch/YOUR_NCI_PROJECT
#PBS -j oe

set -euo pipefail
cd "$PBS_O_WORKDIR"
source ./config.sh
: "${sample:?Supply sample with qsub -v sample=XXX}"

bam="${BAM_DIR}/${sample}.final.bam"
[[ -s "$bam" ]] || { echo "Missing $bam" >&2; exit 1; }
[[ -s "$CATALOG" ]] || { echo "Missing $CATALOG" >&2; exit 1; }
[[ -s "$REF" ]] || { echo "Missing $REF" >&2; exit 1; }
[[ -x "$EH_BIN" ]] || { echo "Missing executable $EH_BIN" >&2; exit 1; }
mkdir -p "${ROOT}/eh/${sample}"

"$EH_BIN" \
  --reads "$bam" \
  --reference "$REF" \
  --variant-catalog "$CATALOG" \
  --output-prefix "${ROOT}/eh/${sample}/${sample}" \
  --threads "${PBS_NCPUS:-4}"
```

**Before submission:** replace `YOUR_NCI_PROJECT` in the PBS resource line with the actual project ID. PBS directives are interpreted **before** `config.sh` is sourced, so shell variables cannot be used in that directive. Adjust memory and walltime after benchmarking. The EH command is based on the upstream CLI; verify supported options with `ExpansionHunter --help` for your pinned version.

### Submit a single sample

```bash
qsub -P YOUR_NCI_PROJECT -v sample=XXX -N EH_XXX jobs/expansionhunter.pbs
```

### Submit all final BAMs

```bash
source ./config.sh
for bam in "${BAM_DIR}"/*.final.bam; do
    [[ -e "$bam" ]] || continue
    sample=$(basename "$bam" .final.bam)
    echo qsub -P "$PROJECT" -v "sample=$sample" -N "EH_$sample" jobs/expansionhunter.pbs
done
```

The loop **prints commands** for inspection; remove `echo` to submit. Validate sample identifiers against `[A-Za-z0-9_.-]+` before using them in scheduler/job paths. Keep `sex` in a separately curated sample manifest if needed; the standard EH CLI does not require a `sex` argument for the command shown here.

### Interpret EH results

Inspect per-locus VCF and JSON for estimated repeat allele sizes, uncertainty, read evidence, genotype/ploidy and warnings. Compare with a **locus-specific** disease threshold, phenotype and expected inheritance. Record catalog build, repeat motif, reference assembly and tool version. Review supporting alignments with REViewer when supported by the locus and catalog.

## 10. Step 3 — genome-wide discovery with EHdn

### Create `jobs/ehdn_profile.pbs`

```bash
#!/bin/bash
#PBS -N EHdn
#PBS -l ncpus=4
#PBS -l mem=24GB
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/YOUR_NCI_PROJECT+scratch/YOUR_NCI_PROJECT
#PBS -j oe

set -euo pipefail
cd "$PBS_O_WORKDIR"
source ./config.sh
: "${sample:?Supply sample with qsub -v sample=XXX}"

bam="${BAM_DIR}/${sample}.final.bam"
[[ -s "$bam" ]] || { echo "Missing $bam" >&2; exit 1; }
[[ -x "$EHDN_BIN" ]] || { echo "Missing executable $EHDN_BIN" >&2; exit 1; }
mkdir -p "${ROOT}/ehdn/${sample}"

"$EHDN_BIN" profile \
  --reads "$bam" \
  --output-prefix "${ROOT}/ehdn/${sample}/${sample}"
```

Verify the `profile` syntax with the pinned EHdn binary's `--help` and upstream documentation. EHdn profiles use repeat-containing read evidence and do not directly produce clinically interpretable diploid repeat genotypes.

### Submit a single sample

```bash
qsub -P YOUR_NCI_PROJECT -v sample=XXX -N EHdn_XXX jobs/ehdn_profile.pbs
```

### Generate cohort submissions

```bash
source ./config.sh
for bam in "${BAM_DIR}"/*.final.bam; do
    [[ -e "$bam" ]] || continue
    sample=$(basename "$bam" .final.bam)
    echo qsub -P "$PROJECT" -v "sample=$sample" -N "EHdn_$sample" jobs/ehdn_profile.pbs
done
```

Review the printed jobs, then remove `echo` to execute. For large cohorts, PBS job arrays may be preferable, with explicit concurrency limits and a versioned sample manifest.

### Cohort analysis: essential second stage

**Generating individual EHdn profiles is not the end of the analysis.** Merge profiles across the cohort and run either:

1. **Outlier analysis:** recommended for heterogeneous unexplained SCD cases, where different cases may have different causal loci.
2. **Case-control analysis:** appropriate when multiple cases are expected to share an expanded locus and controls are adequately matched.

Use the upstream EHdn scripts and documentation for the **exact commands associated with the pinned release**: [case-control](https://github.com/Illumina/ExpansionHunterDenovo/blob/master/documentation/03_Case_control_quickstart.md) and [outlier](https://github.com/Illumina/ExpansionHunterDenovo/blob/master/documentation/04_Outlier_quickstart.md). Do not compare raw EHdn counts from mismatched read lengths, library protocols, coverage or alignment pipelines without careful normalisation and sensitivity analysis. EHdn candidate positions may be approximate and must be localised and genotyped using additional methods.

## 11. Step 4 — prioritisation and validation

A defensible prioritisation framework for SCD cohorts:

1. **Technical validity:** adequate coverage, uniquely anchored reads, stable signal under filtering; inspect BAM and repeat context.
2. **Locus identification:** map EHdn signal to a genomic region, verify motif, check assembly/strand, and determine whether an EH catalog can genotype it.
3. **Population evidence:** use technically matched controls, public repeat resources where appropriate, ancestry-aware frequency estimates and quality filters.
4. **Clinical concordance:** evaluate autopsy findings (structural vs autopsy-negative), ECG, conduction disease, arrhythmia, cardiomyopathy, neuromuscular symptoms and family history.
5. **Mechanistic plausibility:** established gene–disease association, pathogenic repeat motif and range, expected AD/AR/X-linked genotype, transcript effect and other alleles.
6. **Segregation:** test available relatives; phase both alleles in AR disorders.
7. **Orthogonal confirmation:** repeat-primed PCR, capillary electrophoresis, Southern blot or long-read sequencing as appropriate.
8. **Alternative explanations:** review SNVs/indels, CNVs/SVs, mitochondrial variants, other repeat loci, toxicology and non-genetic causes.
9. **Clinical interpretation:** discuss high-priority findings at a multidisciplinary cardiac genetics meeting and confirm in an appropriately accredited laboratory before clinical return.

**Statistical caution:** An enrichment result (e.g. Fisher's exact test for intermediate repeat alleles) does not establish that those alleles are pathogenic or causal in an individual. Evaluate ancestry, technical ascertainment, phenotype, repeat length distribution, multiple testing and effect-size uncertainty.

## 12. Reporting and reproducibility

Minimum recommended columns:

| Field | Example / description |
|---|---|
| Sample ID | Pseudonymised study identifier |
| Clinical category | Autopsy-positive / suspicious / indeterminate, with precise definitions |
| Phenotype | SCD, ECG findings, cardiac pathology, extracardiac features |
| Reference | Exact GRCh38 FASTA accession/build and checksum |
| Data provenance | Sequencer, library prep, read length, coverage, alignment pipeline |
| Tool | EH/EHdn executable version and parameters |
| Catalog | JSON source, date, checksum, reference assembly |
| Locus | Gene, genomic coordinates, repeat motif, transcript |
| Alleles | Estimated repeat counts and uncertainty; indicate unknowns |
| Evidence | Spanning/flanking/IRR counts, QC and EHdn outlier statistic |
| Inheritance | AD/AR/X-linked; sex chromosome ploidy; phase if relevant |
| Interpretation | Known disease range, phenotype concordance, gene validity |
| Confirmation | Method, laboratory, confirmation status |
| Classification | Research candidate / confirmed expansion / clinical interpretation |

Record tool versions, checksums, PBS job IDs, resource use and log files. Preserve only de-identified metadata in public repositories. Obtain ethics/data-governance approval for storage, processing and any data sharing.

## 13. Limitations

- Short-read WGS cannot reliably resolve all large, GC-rich, interrupted or complex expansions; sensitivity depends on read length, insert size, coverage and locus architecture.
- EH only interrogates catalogued loci; catalog completeness and reference compatibility are critical.
- EHdn does not genotype repeat alleles and can produce imprecise genomic locations; robust comparator data are essential.
- PCR amplification, GC bias, alignment, mosaicism and tissue specificity can influence results.
- A negative short-read result **does not exclude** a pathogenic repeat expansion.
- Pathogenic thresholds and genotype–phenotype correlations are **locus-specific**; not every heterozygous expansion is disease-causing.
- Correlation with cardiac phenotypes is not sufficient to attribute an individual SCD to a repeat expansion.
- Software availability on Gadi and PBS project/storage directives must be verified locally before production use.

## 14. References and resources

### Software and HPC

- [Illumina ExpansionHunter — repository and documentation](https://github.com/Illumina/ExpansionHunter)
- [ExpansionHunter usage documentation](https://github.com/Illumina/ExpansionHunter/blob/master/docs/03_Usage.md)
- [Illumina ExpansionHunter Denovo — repository and documentation](https://github.com/Illumina/ExpansionHunterDenovo)
- [EHdn introduction and workflow](https://github.com/Illumina/ExpansionHunterDenovo/blob/master/documentation/00_Introduction.md)
- [EHdn outlier analysis](https://github.com/Illumina/ExpansionHunterDenovo/blob/master/documentation/04_Outlier_quickstart.md)
- [STRipy repeat expansion database](https://stripy.org/database)
- [NCI Gadi HPC systems](https://www.nci.org.au/infrastructure/hpc-systems)
- [NCI Gadi software modules introduction (Sydney Informatics Hub)](https://sydney-informatics-hub.github.io/training.gadi.intro/06-Software/index.html)

### Cardiac and disease mechanisms

- [DMPK 3′ untranslated repeat expansions in unexplained sudden cardiac death in the young (*European Heart Journal*, 2026)](https://academic.oup.com/eurheartj/article/47/33/4697/8586462)
- [Electrocardiographic abnormalities and sudden death in myotonic dystrophy type 1 (*NEJM*, 2008)](https://www.nejm.org/doi/full/10.1056/NEJMoa062800)
- [CTG repeat length underlying cardiac events and sudden death in DM1 (PubMed)](https://pubmed.ncbi.nlm.nih.gov/39391712/)
- [Relation of CTG expansion to conduction abnormalities and sudden death (PubMed)](https://pubmed.ncbi.nlm.nih.gov/14678805/)

### Suggested repository citation

When publishing, cite the **original methods papers** for EH and EHdn in addition to the exact software release, the variant catalog source, the relevant clinical references and the repository commit. This README is a research workflow and is not a validated diagnostic protocol.

---

**Status:** Research-ready template; requires local path/module validation, catalog selection, pilot benchmarking and institutional approval before production use.
