# Third-party software: citation and licence rules

This note explains how third-party software is handled in this repository.
Three distinct obligations apply, and they are often confused:

| Obligation | What it means | Where it is satisfied |
|---|---|---|
| **Citation** (academic credit) | the tool's paper is cited where the work is published | reference list |
| **Attribution** (licence requirement) | the original copyright notice is preserved | alongside the redistributed file, plus this directory |
| **Licence compliance** | the terms of the tool's licence are followed | depends on whether code is *called* or *copied* |

---

## 1. Tools that are only called — no action beyond citation

Most software in this pipeline is executed as an external program. Its source is **not** copied
into this repository, so nothing needs to be relicensed and no licence text needs to be reproduced.

For each such tool, three things are recorded:

1. **the exact version actually run** — in `environment/versions.tsv`
2. **the citation** — see the reference list of the publication that uses this code
3. **the parameters** — in the README of the step that calls the tool

Example, as it appears in `environment/versions.tsv`:

```text
tool	version	purpose	source
HISAT2	2.2.1	RNA-seq read alignment	http://daehwankimlab.github.io/hisat2/
```

A note on copyleft: running a GPL-licensed program from the command line does **not**
place this repository under the GPL. Copyleft obligations attach to distributing the GPL code
itself or a derivative of it — not to invoking the program on your own data. The situation
changes if GPL source is *copied into* this repository (see section 2).

---

## 2. Code written by others that IS redistributed — attribution required

### 2.1 `polish.sh` — T2T-Polish automated polishing

| Field | Value |
|---|---|
| Upstream repository | https://github.com/arangrhie/T2T-Polish |
| Upstream path | `automated_polishing/` |
| File authorship | Ivan Sovic and Ann Mc Cartney, September 2021 |
| Licence | **Public domain** — U.S. Government Work, NHGRI/NIH |
| Licence text | reproduced in `third_party/T2T-Polish-LICENSE.txt` |
| Used for | final polishing of the Gifu and MG20 assemblies |

**Licence status: resolved.** The upstream `LICENSE` is a PUBLIC DOMAIN NOTICE stating that the
software is United States Government Work, cannot be copyrighted, and that restrictions cannot be
placed on its present or future use. Redistribution in this repository is therefore permitted.
The notice asks that the authors be cited.

The file is therefore included here under the following terms:

1. the original author header in the script is unchanged;
2. the PUBLIC DOMAIN NOTICE is reproduced alongside it
   (`third_party/T2T-Polish-LICENSE.txt`);
3. the repository-root MIT licence does **not** apply to it — it is public domain, not ours to license;
4. the authors are cited (see 2.2).

### 2.2 Citation

The upstream README states explicitly:

> The original script used in **McCartney et al, 2021**
> (https://doi.org/10.1101/2021.07.02.450803) has been updated to use the Merfin release
> version and contains minor corrections to Winnowmap alignments.

The citation requested by the authors of this workflow is therefore **Mc Cartney et al.**, not the
human T2T genome paper. The published version of that preprint is:

> Mc Cartney, A. M. et al. Chasing perfection: validation and polishing strategies for
> telomere-to-telomere genome assemblies. *Nature Methods* **19**, 687–695 (2022).
> https://doi.org/10.1038/s41592-022-01440-3

**Which citation to use:** the workflow is maintained in Arang Rhie's repository, but its own
documentation asks users to cite **Mc Cartney et al.** rather than the human T2T genome paper. The
individual tools invoked by the script are Winnowmap2, falconc, Racon, Merfin and bcftools.

### 2.3 Version and parameter differences

The copy of `polish.sh` used in this study **is not identical to either version currently
published upstream**. The three differ as follows:

| Step | This study (`polish.sh`) | Upstream `automated-polishing.sh` | Upstream `automated-polishing-legacy.sh` |
|---|---|---|---|
| Repetitive k-mer collection | `meryl count k=21` | `meryl count k=15` | none (uses a pre-built `bad_mers.txt`) |
| Merfin invocation | `merfin -polish … -peak 106.7` | `merfin -polish … -peak 106.7` | `merfin … -peak 106 -vmer -disable-kstar` (no `-polish`) |
| Consensus input | `*.polish.vcf` | `*.polish.vcf` | `*.filter.vcf` |
| Assembly-to-assembly variant block | absent | absent | present (minimap2 + paftools) |

The study's script matches the **current** upstream structure but retains the **k=21**
repeat-mer setting, which upstream later changed to k=15 as part of the stated
"minor corrections to Winnowmap alignments".

The value actually used in this study is **k=21**, which is what
`scripts/01_assembly/07_polish/README.md` documents. The file itself is
redistributed unmodified, with its original author header intact.

### 2.4 Other scripts in this directory

| Script | Basis for inclusion |
|---|---|
| `01_assembly/06_gapclose/quarTeT.sh` | wraps the third-party quarTeT package |
| `02_annotation/01_repeat_edta/process.sh` | wraps EDTA |
| `05_centromere/04_stainglass/run_stainglass.sh` | copied from the StainedGlass 0.6 workflow |
| `05_centromere/04_stainglass/process.sh` | copied from the StainedGlass 0.6 workflow |
| `04_rDNA/01_locus_and_arrays/result/process.sh` | calls StainedGlass 0.6 scripts directly |

A script that is a thin wrapper setting arguments for an external tool counts as our own
orchestration code and lives in `scripts/`, with the external tool cited. A script containing
copied source is recorded here with its licence.

---

## 3. Third-party code without an established licence

A third-party script whose licence cannot be established is not redistributed here; it is
described instead:

- name the authors and year,
- give the URL of the original source,
- state that it was used unmodified,
- list it in `environment/versions.tsv`.

A script that is essential to reproducing a result is re-implemented by its **logic** in our own
code with a comment pointing at the original, rather than copied as a file.

---

## 4. Third-party tools listed for completeness

The following tools are called by the pipeline and are **not** redistributed here.
Their versions and canonical sources are recorded in `environment/versions.tsv`:

hifiasm, gfatools, purge_dups, minimap2, BLAST+, TGS-GapCloser2, quarTeT, winnowmap, racon,
merfin, meryl, falconc (pbipa), bcftools, samtools, BUSCO, Merqury, HiC-Pro, cooler, cooltools,
EDTA, LTR_retriever, HISAT2, StringTie, isoseq3, TAMA, TransDecoder, BRAKER3, AUGUSTUS, GETA,
Helixer, AGAT, eggNOG-mapper, PfamScan, RSEM, JCVI, OrthoFinder, clusterProfiler, TeloComp,
Barrnap, bedtools, StainedGlass, Tandem Repeats Finder, CD-HIT-EST, MAFFT, MUMmer4, SyRI, pysam,
limma, WGCNA, Mfuzz, Cell Ranger, scVI, Leiden, scTenifoldKnk, SequenceServer, JBrowse 2,
Cirrocumulus.
