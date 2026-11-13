# doh-jdj0303/vaper: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Alpha - [11/27/2023]
- Repo created

## Beta - [12/28/2023]
- Beta version released

## v1.0 - [06/05/2025]
- Version 1.0 released

## v2.0.0-rc.1 - [02/17/2026]
- Version 2.0.0 release candidate

## v2.0.0 - [02/17/2026]
- Version 2.0.0 release candidate

## Unreleased
- Alignments are now randomly subsampled with `rasusa aln` before consensus calling. `--cons_max_depth` sets the Rasusa target depth instead of the `samtools mpileup -d` limit, which kept the first reads in the coordinate-sorted BAM rather than a random subset.
- The full (`<sample>.bam`) and subsampled (`<sample>.subsampled.bam`) alignments are both saved to `<outdir>/<sample>/align/`. Read statistics are still calculated from the full alignment.
