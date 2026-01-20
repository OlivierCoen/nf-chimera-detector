# nf-chimera-detector

[![GitHub Actions CI Status](https://github.com/nf-chimera-detector/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-chimera-detector/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-chimera-detector/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-chimera-detector/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.04.2-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-chimera-detector)

## Introduction

**nf-chimera-detector** is a bioinformatics pipeline that aims at detecting transposable elements in public genomic data from DNA viruses.

> [!TIP]
> Before running the pipelines, be sure to have Nextflow + Apptainer or Nextflow + Micromamba installed.
> See the prerequisites [here](docs/prerequisites.md).

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

> [!IMPORTANT]
> This pipeline may have to handle a very important amount of data, especially when fetching all data from big families. Make sure to have enough disk space available. Alternatively, you can add the `cleanup` profile to your profiles (e.g. `-profile apptainer,cleanup`), which will perform automatic cleanup of intermediate files using the `nf-boost` plugin. This solution should be used with caution, as it will prevent users from using the `-resume` parameter.

## Usage

Prepare a `.txt` file with the following format:

```
Pithoviridae
Ascoviridae
Nudiviridae
Mamonoviridae
```
and a `.fasta` file containing the sequences that you expect to find in your NGS data.

Now, you can run the pipeline using:

```bash
nextflow run OlivierCoen/nf-chimera-detector \
   -latest \
   -profile apptainer \
   --family_file <family file> \
   --target_db <target DB Fasta file> \
   --outdir <output directory>
```

>[!TIP]
>More detailed usage documentation is available [here](docs/usage.md).
>In particular, if you want to run the pipeline with your own data, see [Use local data only](docs/usage.md#2---use-local-data-only).

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

nf-chimera-detector was originally written by Olivier Coen.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-chimera-detector for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
