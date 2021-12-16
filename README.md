Code for [Walker, Saunders, Rai et al., (2021).](#)

If you use the code or data from this repository, please cite our manuscript as described
in `CITATION.cff` or using GitHub's Cite this repository link shown on the right.

## Instructions

Due to complexity of many unique analyses and sensitivity of raw data, this repository
cannot be run in a "one-click-all" fashion. If you wish to reproduce results from our
manuscript, please first request access to the raw data from the European Genome-Phenome
Archive (EGA) submission.

Once raw data is available, there are three broad topics of analyses:
- Bulk RNA-seq
- Multiome RNA and ATAC
- Integration of the data modalities

These analyses are done in the following directories:

```
.
├── CITATION.cff
├── multiome_mapping-qc
│   ├── main --> Processing of "main" batch of Multiome
│   └── pilot --> Processing of the "pilot" batch of Multiome
├── rnaseq --> Processing of the bulk RNA-seq data
│   ├── bin
│   ├── control
│   ├── Rakefile
│   ├── scripts
│   └── src
```

### Directory organization

Most of the analysis directories are managed using,
[`makebio`](https://github.com/raivivek/makebio), a tool for managing computational
biology projects.

For example, the `rnaseq` analysis directory has the following structure:

```
rnaseq
├── bin
├── control
├── Rakefile
├── scripts
└── src
```

1. `src/` - source code corresponding to external packages, if used
1. `control/` - contains distinct "sets" of analyses
1. `scripts/` - Rmd notebooks, scripts
1. `bin/` - Executable scripts used by different components of the pipeline
1. `Rakefile` - Driver script


In order to run an analyses, please inspect the `Rakefile` driver script first. That
should give you an idea of how to execute each analyses. Within `control/` directory,
different analyses are contained in distinct directories with each having similar
structure. Please ensure that appropriate data paths and variables are set before
running any code.

## Data Access

Raw data submission to European Genome-Phenome Archive (EGA) is _in progress_
(RRID:SCR\_004944).

## Contact

Please reach out to [Steve Parker](mailto:scjp@umich.edu) for any data and code related
questions.
