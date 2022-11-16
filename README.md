# tkp-finder

[![PyPI - Version](https://img.shields.io/pypi/v/tkp-finder.svg)](https://pypi.org/project/tkp-finder)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/tkp-finder.svg)](https://pypi.org/project/tkp-finder)

`tkp-finder` is a CLI tool to discover and annotate tandem protein kinases.

It's based on [lXtractor]() -- a general-purpose library for data mining from sequences and structures.
The latter is under active development, so bugs are possible.

-----

**Table of Contents**

- [Installation](#installation)
- [License](#license)
- [Usage](#usage)

## Installation

```console
pip install tkp-finder
```

## License

`tkp-finder` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.

## Usage

The installation should make the script `tkp-finder` globally available.
The interface has two commands:

The `setup` command will download and prepare HMM models for annotation.

```bash
→ tkp-finder setup --help

Usage: tkp-finder setup [OPTIONS]

  Command to initialize the HMM data needed for TKPs' annotation.

Options:
  -H, --hmm_dir DIRECTORY  Path to a directory to store hmm-related data.
                           [required]
  -d, --download           If True, download the Pfam data from interpro.
  -q, --quiet              Disable verbose output.
  --path_pfam_a FILE       A path to downloaded Pfam-A HMM profiles. By
                           default, if `download` is ``False``,will try to
                           find it within the `hmm_dir`.
  --path_pfam_dat FILE     A path to downloaded Pfam-A (meta)data file. By
                           default, if `download` is ``False``,will try to
                           find it within the `hmm_dir`.
  -h, --help               Show this message and exit.
```

For the first-time usage, invoke

```bash
→ tkp-finder setup -H hmm -d
```

This will download Pfam-A HMMs and accompanying metadata, and split the models into categories.
The resulting directory:

```bash
→ tree -L 2 hmm

hmm
├── PF00069.hmm
├── Pfam-A.hmm
├── Pfam-A.hmm.dat
├── pfam_entries.tsv
└── profiles
    ├── Coiled-coil
    ├── Disordered
    ├── Domain
    ├── Family
    ├── Motif
    ├── Repeat
    └── unknown
```

To dicover and annotate TKPs, refer to `tkp-finder find` command:

```bash
→ tkp-finder find --help

Usage: tkp-finder find [OPTIONS] [FASTA]...

Options:
  -H, --hmm_dir DIRECTORY    Directory with HMM profiles. Expected to contain
                             `profiles` dir and target PK profile
                             (PF00069.hmm). See `tkp-finder setup` on how to
                             prepare this dir.
  -t, --hmm_type TEXT        Which HMM types to use for annotating the
                             discovered TKPs. The names must correspond to
                             folders within he `hmm_dir`.  [default: Family,
                             Domain, Motif]
  -p, --pk_profile FILE      A path to the PK HMM profile. By default, will
                             try to find it within the `hmm_dir`.
  -m, --motif TEXT           A motif to discriminate between PKs and pseudo
                             PKs. This corresponds to the following conserved
                             elements::  (1) b3-Lys(2) aC-helix Glu(3-4-5) HRD
                             motif(6-7-8) DFG motif  [default: KEXXDDXX]
  -o, --output DIRECTORY     Output directory to store the results. Be
                             default, will store within `./tkp-finder`.
  -n, --num_proc INTEGER     The number of cpus for data parallelism: each
                             input fasta will be annotated within separate
                             process. HINT: one may split large fasta files
                             for faster processing.
  -q, --quiet                Disable logging and progress bar
  --pk_map_name TEXT         Use this name for the protein kinase domain.
                             [default: PK]
  --ppk_map_name TEXT        Use this name for pseudo protein kinases.
                             [default: PPK]
  --min_domain_size INTEGER  The minimum number of amino acid residues within
                             a PK domain.  [default: 150]
  --min_domains INTEGER      The number of domains to classify a protein as
                             TKP.
  -h, --help                 Show this message and exit.
```
