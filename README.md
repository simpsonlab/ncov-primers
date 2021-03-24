# ncov-primers

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


`ncov-primers` provides a set of tools for processing a primer BED file and
generating an amplicon BED file for further processing.

## Installation
```
pip install .
```


## Usage
After installing the package, you can run `bin/create_amplicons.py`:

```
$ create_amplicons.py --file /path/to/primer/bed
                      --output <name of output bed file>
                      --type <type of amplicons to generate>
```

The `--type` argument can be one of:

* unique_amplicons (unique regions removing primers)

* full (includes primers)

* no_primers (primers are removed)

