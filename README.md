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
                      --idloc <index of id location in primer name>
                      --keep
```

The `--type` argument can be one of:

* unique_amplicons (unique regions removing primers)

* full (includes primers)

* no_primers (primers are removed)

The `--keep` flag is used to keep a temporary directory used to construct
the amplicon file.  This can be used for debugging purposes.


## Primer names
The primer names within the BED file require a specific format.  The name is
typically found in column 4 of a standard BED file.  The name must contain:

```
* a unique integer ID

* LEFT/RIGHT to indicate the direction of the primer

* components separated by _ (underscore)
```

Examples of valid primer names are: `ncov_1_RIGHT` and `ncov_22_LEFT`.  In this
example, the `--idloc` argument should be set to `1` to indicate the location
within name for the id (Note: 0-based index).

