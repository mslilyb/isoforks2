Isoforms
========

## Notes ##

Redesign the isoform generators

- suffix tree algorithm
- scoring during generation
- option to keep N-best isoforms
- probabilities at the end
- self gff output

## Programs ##

### Python

Python-based applications use `isoform.py` for common functions.

+ `cmpiso` - compares collections of isoforms in GFF
+ `geniso` - generates isoforms and their probabilities
+ `modelbuilder` - creates the various model files
+ `optiso` - optimizes model parameters with a genetic algorithm
+ `run_apc` - builds Makefile to run `isoformer` and `optiso` on apc set

### C

The genomikon repo contains a couple of faster implementations in the
`isoformer` directory.

+ `isoformer` - this is the same as `geniso` but ~100x faster
+ `isocounter` - as above, but only counting, not calculating probabilities
+ `isorandom` - counting isoforms in random sequences

## Utilities ##

+ `conformity.py` - compares outputs of `geniso` and `isoformer`
+ `optiso-mp` - multi-processing version with some odd bugs
+ `speedo.py` - compares speeds of `geniso` and `isoformer`
+ `summary.py` - creates TSV of the apc set

## Data ##

Data collection is described in `datacore2024/project_splicing`. The 665 genes
of the APC set are in the `apc` directory.

## Models ##

See the `models` directory for standard models and `modelbuilder` for how to
build the models.
