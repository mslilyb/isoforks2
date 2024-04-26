Isoforms
========

## Data ##

Data collection is described in `datacore2024/project_splicing`. The 665 genes
of the APC set are in the `apc` directory.

## Models ##

See the `models` directory for standard models and `modelbuilder` for how to
build the models.

## Sources ##

### isoform.py

Python-based applications use `isoform.py` for common functions.

+ `cmpiso` - compares collections of isoforms in GFF
+ `geniso` - generates isoforms and their probabilities
+ `modelbuilder` - creates the various model files
+ `optiso` - optimizes model parameters with a genetic algorithm
+ `optiso-mp` - multi-processing version (bugged currently)
+ `run_optiso` - creates a build directory and makefile for running many jobs

### genomikon

The genomikon library contains a couple of faster implementations in the
`isoformer` directory.

+ `isoformer` - this is the same as `geniso`
+ `isocounter` - as above, but only counting, not calculating probabilities
+ `isorandom` - counting isoforms in random sequences

## To Do ##

+ `optiso-mp` is bugged for no unknown reasons
+ `optiso`
	+ produces negative weights... do we want this?
	+ is volatile, many good solutions
	+ icost might need to be bounded
	+ maybe any very good solutions should not attempt fine tuning
