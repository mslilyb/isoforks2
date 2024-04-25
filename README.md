Isoforms
========

## Data ##

Data collection is described in `datacore2024/project_splicing`.

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

### genomikon

The genomikon library contains a couple of faster implementations in the
`isoformer` directory.

+ `isoformer` - the standard algorithm for all reasonable isoforms
+ `isocounter` - as above, but only counting, not calculating
+ `isorandom` - counting isoforms in random sequences


## To Do ##

+ Length models need to be done better
+ Isoform comparison function shouldn't need `deepcopy`
+ `optiso-mp` is bugged for no unknown reasons
