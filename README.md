Isoforms
========

## Data ##

Data collection is described in `datacore2024/project_splicing`. The 665 genes
of the APC set are in the `apc` directory.

## Models ##

See the `models` directory for standard models and `modelbuilder` for how to
build the models.

## Sources ##

`conformity.py` compares pairs of implementations against each other to see if
they get the same results. There are minor math differences (1e-5) between
`geniso` and `isoformer`.

`speedo.py` compares the speeds of `geniso` and `isoformer`.

### isoform.py

Python-based applications use `isoform.py` for common functions.

+ `cmpiso` - compares collections of isoforms in GFF
+ `geniso` - generates isoforms and their probabilities
+ `modelbuilder` - creates the various model files
+ `optiso` - optimizes model parameters with a genetic algorithm
+ `optiso-mp` - multi-processing version (bugged currently)
+ `run_apc` - runs `isoformer` and `optiso` on the whole apc set

### genomikon

The genomikon library contains a couple of faster implementations in the
`isoformer` directory.

+ `isoformer` - this is the same as `geniso` but 40-90x faster
+ `isocounter` - as above, but only counting, not calculating probabilities
+ `isorandom` - counting isoforms in random sequences

## To Do ##

+ `optiso-mp` is bugged for unknown reasons
+ `optiso` may have many good solutions
