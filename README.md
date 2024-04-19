Isoforms
========

Consolidating the various splicing things in one place.

+ Data
	+ Genes
		+ Exons
		+ Introns
		+ Splice site usage from RNA-seq
		+ Exon coverage from RNA-seq
		+ Whole isoforms?
	+ Organisms
		+ A.thaliana
		+ C.elegans
		+ D.melanogaster
		+ Others?
+ Models
	+ PWM
	+ Markov/k-mer
	+ Length
	+ Others? WAM, DT
+ Algorithms
	+ APC
	+ BLI
	+ FB
	+ SV
+ Implementations
	+ Model files
	+ C genomikon
	+ Go Alan's
	+ Python (Ian, Ismael)
+ Results

## Data ##

## Models ##

The model file formats are specified as follows.

### PWM

### MM

### Length

## Algorithms ##

## Implementations ##

### genomikon

The genomikon library contains a couple of implementations in the `isoformer`
directory.

+ `isoformer` - the standard algorithm for all reasonable isoforms
+ `isocounter` - as above, but only counting, not calculating
+ `isorandom` - counting isoforms in random sequences

## Results ##
