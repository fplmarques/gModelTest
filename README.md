# gModelTest
## Description:
This program uses Garli to select substitution models for phylogenetic analyeses usin Maximum Likelihood \
as optimality criterium. As such it requires [Garli](https://code.google.com/archive/p/garli/).

## Usage:
``` $ gmodel_test.pl [option] aligned_sequences.fas ```
### Options:
-fixed|-fix|-f  -> uses fixed topology to compute ML scores for all 88 substitution model.
-parallel|-par|-p -> implement parallel execution using 8 cores.
