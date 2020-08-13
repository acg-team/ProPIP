## 1: Inferring an MSA using an indel-aware substitution model from nucleotide input sequences

### Prepare the configuration file

We list all the parameters required by ProPIP in a text file named `indel_aware.txt` 
(the order of the parameters is not relevant). For more information on the syntax of the parameter.

```
analysis_name = TEST
model_description = JC69+PIP
alphabet=DNA
seed=1
alignment=true
alignment.version=ram
input.sequence.file=../tests/input/test_07/seqs.fasta
input.sequence.sites_to_use=all
init.tree=user
input.tree.file=../tests/input/test_07/tree.newick
model=PIP(model=JC69,lambda=15.8,mu=0.06)
rate_distribution=Constant
optimization=None
output.msa.file=../tests/output/test_07/msa.fasta
output.tree.file=../tests/output/test_07/tree.nwk
output.estimates.file=../tests/output/test_07/estimates.log
output.estimates.format=json

```

### Execute the analysis
```
$ ProPIP params=/directory/of/the/indel_aware.txt
```

## 2: Inferring an MSA using an indel-aware substitution model from amino-acids input sequences

### Prepare the configuration file

We list all the parameters required by ProPIP in a text file named `indel_aware.txt` 
(the order of the parameters is not relevant). For more information on the syntax of the parameter.

```
analysis_name = TEST
model_description = WAG01+PIP
alphabet=Protein
seed=1
alignment=true
alignment.version=ram
input.sequence.file=../tests/input/test_08/seqs.fasta
input.sequence.sites_to_use=all
init.tree=user
input.tree.file=../tests/input/test_07/tree.newick
model=PIP(model=WAG01,lambda=15.8,mu=0.06)
rate_distribution=Constant
optimization=None
output.msa.file=../tests/output/test_08/msa.fasta
output.tree.file=../tests/output/test_08/tree.nwk
output.estimates.file=../tests/output/test_08/estimates.log

```

### Execute the analysis
```
$ ProPIP params=/directory/of/the/indel_aware.txt
```

## 3: Inferring an MSA using an indel-aware substitution model from amino-acids input sequences. The tree and the indel rates are not provided and therefore inferred from the input sequences.  


The initial values of insertion and deletion rates of the PIP model are inferred from the unaligned 
input sequences, as an average of all rates estimated from pairwise alignments.  
The pairwise alignments are computed using the Needleman-Wunsch algorithm with gap opening and extension 
penalties for nucleotide sequences and a Grantham distance-based scoring method for amino acids. From the 
pairwise alignments, the indel rates are calculated by minimizing a Least Squares system (non-linear 
least-squares Marquardt-Levenberg algorithm.
Therefore the initial gap patterns are enforced through a system of equations with respect to the parameters 
of the PIP models, i.e., the insertion and deletion rates and the input sequences length.

Providing a realistic initial guide tree and indel rates helps to make the MSA inference more accurate. 
These can be provided by the user when known. If the guide tree is not provided then the tool computes a 
distance matrix from the pairwise alignments and infers a rooted guide tree using the BioNJ algorithm.

### Prepare the configuration file

We list all the parameters required by ProPIP in a text file named `indel_aware.txt` 
(the order of the parameters is not relevant). For more information on the syntax of the parameter.

```
alphabet=Protein
alignment=true
alignment.version=ram
input.sequence.file=../tests/input/test_10/seqs.fasta
input.sequence.sites_to_use=all
init.tree=distance
init.distance.method=infere_distance_matrix
model=PIP(model=WAG01)
rate_distribution=Constant
optimization=None
output.msa.file=../tests/output/test_10/msa.fasta
output.tree.file=../tests/output/test_10/tree.nwk
output.estimates.file=../tests/output/test_10/estimates.log
```

### Execute the analysis
```
$ ProPIP params=/directory/of/the/indel_aware.txt
```

