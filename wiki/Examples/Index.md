## 1: Inferring influenza tree using an indel-aware substitution model

- **Dataset:** example_dataset_nuc.fa (influenza real nucleotide dataset)
- **Initial tree:** distance based (bioNJ)
- **Substitution model:** PIP(with initial rates estimated from the data)+GTR(with initial rates estimated from the data)
- **Optimisation (numerical parameters):** D-BFGS(derivatives=BFGS)
- **Optimisation (topology):** Mixed strategy (NNI + F/VF-NNI + SPR) from random starting nodes (n=50) and branch lenght optimisation during tree search perfomed using BFGS. The process will span automatically on 10 threads.

### Prepare the configuration file

We list all the parameters required by Castor in a text file named `indel_aware.txt` (the order of the parameters is not relevant). For more information on the syntax of the parameter, please check the [documentation](https://bitbucket.org/lorenzogatti89/castor/wiki/Features/Index).

```
analysis_name = gene_HA
model_description = GTR+PIP
input_folder = /fullpath-input-directory/
output_folder = /fullpath-output-directory/

alphabet=DNA
alignment=false
input.sequence.file=$(input_folder)/example_dataset_nuc.fa
input.sequence.sites_to_use=all
init.tree=distance
init.distance.method=bionj
model=PIP(model=GTR(initFreqs=observed),initFreqs=observed)
rate_distribution=Constant
optimization=D-BFGS(derivatives=BFGS)
optimization.max_number_f_eval=5000
optimization.tolerance=0.001
optimization.final=bfgs
optimization.topology=true
optimization.topology.algorithm=Mixed(coverage=best-search,starting_nodes=Hillclimbing(n=50),max_cycles=100,tolerance=0.001,brlen_optimisation=BFGS,threads=10)
output.msa.file=$(output_folder)/$(analysis_name).$(model_description).out.align.fasta
output.tree.file=$(output_folder)/$(analysis_name).$(model_description).out.tree.tree
output.estimates.file=$(output_folder)/$(analysis_name).$(model_description).out.estimates.log
output.estimates.format=json
support=none

```

### Execute the analysis
```
$ Castor params=/directory/of/the/indel_aware.txt

```


### Results

The analysis should end with the following results. 

    Output alignment to file.............:   /fullpath-output-directory/gene_HA.GTR+PIP.out.align.fasta
    Output tree file ......................: /fullpath-output-directory/gene_HA.GTR+PIP.out.tree.tree
    Output tree format ....................: Newick
    Final Log likelihood...................: -43363.9641982363
    PIP.GTR.a..............................: 0.988982
    PIP.GTR.b..............................: 0.213126
    PIP.GTR.c..............................: 0.148805
    PIP.GTR.d..............................: 0.413321
    PIP.GTR.e..............................: 0.168315
    PIP.GTR.theta..........................: 0.417233
    PIP.GTR.theta1.........................: 0.487393
    PIP.GTR.theta2.........................: 0.533194
    PIP.lambda.............................: 63.5404
    PIP.mu.................................: 0.0421603
    PIP.intensity..........................: 2.67888
    WARNING!!! This parameter has a value close to the boundary: BrLen20(1e-06).
    WARNING!!! This parameter has a value close to the boundary: BrLen22(1e-06).
    WARNING!!! This parameter has a value close to the boundary: BrLen28(1e-06).
    WARNING!!! This parameter has a value close to the boundary: BrLen30(1e-06).
    WARNING!!! This parameter has a value close to the boundary: BrLen82(1e-06).
    Output estimates format................: json
    Output estimates to file...............: /fullpath-output-directory/gene_HA.GTR+PIP.out.estimates.log
    Total execution time: 0.000000d, 10.000000h, 13.000000m, 35.000000s.




## 2: Inferring influenza tree using a substitution model

- **Dataset:** example_dataset_nuc.fa (influenza real nucleotide dataset)
- **Initial tree:** distance based (bioNJ)
- **Substitution model:** GTR(with initial rates estimated from the data)
- **Optimisation (numerical parameters):** D-BFGS(derivatives=BFGS)
- **Optimisation (topology):** Mixed strategy (NNI + F/VF-NNI + SPR) from random starting nodes (n=50) and branch lenght optimisation during tree search perfomed using BFGS. The process will span automatically on 10 threads.

### Prepare the configuration file

We list all the parameters required by Castor in a text file named `indel_non_aware.txt` (the order of the parameters is not relevant). For more information on the syntax of the parameter, please check the [documentation](https://bitbucket.org/lorenzogatti89/castor/wiki/Features/Index).

```
analysis_name = gene_HA
model_description = GTR
input_folder = /fullpath-input-directory/
output_folder = /fullpath-output-directory/

alphabet=DNA
alignment=false
input.sequence.file=$(input_folder)/example_dataset_nuc.fa
input.sequence.sites_to_use=all
init.tree=distance
init.distance.method=bionj
model=GTR(initFreqs=observed)
rate_distribution=Constant
optimization=D-BFGS(derivatives=BFGS)
optimization.max_number_f_eval=5000
optimization.tolerance=0.001
optimization.final=bfgs
optimization.topology=true
optimization.topology.algorithm=Mixed(coverage=best-search,starting_nodes=Hillclimbing(n=50),max_cycles=100,tolerance=0.001,brlen_optimisation=BFGS,threads=10)
output.msa.file=$(output_folder)/$(analysis_name).$(model_description).out.align.fasta
output.tree.file=$(output_folder)/$(analysis_name).$(model_description).out.tree.tree
output.estimates.file=$(output_folder)/$(analysis_name).$(model_description).out.estimates.log
output.estimates.format=json
support=none

```

### Execute the analysis
```
$ Castor params=/directory/of/the/indel_non_aware.txt
```


### Results

The analysis should end with the following results. 

    Output alignment to file.............:   /fullpath-output-directory/gene_HA.GTR.out.align.fasta
    Output tree file ......................: /fullpath-output-directory/gene_HA.GTR.out.tree.tree
    Output tree format ....................: Newick
    Final Log likelihood...................: -35004.4671436324
    GTR.a..................................: 1.21461
    GTR.b..................................: 0.2429
    GTR.c..................................: 0.178708
    GTR.d..................................: 0.422956
    GTR.e..................................: 0.184716
    GTR.theta..............................: 0.422552
    GTR.theta1.............................: 0.56559
    GTR.theta2.............................: 0.55595
    WARNING!!! This parameter has a value close to the boundary: BrLen108(1e-06).
    WARNING!!! This parameter has a value close to the boundary: BrLen112(1e-06).
    Output estimates format................: json
    Output estimates to file...............: /fullpath-output-directory/gene_HA.GTR.out.estimates.log
    Total execution time: 0.000000d, 8.000000h, 36.000000m, 16.000000s.


### Dataset

You can download the example dataset [here](https://bitbucket.org/lorenzogatti89/castor/downloads/example_dataset_nuc.fasta)