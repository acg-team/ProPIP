[Back](../Index.md) | [Home](https://github.com/acg-team/ProPIP/blob/master/ProPIP.wiki/ProPIP-Progressive-Multiple-Sequence-Alignment-with-Poisson-Indel-Process.md)

---
#  Topology optimisation options
---


To perform a tree-search optimisation step, one should add the following arguments:

    optimisation.topology=<bool>	                            Enable the tree-search
    optimisation.topology.algorithm={Swap(),Phyml(),Mixed()}	Method to use during the tree-search



miniJATI inherits the following algorithms to perform tree-search from the Tree-Search-Heuristics Library (TSHLib).

    Swap(...)
    Phyml(...)
    Mixed(...)

The available options for the methods above mentioned are:


    coverage={nni-search|spr-search|tbr-search|best-search}                     Tree-space coverage options
    starting_nodes={Greedy(),Hillclimbing(n=4),SwarmParticle(cores=100)}        Enumerating method to define candidate moves
    max_cycles=<int> (i.e. 50)                                                  Max number of tree-search cycles
    tolerance=<double> (i.e 0.01)                                               Likelihood tolerance to limit the number of Tree-Searchcycles
    brlen_optimisation={Newton|Gradient|BFGS|Brent}                             Branch lenght optimisation method during Tree-Search


The available options for the `starting_nodes` argument are the following:

    Greedy()                                All the nodes on the topology are used to enumerate candidate tree-rearrangement moves
    Hillclimbing(n=4)                       Moves are enumerated starting from a limited set of starting nodes (i.e. n=4)
    SwarmParticle(cores=100)


Tree rearrangment operations to apply during tree-search are encoded in the available options for the `coverage` argument, as follows:

    nni-search
    spr-search
    tbr-search
    best-search

Method to use to optimise the Branch Lengths during the topology optimization step can be set using the available options for the `brlen_optimisation` argument, as follows:

    Newton
    Gradient
    BFGS
    Brent

Automatically the algorithm instantiates analytical or numerical derivatives to guide the optimization process.


### Examples

optimization.topology=true optimization.topology.algorithm=Swap(coverage="best-search",starting_nodes=Greedy(),max_cycles=50,tolerance=0.1,brlen_optimisation=Brent)
