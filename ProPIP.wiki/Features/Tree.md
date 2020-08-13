[Back](./Index.md) | [Home](../Home.md)

---
#  Initial tree options
---

    init.tree={user|random|distance}                        Set the method for the initial tree to use.
                                                            The user option allows you to use an existing file passed via input.tree.file
                                                            This file may have been built using another method like neighbor joining or
                                                            parsimony for instance. The random option picks a random tree, which is handy
                                                            to test convergence.  This may however slows down significantly the optimization
                                                            process.
    init.distance.matrix.file={path}                        A distance matrix can be supplied instead of being computed from the alignment.
    init.distance.method={wpgma|upgma|nj|bionj|
                              infere_distance_matrix}       When distance method is required, the user can specify which algorithm to use.
                                                            With the option 'infere_distance_matrix' the algorithm calculates the distance 
                                                            matrix from the pairwise alignments and then with bionj algorithm infers the tree.
    init.brlen.method={method description}                  Set how to initialize the branch lengths. Available methods include:

    Input(midpoint_root_branch={boolean})                   Keep initial branch lengths as is. Additional argument specifies if the root
                                                            position should be moved to the midpoint position of the branch containing it.



    Equal(value={float>0})                                  Set all branch lengths to the same value, provided as argumemt.

    Clock                                                   Coerce to a clock tree.

    Grafen(height={{real>0}|input}, rho = {real>0})         Uses Grafen's method to compute branch lengths.
                                                            In Grafen's method, each node is given a weight equal to the number of underlying
                                                            leaves. The length of each branch is then computed as the difference of the weights
                                                            of the connected nodes, and further divided by the number of leaves in the tree.
                                                            The height of all nodes are then raised to the power of 'rho', a user specified value.
                                                            The tree is finally scaled to match a given total height, which can be the original
                                                            one (height=input), or fixed to a certain value (usually height=1). A value of
                                                            rho=0 provides a star tree, and the greater the value of rho, the more recent the
                                                            inner nodes.

    input.tree.check_root = {boolean}                       Tell if the input tree should be checked regarding to the presence of a root.




If the `init.tree metod=user`, then refer to the option you find [here](./Input.md)  
