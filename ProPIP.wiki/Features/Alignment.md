[Back](./Index.md) | [Home](https://github.com/acg-team/ProPIP/wiki/ProPIP:-Progressive-Multiple-Sequence-Alignment-with-Poisson-Indel-Process.md)

---
#  Alignment options
---

    alignment = {boolean}                                   If set to true then it inferred the multiple sequence alignment.
    
    alignment.version = {ram|sb}                            RAM: version of Dynamic Programming that uses less RAM memory. 
                                                            SB (Stochastic Backtracking version): this version generates a 
                                                            distribution of suboptimal MSAs at each node according to their 
                                                            probability.  

    alignment.sb_temperature = {real}                       SB is parameterised by a temperature 'T' which tunes the deviation 
                                                            from the optimal alignment. For 'T = 0.0' SB returns the optimal 
                                                            alignment, falling back to classical DP. By setting 'T > inf', each 
                                                            alignment becomes equiprobable and the solution is therefore purely 
                                                            random. In the range '0 < T < inf' the parameter controls the devia-
                                                            tion from the optimal alignment allowing, gradually, the generation 
                                                            of sub-optimal alignments.  

    alignment.sb_solutions = {integer}                      Option used only with the 'alignment.version = sb' version. 
                                                            Tells the number of suboptimals the algorithm generates at each 
                                                            node.

    seed = {real}                                           Sets the seed value of the random number generator.



