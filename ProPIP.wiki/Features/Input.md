[Back](./Index.md) | [Home](../ProPIP-Progressive-Multiple-Sequence-Alignment-with-Poisson-Indel-Process.md)

---
#  Input options
---


## Reading sequences

    input.sequence.file={path}                              The sequence file to use. (These sequences can also be not aligned).
    input.sequence.format={format}                          The sequence file format.
    input.sequence.sites_to_use={all|nogap|complete}        Tells which sites to use
    input.sequence.remove_stop_codons={boolean}             Removes the sites where there is a stop codon (default: 'yes')
    input.sequence.max_gap_allowed=100%                     It specifies the maximum amount of gap allowed per site.
    input.sequence.max_unresolved_allowed=100%              It specifies the maximum amount of unresolved states per site.
    input.site.selection={list of integers}                 Will only consider sites in the given list of positions, in extended format :
                                                            positions separated with ",", and "i:j" for all positions between i and j,
                                                            included.
    input.site.selection = {Sample(n={integer} [, replace={true}])}
                                                            Will consider {n} random sites, with optional replacement.

*The following formats are currently supported:*

    Fasta(extended={bool}, strictNames={bool})              The fasta format. The argument extended, default to 'no' allows to enable the HUPO-PSI
                                                            extension of the format. The argument strict_names, default to 'no', specifies that
                                                            only the first word in the fasta header is used as a sequence names, the rest of the
                                                            header being considered as comments.
    Mase(siteSelection={chars})                             The Mase format (as read by Seaview and Phylo_win for instance), with an optional site
                                                            selection name.
    Phylip(order={interleaved|sequential}, type={classic|extended}, split={spaces|tab})
                                                            The Phylip format, with several variations. The argument order distinguishes between
                                                            sequential and interleaved format, while the option type distinguished between the
                                                            plain old Phylip format and the more recent extention allowing for sequence names
                                                            longer than 10 characters, as understood by PAML and PhyML. Finally, the split
                                                            argument specifies the type of character that separates the sequence name from the
                                                            sequence content. The conventional option is to use one (classic) or more (extended)
                                                            spaces, but tabs can also be used instead.
    Clustal(extraSpaces={int})                              The Clustal format.
                                                            In its basic set up, sequence names do not have space characters, and one space splits
                                                            the sequence content from its name. The parser can however be configured to allow
                                                            for spaces in the sequence names, providing a minimum number of space characters is
                                                            used to split the content from the name. Setting extraSpaces to 5 for instance, the
                                                            sequences are expected to be at least 6 spaces away for their names.
    Dcse()                                                  The DCSE alignment format. The secondary structure annotation will be ignored.
    Nexus()                                                 The Nexus alignment format. (Only very basic support is provided)
    GenBank()                                               The GenBank not aligned sequences format.
                                                            Very basic support: only retrieves the sequence content for now, all features are
                                                            ignored.

### Reading trees

    input.tree.file={path}                                  The phylogenetic tree file to use.
    input.tree.format={Newick|Nexus|NHX}                    The format of the input tree file.
