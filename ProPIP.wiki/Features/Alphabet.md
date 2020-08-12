[Back](./Index.md) | [Home](../Home.md)

---
#  Alphabet options
---


    alphabet={DNA|RNA|Protein)|Codon(letter={DNA|RNA},type={Standard|EchinodermMitochondrial|InvertebrateMitochondrial|VertebrateMitochondrial})}
                                                            The alphabet to use when reading sequences. DNA and RNA alphabet can in addition take
                                                            an argument: **bangAsgap={bool}**
                                                            Tell is exclamation mark should be considered as a gap character. The default
                                                            is to consider it as an unknown character such as 'N' or '?'.
    genetic_code={translation table}                        Where the translation table specifies the code to use, either as a text description,
                                                            or as the NCBI number. The following table give the currently implemented codes
                                                            with their corresponding names:
                                                            Standard                    1
                                                            VertebrateMitochondrial     2
                                                            YeastMitochondrial          3
                                                            MoldMitochondrial           4
                                                            InvertebrateMitochondrial   5
                                                            EchinodermMitochondrial     9
                                                            AscidianMitochondrial       13
                                                            The states of the alphabets are in alphabetical order.
