[Back](./Home) | [Home](./Home)

---
# Releases
---


#### Version 1.0.8 / 20.07.2018

> #### Task
- [[MIN-45](https://devops.lorenzogatti.me/browse/MIN-45)] - Update codebase to support BPP libraries version 2.4.0

> #### Bug
- [[MIN-46](https://devops.lorenzogatti.me/browse/MIN-46)] - Codon models MX exit code 11
- [[MIN-47](https://devops.lorenzogatti.me/browse/MIN-47)] - Bug affecting double number comparison in function ComparisonUtils::areLogicallyEqual(double a, double b)
- [[MIN-48](https://devops.lorenzogatti.me/browse/MIN-48)] - PIP initial parameter estimates from data cli argument should be compatible with bpp standards


#### Version 1.0.7 / 02.07.2018

> #### Task
- [[MIN-43](https://devops.lorenzogatti.me/browse/MIN-43)] -         Add ML/wML to Initial tree topology estimator schemes (pairwise distances)
- [[MIN-44](https://devops.lorenzogatti.me/browse/MIN-44)] -         Pairwise distance estimation should be computed under nested models (+PIP)
    
> #### Bug
- [[MIN-40](https://devops.lorenzogatti.me/browse/MIN-40)] -         Initial likelihood value does not match with MSA score
- [[MIN-41](https://devops.lorenzogatti.me/browse/MIN-41)] -         Initial tree topology cannot be performed under PIP model. Distance measure is missing
- [[MIN-42](https://devops.lorenzogatti.me/browse/MIN-42)] -         Optimisation of branches during initial tree estimation cannot be performed under PIP
            
> #### Story
- [[MIN-8](https://devops.lorenzogatti.me/browse/MIN-8)] -         Design new parameter optimization routine to support more complex scenarios (using lk tolerance)


#### Version 1.0.5 / 08.06.2018
> #### Story
- Testing and tuning model performances to infer phylogenetic trees under nucleotide models (+PIP)
  
#### Version 1.0.4 / 08.05.2018
> #### Bugs
- [[MIN-17]](https://devops.lorenzogatti.me/browse/MIN-17) -         Bug affecting exchangeability and generator matrix filling. 

> #### Sub-task
- [[MIN-29](https://devops.lorenzogatti.me/browse/MIN-29) -         Update newick exported to support multiple node attribute and internal node name
- [[MIN-31]](https://devops.lorenzogatti.me/browse/MIN-31) -         Implementing new measure nh/ng according to ng/nh indices

> #### Story
- [[MIN-1]](https://devops.lorenzogatti.me/browse/MIN-1) -         Extending PIP model to support codon models
- [[MIN-2]](https://devops.lorenzogatti.me/browse/MIN-2) -         Adding parametric frequency rates (+F, +FQ, +F1X4 +F3X4) sets to codon models under PIP
- [[MIN-3]](https://devops.lorenzogatti.me/browse/MIN-3) -         Activating phylogenetic+alignment inference under codon models
- [[MIN-4]](https://devops.lorenzogatti.me/browse/MIN-4) -         Testing integration of codon models in miniJATI pipeline
- [[MIN-23]](https://devops.lorenzogatti.me/browse/MIN-23) -         Add JSON exporter to miniJATI
- [[MIN-30]](https://devops.lorenzogatti.me/browse/MIN-30) -         Implement new measure to detect clades with high-level of divergency based on InDel patterns