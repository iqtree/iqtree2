Identifying equally scoring trees in phylogenomics with incomplete data using Gentrius
=======

Gentrius is a deterministic algorithm to generate trees from unrooted incomplete subtrees. It phylogenetics it can be applied to trees inferred with common methods such as supermatrix and supertrees to generate a set of all trees that have the same induced subtrees as an inferred tree. This set is called a stand. Stands are extremely valuable, since under many scoring functions all trees from the same stand have identical analytical score. Read more in our manuscript.

Preprint is available from bioRxiv
* O. Chernomor, C. Elgert, A. von Haeseler (2023) Identifying equally scoring trees in phylogenomics with incomplete data using Gentrius

Analysis with Gentrius
-----------------------
### Basic Input
-----
There are several ways and inputs one can use to start the analysis. Here are the examples of command lines:

1.  __iqtree2 -gentrius <file 1> -pr_ab_matrix <file 2>__
2.  __iqtree2 -gentrius <file 1> -s <file 3> -p <file 4>__
3.  __iqtree2 -gentrius <file 5>__

    - File 1: a phylogenetic tree to be analysed

    - File 2: a species per locus presence-absence matrix. First line is a number of species and the number of loci. Other lines contain species name followed by 0's and 1's (with space separator) to indicate presence/absence of gene sequences for corresponding loci. 

            Example: 
            ---------------
            5 3  
            species_1 0 1 1  
            species_2 1 0 1  
            species_3 1 1 1  
            species_4 1 1 0  
            species_5 1 1 1  

    - File 3: an alignment file with concatenated loci
    - File 4: a partition info file for the alignment in file 3
    - File 5: contains a set of subtrees to be analysed

### Stopping Rules
-----
As the underlying problems are computationally hard, we use three different rules to stop the analysis. The thresholds can be changed by the corresponding options.

- __Rule 1:__ stop the analysis after generating NUM of trees. Default value: 1MLN trees.   
         __-g_stop_t NUM__  - To change the threshold  
         __-g_stop_t 0__  - To turn off the rule  

- __Rule 2:__ stop the analysis after NUM of intermediate trees were visited. Default value: 10MLN trees.       
         __-g_stop_i NUM__  - To change the threshold  
         __-g_stop_i 0__  - To turn off the rule  

- __Rule 3:__ stop the analysis after NUM hours (CPU time). Default: 168 hours (7 days).   
         __-g_stop_h NUM__  - To change the threshold  
         __-g_stop_h 0__  - To turn off the rule  

- To turn off all stopping rules use:  
         __-g_non_stop__ 


### Output
-----
If no other option is used, the ouput is just the log file. The following information is provided about the analysis:

    INFORMATION about input dataset:  
    Number of taxa: NUM  
    Number of partitions: NUM  
    Number of special taxa (row sum = 1): NUM  
    % of missing entries in pr_ab_matrix: NUM  
    Number of taxa on initial tree: NUM  
    Number of taxa to be inserted: NUM  

    SUMMARY:  
    Number of trees on stand: NUM  
    Number of intermediated trees visited: NUM  
    Number of dead ends encountered: NUM  


Here are additional options to output more information:

__-g_print__    - Write all generated trees into a file.  
>WARNING: There might be millions of trees! Use the next option to set the limit on the number of trees to be written into a file.  

__-g_print_lim NUM__  - Limit on the number of trees to be written to a file  
>NOTE: The program will continue generation, but will stop writing trees.  
>WARNING: Do not use summary statistics just on a fraction of trees! Due to the way the generation is performed, consecutive trees have more similar topologies. Therefore, the summary statistics for a fraction of generated trees might be misleading.

__-g_print_induced__  - Write induced partition subtrees  
__-g_print_m__ - Write corresponding presence-absence matrix  


### Additional Analyses
-----

* __Alternative approach: alternative initial trees__   
Due to complexity there might be cases, when Gentrius fails to generate any tree in reasonable time. In this case one can try generating trees using alternative approach: remove NUM of leaves, generate trees by re-inserting these removed leaves. This technique cannot guarantee generating all trees from the stand (i.e. trees compatible with a set of subtrees). However, by using this approach on simulated datasets it was possible to generate millions of trees for complex examples providing a lower bound on stand size.

    __-g_rm_leaves NUM__   - To use alternative approach

* __Test set of trees__   
It is also possible to check, if some trees from a query set have identical set of subtrees as a representative tree (if they belong to the same stand). 

    __-g_query FILE__  - To perform check for a set of trees


User support
------------
If you have further questions, feedback, feature requests, and bug reports, please sign up the following Google group (if not done yet) and post a topic to the 

<https://groups.google.com/d/forum/iqtree>

_The average response time is one working day._

Citations
---------
Preprint is available from bioRxiv
* O. Chernomor, C. Elgert, A. von Haeseler (2023) Identifying equally scoring trees in phylogenomics with incomplete data using Gentrius
