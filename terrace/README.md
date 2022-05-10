Gentrius - Generating Trees from Incomplete Unrooted Subtrees
=======

Introduction
------------
In phylogenomics one infers a species-tree for a group of species using genetic information from multiple genes. Supertree methods typically infer gene trees separately and then summarize them on a single topology. While supermatrix approaches concatenate gene alignments and infer trees from the resulting alignment, typically assuming a partition model.

Both methods are affected by missing data, i.e. when some species do not have sequences for some genes. Identifying, whether incomplete gene trees are compatible, i.e. can be represented on a single tree topology, is computationally hard. Furthermore, depending on the approach, amalgamating incomplete gene trees can lead to different species trees in supertree methods. On the other hand, sparse concatenated alignments in supermatrix approaches might result into huge collections of trees with identical scores, termed Phylogenetic Terraces.

Gentrius is a deterministic algorithm that given a set of incomplete unrooted subtrees identifies, if they are compatible, and, if yes, generates a collection of all corresponding compatible species trees. The same algorithm tackles the problem of enumerating trees on the same Phylogenetic Terrace. Hence, it provides means to ascertain existence (supertree) and uniqueness (supertree and supermatrix) of species trees in the presence of missing data.


Analysis with Gentrius
-----------------------
### Basic Input
-----
There are several ways and inputs one can use to start the analysis. Here are the examples of command lines:

1.  __iqtree2 -gentrius <file 1> -pr_ab_matrix <file 2>__
2.  __iqtree2 -gentrius <file 1> -s <file 3> -p <file 4>__
3.  __iqtree2 -gentrius <file 5>__

    - File 1: a species-tree to be analysed

    - File 2: a presence-absence matrix of loci coverage. First line is a number of species and the number of loci. Other lines contain species name followed by 0's and 1's (with space separator) to indicate presence/absence of gene sequences for corresponding loci. 

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
As the underlying problem is computationally hard, we use three different rules to stop the analysis. The thresholds can be changed by the corresponding options.

- __Rule 1:__ stop the analysis after generating NUM of species-trees. Default value: 1MLN trees.   
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
    % of missing entries in supermatrix: NUM  
    Number of taxa on initial tree: NUM  
    Number of taxa to be inserted: NUM  

    SUMMARY:  
    Number of trees on terrace: NUM  
    Number of intermediated trees visited: NUM  
    Number of dead ends encountered: NUM  


Here are additional options to output more information:

__-g_print__    - Write all generated species-trees into a file.  
>WARNING: There might be millions of trees! Use the next option to set the limit on the number of trees to be written into a file.  

__-g_print_lim NUM__  - Limit on the number of species-trees to be written to a file  
>NOTE: The program will continue generation, but will stop writing trees.  
>WARNING: Do not use summary statistics just on a fraction of trees! Due to the way the generation is performed, consecutive trees have more similar topologies. Therefore, the summary statistics for a fraction of generated trees might be misleading.

__-g_print_induced__  - Write induced partition subtrees  
__-g_print_m__ - Write corresponding presence-absence matrix  


### Additional Analyses
-----

* __Backward approach__   
Due to complexity there might be cases, when Gentrius fails to generate any species-tree. In this case one can try to generate species-trees using backward approach: remove NUM of leaves, generate species-trees by re-inserting these removed leaves. This technique cannot guarantee generating all species-trees compatible with a set of subtrees. However, by using this approach on simulated datasets it was possible to generate millions of species-trees for complex examples.

    __-g_rm_leaves NUM__   - To use backward approach

* __Test set of species-trees__   
It is also possible to check, if some species-trees from a query set have identical set of subtrees as a representative species-tree (if they belong to the same phylogenetic terrace). 

    __-g_query FILE__  - To perform check for a set of species-trees


User support
------------
If you have further questions, feedback, feature requests, and bug reports, please sign up the following Google group (if not done yet) and post a topic to the 

<https://groups.google.com/d/forum/iqtree>

_The average response time is one working day._

Citations
---------

The manuscript is currently in preparation:
* O. Chernomor, A. von Haeseler (in prep.) Gentrius: generating trees from incomplete unrooted subtrees.
