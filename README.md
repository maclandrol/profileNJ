python-tree-processing
=======================

Utility package for tree processing written in python and based on the [ete2](https://pythonhosted.org/ete2/index.html) toolkit.


## TreeClass

Bases from the TreeNode class of the ete2 package, TreeClass is a tree representation class. A tree consists of a collection of TreeClass instances connected in a hierarchical way. A TreeClass object can be loaded from the New Hampshire Newick format (newick).
TreeClass add specific functions for tree processing not present in the ete2 TreeNode.

run pydoc for minimum documentation.

## TreeUtils 

TreeUtils offer several static functions related to phylogeny tree. With You can fetch ensembl genetree and reconcile a gene tree to its species tree.


+ **fetch_ensembl_genetree_by_id** : fetch an ensembl tree using the tree id.
+ **fetch_ensembl_genetree_by_member** : fetch an ensembl tree using a member id.
+ **lcaMapping** : Map a genetree to a specietree.
+ **reconcile** : Reconcile the genetree with the specietree

## ClusterUtils 

ClusterUtils is an implementation of UPGMA (Unweighted Pair Group Method with Arithmetic Mean) and NJ (Neighbor-Joining), two clustering distance-based method for tree construction.


## Multipolysolver

Multipolysolver is a module for polytomy correction in genetree using the duplication-lost criterion. Mutltipolysolver output all the binary solution for a non-binary gene tree that minimize the duplication-lost score. If the input tree is considered unrooted, Multipolysolver can test every possible root and return the binary tree with the lowest duplication-lost or all rooted binary tree.


## PolytomySolver

PolytomySolver is the executable version of **Multipolysolver**.

#### Command line arguments

+  *-h, --help*            
        Show the help message and exit.

+  *-s SPECIETREE, --sFile SPECIETREE* 
        Name of the file containing the species newick tree (default: None).

+  *S SMAP, --sMap SMAP*
        Gene to species map. Use the standard format. (default: None)

+  *-g GENETREE, --gFile GENETREE*
        Name of the file containing the gene newick tree (default: None).

+  *-d DISTFILE, --dist DISTFILE*
        Name of the file containing the distances between each pair of genes (The gene set should be the same for the leaf set of the genetree). (default: None)

+  *-o OUTFILE, --output OUTFILE*
        Output file with the corrected tree. The genetree is printed on stdout if omitted. (default: None)

+  *-gl GLINE, --gLine GLINE*
        Index of the line in the gene tree file that corresponds to the current gene tree. Indexing start with "1" (default: 1)

+  *-r {none,all,best}, --reroot {none,all,best}*
        Enable/Disable root mode. (default: none)
          	* none: disable reroot mode, correct the input polytomies and return the result.
          	* all: enable reroot mode, reroot the genetree at each node and return all polytomy corrected version for each rooted tree.
          	* best: enable reroot mode, rerrot the genetree at each node and return all polytomy corrected version for the rooted tree with the smallest Dup-Lost score (First one if not unique).
                        
+  *-nf, --nnflag*         
        Treat negative distances as large distances (default:False).

+  *--sep GENE_SEP*
        Gene-Specie separator for each leaf name in the genetree. PolytomySolver will guess by default in a very limited list of special character. (default: None) 
        **';;' is not a valid separator for the newick format! IF YOUR SEPARATOR IS ";;", DON'T USE THIS FLAG. THE PROGRAM WILL AUTOMATICALLY GUESS. ** 

+  *--mValue MVAL*        
        Set largest value in the distance matrix. Entries on the main diagonal and negative values will be replaced by mValue. (default: 1e+305)

+  *-c {nj,upgma}, --cluster {nj,upgma}*
        Set the clustering methods. (default: nj) 
            * upgma: UPGMA (Unweighted Pair Group Method with Arithmetic Mean) clustering algo.
            * nj: neighbor joining clustering method, (slower).


+  *--slimit SOL_LIMIT*    
        Set the max number of solution per genetree. Possible values are -1 to return all the solution or n, n>0 for a specific number of solution. Setting this argument to -1 is computationally expensive. (default: 30)

+  *--plimit PATH_LIMIT*   
        Set the max number of solution for each polytomy in the genetree. Possible values are -1 to explore all the solution or n, n>0 for a specific number of  solution. Setting this argument to -1 is also computationally expensive. (default: -1)

+  *-v*                    
        Output verbosity (default: False)

+  *--cap*                 
        Capitalize the species name of the genetree leaves to  match each species. Almost all functions are case sensitive. (default: False)

#### File formats
  see [polytomy-solver-distance] (https://github.com/UdeM-LBIT/polytomy-solver-distance#file-formats)


## NCBI_tree_of_life 

*see specie.tar.gz file*

Script to reconstruct the tree of life using the ncbi taxonomy. The current newick file (**tree.nw**) is obtained with the latest ncbi taxonomy release.
