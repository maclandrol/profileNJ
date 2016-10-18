[![Build Status](https://travis-ci.org/UdeM-LBIT/profileNJ.svg?branch=master)](https://travis-ci.org/UdeM-LBIT/profileNJ)

profileNJ
=======================

Utility package for gene tree correction and reconciliation with a species tree, written in python and based on the [ete3](http://etetoolkit.org/) toolkit.

## Installation 
#### Dependencies
profileNJ depends on the following package:
 - ete3
 - numpy>=1.8
 - lxml
 - PyQt4 (Not mandatory but required by `reconcile` to print trees)
 
[PyQt4](http://pyqt.sourceforge.net/Docs/PyQt4/installation.html) can be tricky to install, I recommand using either your distribution version or install it with conda (`conda install -c anaconda pyqt=4`). 

To install profileNJ, download the package on github and install with :

`python setup.py install`

or use pip : `pip install profileNJ`

You may need sudo privileges. You can also install a local version by using the *'--user'* flag. 


## Minimum Documentation

### profileNJ

_profileNJ_ correct genetree by contracting weak branches and resolving them to have binary trees with a minimum reconciliation cost to their specietree. _profileNJ_ use NJ in order to keep sequence information as much as possible an can output multiple solutions. If the input tree is considered unrooted, profileNJ can test every possible root and return the binary tree with the lowest duplication-lost or all rooted binary tree. A detailed description of the algorithm can be found in our paper : 

        Noutahi E, Semeria M, Lafond M, Seguin J, Boussau B, et al. (2016) Efficient Gene Tree Correction Guided by Genome Evolution. PLoS ONE 11(8): e0159559. doi: 10.1371/journal.pone.0159559


+  *-h, --help*            
        Show the help message and exit.
+  *-s SPECIETREE, --sFile SPECIETREE*          
       Name of the file containing the species newick tree (default: None).
+  *-S SMAP, --sMap SMAP*        
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
          	- **none** : disable reroot mode, correct the input genetree and return the result.        
          	- **all** : enable reroot mode, reroot the genetree at each node and return all solution for each rooted tree.   
          	- **best** : enable reroot mode, reroot the genetree at each node and return any solution that minimize the reconciliation score.
+  *-nf, --nnflag*         
        Treat negative distances as large distances (default:False).
+  *--sep GENE_SEP*
        Gene-Specie separator for each leaf name in the genetree. PolytomySolver will guess by default in a very limited list of special character. (default: None) 
        **';;' is not a valid separator for the newick format! IF YOUR SEPARATOR IS ";;", DON'T USE THIS FLAG. THE PROGRAM WILL AUTOMATICALLY GUESS. ** 
+  *--spos {prefix,postfix}*
        The position of the specie name according to the separator. Supported option are prefix and postfix (default: prefix)
+  *--mValue MVAL*        
        Set largest value in the distance matrix. Entries on the main diagonal and negative values will be replaced by mValue. (default: 1e+305)
+  *-c {nj,upgma}, --cluster {nj,upgma}*
        Set the clustering methods. (default: nj)       
            - **upgma** : UPGMA (Unweighted Pair Group Method with Arithmetic Mean) clustering algo.    
            - **nj** : neighbor joining clustering method, (slower).
+  *--slimit SOL_LIMIT*    
        Set the max number of solution per genetree. Possible values are -1 to return all the solution or n, n>0 for a specific number of solution.
        Setting this argument to -1 is computationally expensive. (default: 30)
+  *--plimit PATH_LIMIT*   
        Set the max number of solution for each polytomy in the genetree. Possible values are -1 to explore all the solution or n, n>0 for a specific number of  solution.
        Setting this argument to -1 is also computationally expensive. (default: -1)
+  *-v*                    
        Output verbosity (default: False)
+  *--cap*                 
        Capitalize the species name of the genetree leaves to  match each species. Almost all functions are case sensitive. (default: False)
+  *--batch*  
        Use this flag to enable batch mode. In batch mode, gLine value is discarded, --dist should be a file whose line link to the distance matrix file of the genetree at the same line number in your genetree file
+  *--parallelize*  
        Use parallelization (default False)
+  *--firstbest*  
        Only output solution for the first root with the best dl score encountered
+  *--cost D L*  
        Change the cost of duplications (D) and losses(L). 
        D L : 2 float values, duplication and loss cost in this order (default:  D=1 and L=1 )
+  *--seuil*  
        Branch contraction threshold, when the tree is binary. Use only when the tree is binary.


#### File formats
  see [polytomy-solver-distance] (https://github.com/UdeM-LBIT/polytomy-solver-distance#file-formats)


### reconcile

_reconcile_  compute  and output the reconcilied gene tree, and it's cost, between a binary genetree and a binary species tree. Two mode are possible : in the run mode, you can compute the reconcilied gene tree, whereas in the smap mode, _reconcile_ return an automatic map between the genes in the gene trees and the species.

optional arguments:
+  *-h, --help*
        show this help message and exit
+  *-s SPECIETREE, --sFile SPECIETREE*
        Specie tree in newick format
+  *-S SMAP, --sMap SMAP*
        Gene to species map. Use the standard format.
+  *-g GENETREE, --gFile GENETREE*
        Gene tree in newick format.
+  *--sep GENE_SEP*
        Gene-Specie separator for each leaf name in the genetree. The program will guess by default. But you should provide it
+  *--spos {prefix,postfix}*
        The position of the specie name according to the separator. Supported option are prefix and postfix (default: prefix)
+  *--cap*
        Capitalize the species name of the genetree leaves to
        match each species. Almost all functions are case sensitive.
+  *--display_losses*
        Display losses in the reconciliation
+  *--output OUTPUT, --out OUTPUT*
        Output an image of the reconciliated tree.
+  *--outform OUTFORM*
        Accepted format are svg|pdf|png
+  *--export_orthoxml*
        Export reconciliated tree to export_orthoxml. Losses are not exported
+  *--show_branch_tag*
        Show branch length and support
+  *--verbose, -v*
        Verbosity
+  *--reroot, -r*
        Reroot in order to display all root reconciliation
        result. Not Available in batch mode !!
+  *--debug*
        Debugging ( test, will be removed )
+  *--batch*
        Batch mode, use profileNJ file output


### polytomySolver

_polytomySolver_ is a new algorithm for resolving gene trees with polytomies in linear time. _polytomySolver_ support both unit and weighted duplication and loss cost. It's an improved version of the quadratic algorithm described by Lafond and al. in 2012 (M. Lafond, K.M. Swenson, and N. El-Mabrouk. An optimal reconciliation algorithm for gene trees with polytomies. In LNCS, volume 7534 of WABI, pages 106-122, 2012.), using the compressed species tree idea of Zheng and Zhang (Y. Zheng and L. Zhang. Reconciliation with non-binary gene trees revisited. In Lecture Notes in Computer Science, volume 8394, pages 418-432, 2014. Proceedings of RECOMB.). _polytomySolver_ is faster than Notung and thus can be used on large trees.

+  *-h, --help*
        show this help message and exit
+  *-s SPECNW, --spectree SPECNW*
        Name of the file containing the species newick tree.
+  *-S SMAP, --sMap SMAP*
        Gene to species map. Use the standard format.
+  *-g GENENW, --genetree GENENW*
        Name of the file containing the gene newick tree.
+  *--sep GENE_SEP*
        Specify a gene separator if you're are not using a smap
+  *--losscost LOSSCOST*
        Specify the losses cost
+  *--dupcost DUPCOST*
        Specify the duplication cost
+  *--nsol NSOL*
        Number of solution to output, Our implementation of Notung algorithms and both linear and dynamic versions of the Zheng and Zhang's algorithm only output one solution
+  *--spos SPOS*
        Gene position when you have specified a separator. Possible values are "prefix" and "postfix". Default value is prefix.
+  *-o OUTFILE, --output OUTFILE*
        Name of your output files with the corrected tree. The resolutions are printed on stdout if omitted.
+  *--mode {psolver,linzz,dynzz,notung,dynzz2}*
        Algorithm to use. psolver is the default algorithm and correspond to _polytomySolver_. "linzz" and "dynzz" are respectively the linear and dynamic version of the Zheng and Zhang's algorithm. 
+  *--showcost*
        Use this to show the reconciliated cost at the end. By default, only the resolved tree is shown


## Reusable modules

### TreeClass

Bases from the TreeNode class of the ete package, TreeClass is a tree representation class. A tree consists of a collection of TreeClass instances connected in a hierarchical way. A TreeClass object can be loaded from the New Hampshire Newick format (newick).
TreeClass add specific functions for tree processing not present in ete's TreeNode.

run pydoc for minimum documentation.

### TreeUtils 

TreeUtils offer several static functions related to phylogeny tree. With You can fetch ensembl genetree and reconcile a gene tree to its species tree.


+ **fetch_ensembl_genetree_by_id** : fetch an ensembl tree using the tree id.
+ **fetch_ensembl_genetree_by_member** : fetch an ensembl tree using a member id.
+ **lcaMapping** : Map a genetree to a specietree.
+ **reconcile** : Reconcile the genetree with the specietree

### ClusterUtils 

ClusterUtils is an implementation of UPGMA (Unweighted Pair Group Method with Arithmetic Mean) and NJ (Neighbor-Joining), two clustering distance-based method for tree construction.

## NCBI_tree_of_life 

*see SPECIES_TREE*

Script to reconstruct the tree of life using the ncbi taxonomy. The current newick file (**tree.nw**) is obtained with the latest ncbi taxonomy release.


## How to cite.

If you use profileNJ or polytomySolver, please cite the following papers:

+ Noutahi E, Semeria M, Lafond M, Seguin J, Boussau B, et al. (2016) Efficient Gene Tree Correction Guided by Genome Evolution. PLoS ONE 11(8): e0159559. doi: 10.1371/journal.pone.0159559

+ Lafond M, Noutahi E and EL-Mabrouk N. (2016) Efficient Non-Binary Gene Tree Resolution with Weighted Reconciliation Cost. 27th Annual Symposium on Combinatorial Pattern Matching (CPM 2016)
