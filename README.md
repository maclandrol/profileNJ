python-tree-processing
======================

[ete2](https://pythonhosted.org/ete2/index.html) TreeNode class extension and a bunch of script/class to process phylogeny tree in python.

Testing scripts can be found in the test/ folder.


### TreeClass

TreeClass is a python class derived from the TreeNode class of the ete2 package.
TreeClass add additional specific function and keep trace of duplication, speciation and lost at each node.

run **pydoc TreeClass.py** on your terminal for documentation.

### TreeUtils 

TreeUtils offer several functions related to phylogeny tree. The resulting tree are returned as a TreeClass object

+ **fetch_ensembl_genetree_by_id** : fetch an ensembl tree using the tree id.
+ **fetch_ensembl_genetree_by_member** : fetch an ensembl tree using a member id.
+ **lcaMapping** : Map a genetree to a specietree.
+ **reconcile** : Reconcile the genetree with the specietree
+ **getTreeFromPhyloxml** : Extract tree/list of tree from a phyloxml file

### NCBI_tree_of_life

*see specie.tar.gz file*

Script to reconstruct the tree of life using the ncbi taxonomy. The current newick file (**tree.nw**) is obtained with the latest ncbi taxonomy release.
