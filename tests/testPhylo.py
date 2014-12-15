# Loads a gene tree and its corresponding species tree. Note that
# species names in sptree are the 3 firs letters of leaf nodes in
# genetree.
from pprint import pprint
from TreeLib import *
from ete2 import PhyloTree
from time import time
import timeit

def standard():
	gene_tree_nw = '((Dme_001,Dme_002),(((Cfa_001,Mms_001),((Hsa_001,Ptr_001),Mmu_001)),(Ptr_002,(Hsa_002,Mmu_002))));'
	species_tree_nw = "((((Hsa, Ptr), Mmu), (Mms, Cfa)), Dme);"
	genetree = PhyloTree(gene_tree_nw)
	sptree = PhyloTree(species_tree_nw)
	print "\n--->Original Tree\n", genetree
	# Let's reconcile our genetree with the species tree
	recon_tree, events = genetree.reconcile(sptree)
	print "Tree after reconciliation\n", recon_tree.get_ascii(attributes=[], show_internal=False)
	recon_tree.describe()

def MyReconcile(perte="no"):
	gene_tree_nw = '((Dme_001,Dme_002),(((Cfa_001,Mms_001),((Hsa_001,Ptr_001),Mmu_001)),(Ptr_002,(Hsa_002,Mmu_002))));'
	species_tree_nw = "((((Hsa, Ptr), Mmu), (Mms, Cfa)), Dme);"
	genetree = TreeClass(gene_tree_nw)
	specietree = TreeClass(species_tree_nw)
	genetree.set_species(sep='_', pos="prefix")
	print "\n--->Original Tree\n", genetree
	mapping=TreeUtils.lcaMapping(genetree, specietree)
	TreeUtils.reconcile(genetree, mapping, perte)
	print "Tree after reconciliation\n", genetree.get_ascii(attributes=["species", "type"], show_internal=False)
	n_leaf=0
	for leaf in genetree:
		n_leaf=n_leaf+len(set(leaf.get_species()))
	print n_leaf

def get_species_name(node_name_string):
	return node_name_string.split('_')[-1]

if __name__ == '__main__':

	#python -mtimeit -s 'import testPhylo' 'testPhylo.MyReconcile("yes")'
	#python -mtimeit -s 'import testPhylo' 'testPhylo.MyReconcile("no")'
	#python -mtimeit -s 'import testPhylo' 'testPhylo.standard()'

	print "****TreeUtils Reconciliation\n"
	MyReconcile("yes")
	print "n\n****PhyloTree Reconciliation"
	standard()
	print "\n\n****Ensembl test\n"
	ensemblTree=TreeUtils.fetch_ensembl_genetree_by_id(treeID="ENSGT00390000003602", nh_format="display_label_composite")
	phylo_genetree=PhyloTree(ensemblTree.write(features=[],format_root_node=True))
	phylo_genetree.set_species_naming_function(get_species_name)

	ensemblTree.set_species(pos="postfix")
	#print ensemblTree.get_ascii(show_internal=False, attributes=list(ensemblTree.get_all_features()))
	specie_list=[]
	for leaf in ensemblTree:
		specie_list.extend(leaf.get_species())
	
	specie_list=list(set(specie_list))
	specieTree= TreeClass()
	specieTree.populate(len(specie_list),names_library=specie_list)
	phylo_specietree= PhyloTree(specieTree.write(features=[],format_root_node=True))
	print "\n---> Ensembl Tree", ensemblTree
	#Reconciliation
	start=time()
	mapping=TreeUtils.lcaMapping(ensemblTree, specieTree)
	TreeUtils.reconcile(ensemblTree, mapping, "yes")
	print "-->TreeUtils reconciliation Time :", time()-start
	#print ensemblTree.get_ascii(attributes=["type"], show_internal=True)
	start=time()
	recon_tree, events = phylo_genetree.reconcile(phylo_specietree)
	print "-->PhyloTree reconciliation Time :", time()-start
	phylo_recontree=TreeClass(recon_tree.write(features=[],format_root_node=True))


	#rf, max_rf, common_leaves, parts_t1, parts_t2=	phylo_recontree.robinson_foulds(ensemblTree)
	#print "RF distance is %s over a total of %s" %(rf, max_rf)
	#print "Partitions in tree2 that were not found in tree1:", parts_t1 - parts_t2
	#print "Partitions in tree1 that were not found in tree2:", parts_t2 - parts_t1

	#Ensembl sequence retrieving test
	ensemblTree=TreeUtils.fetch_ensembl_genetree_by_id(treeID="ENSGT00390000003602", sequence='cdna', aligned=1, output="phyloxml")
	print ensemblTree.get_all_features()
	print ensemblTree.get_ascii(show_internal=False, attributes=['name', 'alt_name', 'support'])
	print ensemblTree.write()
