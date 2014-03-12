from TreeClass import TreeClass
from ete2 import PhyloTree

from TreeUtils import *

genetree = TreeClass("((((a,a,b)node1,a)node2, (b,a)node3)node_6, ((d,e)node_4,c)node_5)root;", format=1)
specietree = TreeClass("((A,B)cat, (C,D,E)dog)mammifere;", format=1)

#print the geneTreee and the specieTree
print genetree.get_ascii(attributes=["species", "type"], show_internal=True)
print specietree.get_ascii(show_internal=True)

# Use a map of gene to specie for reconcilliation
s_map={"a":"A", "b":"B", "c":"C", "d":"D","e":"E"}

#set the specie of the gene, using the map. There are a lot of way to do that, check the pydoc of set_species
genetree.set_species(speciesMap=s_map)

#lcaMapping between genetree and specietree, the return dict will be use for reconciliation
mapping=lcaMapping(genetree, specietree)

#we are doing reconcilliation here, "yes", is to infere gene lost, default="no", see method doc
reconcile(genetree, mapping, "yes")
#After reconciliation
print genetree.get_ascii(attributes=["type", "species"], show_internal=True)

#restrict genetree to those species
genetree.restrictToSpecies(species=["A", "E", "D"])

#a simple way to do this to a specie tree, is to get a list of node with the wanted specie name, using the search_nodes function
#and then use the original prune function of TreeNode
print genetree.get_ascii(attributes=["type", "species"], show_internal=True)



#an example of fetching a genetree from ensembl with its ensembl ID
ensemblTree=fetch_ensembl_genetree_by_id(treeID="ENSGT00390000003602", output="phyloxml", aligned=1, sequence="protein")
#save genetree and all the feature to a file
ensemblTree.write(features=[],outfile="test.txt", format=1)
print ensemblTree

#an example of fetching a genetree from ensembl with a member ID
ensemblT=fetch_ensembl_genetree_by_member(memberID="ENSG00000157764")

#testing show function
ensemblT.show()

#This is another test of reconciliation using events on a specific node
gene_tree_nw = '((Dme_001,Dme_002),(((Cfa_001,Mms_001),((Hsa_001,Ptr_001),Mmu_001)),(Ptr_002,(Hsa_002,Mmu_002))));'
t=TreeClass(gene_tree_nw)
print t
t.set_species()
matches = t.search_nodes(name="Hsa_001")
human_seq = matches[0]
# Obtains its evolutionary events
events = human_seq.get_my_evol_events()
# Print its orthology and paralogy relationships
print "Events detected that involve Hsa_001:"
for ev in events:
    if ev.etype == "S":
        print '   ORTHOLOGY RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs)
    elif ev.etype == "D":
        print '   PARALOGY RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs)

# # Loads an example tree
# nw = """
# ((Dme_001,Dme_002),(((Cfa_001,Mms_001),((Hsa_001,Ptr_001),Mmu_001)),
# (Ptr_002,(Hsa_002,Mmu_002))));
# """
# t = PhyloTree(nw)
# print t
# #                    /-Dme_001
# #          /--------|
# #         |          \-Dme_002
# #         |
# #         |                              /-Cfa_001
# #         |                    /--------|
# #---------|                   |          \-Mms_001
# #         |          /--------|
# #         |         |         |                    /-Hsa_001
# #         |         |         |          /--------|
# #         |         |          \--------|          \-Ptr_001
# #          \--------|                   |
# #                   |                    \-Mmu_001
# #                   |
# #                   |          /-Ptr_002
# #                    \--------|
# #                             |          /-Hsa_002
# #                              \--------|
# #                                        \-Mmu_002
# #
# # To obtain all the evolutionary events involving a given leaf node we
# # use get_my_evol_events method
# matches = t.search_nodes(name="Hsa_001")
# human_seq = matches[0]
# # Obtains its evolutionary events
# events = human_seq.get_my_evol_events()
# # Print its orthology and paralogy relationships
# print "Events detected that involve Hsa_001:"
# for ev in events:
#     if ev.etype == "S":
#         print '   ORTHOLOGY RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs)
#     elif ev.etype == "D":
#         print '   PARALOGY RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs)
