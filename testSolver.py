#Execute Multpolysolver on a no-star genetree with rerooting
from Multipolysolver import *
from TreeLib import *

corrected_tree= TreeClass("testinput/famille_1.corrected")
c_tree= TreeClass("testinput/famille_1.corrected")
corrected_tree.set_species(sep='__', capitalize=True, pos="postfix")
TreeUtils.reset_node_name(corrected_tree, '__')
TreeUtils.reset_node_name(c_tree, '__')

taille=len(corrected_tree.get_descendants())
arbre_species= TreeClass('testinput/speciestree.newick')
mapping=TreeUtils.lcaMapping(corrected_tree, arbre_species)
#print corrected_tree.get_ascii(attributes=['name', 'species'], show_internal=False)

TreeUtils.reconcile(corrected_tree, mapping, "yes")
c_score=TreeUtils.ComputeDupLostScore(corrected_tree)
#"""
origenetree, specietree, distance_matrix, node_order = TreeUtils.polySolverPreprocessing('testinput/nostar_genetree.tree', 'testinput/speciestree.newick', 'testinput/nostar_famille_1.dist', capitalize=True, gene_sep = None)
tree_list=[origenetree]
tree_list.extend(origenetree.reroot())
for genetree in tree_list:
	polysolution = solvePolytomy(genetree, specietree, distance_matrix, node_order, poly_limit=1, method='nj')
	print "Nombre de solution: ", len(polysolution)
	#we are doing reconcilliation here, "yes", is to infere gene lost, default="no", see method doc
	#specietree=TreeClass("(((a,b)e,g)h,(c,d)f)i;", format=1)
	#print specietree.get_ascii(show_internal=True)
	for tree in polysolution:
		#tree=a.copy("newick-extended")
		#tree=TreeClass(a.write(features=["species","name", "cost","dist" ], format=1))
		m_score_sum=0
		#print tree.get_ascii(attributes=['name', 'species'], show_internal=False)
		for node in tree.traverse():
			if node.has_feature('type'):
				node.del_feature('type')
			if node.has_feature('species') and node.is_internal():
				node.del_feature('species')
			if node.has_feature('cost'):
				m_score_sum+=float(node.cost)
		print "***************************************************"
		rf, max_rf,common_leaves, parts_t1, parts_t2= tree.robinson_foulds(c_tree)
		#TreeUtils.customTreeCompare(TreeClass('nostar_genetree.tree'), c_tree, tree)
		print '**tree has', len(tree.get_descendants())+1, "nodes and corrected tree has", taille+1 , 'nodes'
		print tree.get_ascii(attributes=['name', 'species'], show_internal=False)
		#tree.set_species(sep='_', capitalize=False, pos="prefix")
		mapping=TreeUtils.lcaMapping(tree, arbre_species)
		TreeUtils.reconcile(tree, mapping, "yes")
		score=TreeUtils.ComputeDupLostScore(tree)
		#TreeUtils.CleanFeatures(tree, ['type'])
		print "**m_score= ", m_score_sum, 'r_score= ', score, "c_score= ",c_score, 'reconcile tree node= ', len(tree.get_descendants())+1, " leaves: ", len(tree), " reconcile corrected tree node= ", len(corrected_tree.get_descendants())+1, "leaves: ", len(corrected_tree)
		print "**rf=", rf, "max_rf=", max_rf
		print "**tree still has polytomies= ", tree.has_polytomies()
