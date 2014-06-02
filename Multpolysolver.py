"""Python Polysolver"""

from TreeClass import TreeClass
import TreeUtils
import numpy
import collections
import random
from pprint import pprint
import copy
import Upgma
import Utils

"""
Gene matrix are represented by a numpy array
"""

DUP='d'
LOST='l'
SPEC='s'
PARTIAL_RESOLUTION_ITERATOR=1
dupcost,lostcost=1,1
numpy.set_printoptions(threshold='nan')

def polySolver(genetree, specietree, gene_matrix, node_order, limit=-1, verbose=0):
	"""This assume we that we are using the correct specietree for this genetree
	the specie tree root is the latest common ancestor of all the specie in genetree"""
	count=getSpecieCount(genetree) # number of specie in the genetree
	max_y=max(count.values())+1
	# assigning a correspondance between each row and a node
	# the assignment is done in level order

	polytomy_specie_set, row_node_corr= findMaxX(genetree, specietree)

	max_x=len(polytomy_specie_set)
	cost_table= numpy.zeros((max_x,max_y)) # cost cost_table to fill
	path_table=numpy.ndarray((max_x,max_y), dtype='object') # table to save the possible path

	# fill the cost_table and the path_table
	for n in range(0,max_x):
		node=row_node_corr[n]
		zeropos=count[node.name]-1  
		# We have zeropos when the number of node from a specie is the same as the column number
		# Fill the table, using the next/previous case cost
		# The node is a leaf, just fill with dupcost and lostcost
		if(node.is_leaf()):
			# find the column with a cost of zero (zeropos) and fill the table
			# according to this position
			# by default, all the position in the table are 0
			i=zeropos-1
			while(i>=0):
				cost_table[n,i]=cost_table[n,i+1]+dupcost
				path_table[n,i]=DUP
				i-=1
			i=zeropos+1
			while(i<max_y):
				cost_table[n,i]=cost_table[n,i-1]+lostcost
				path_table[n,i]=LOST
				i+=1
			# We should take into account the special case here
		# Here we have an internal node (not a leaf in the genetree)
		else:
			l_child_id=[key for key, value in row_node_corr.iteritems()  if(value==node.get_child_at(0))][0]
			r_child_id=[key for key, value in row_node_corr.iteritems()  if(value==node.get_child_at(1))][0]
			# Fill the table using only the speciation cost(sum of the children's cost of this node)
			for k in range(0, max_y):
				if(k>zeropos):
					cost_table[n,k]= cost_table[l_child_id,k-zeropos-1]+cost_table[r_child_id, k-zeropos-1]
					path_table[n,k]=SPEC
				else:
					cost_table[n,k]=numpy.inf

			# Find all the min score position and try to minimize the score of its
			# neighborhood by lost/dup cost
			minpos= numpy.where(cost_table[n,:]== cost_table[n,:].min())
			for pos in minpos[0]:
				i=pos-1
				while(i>=0):
					if(cost_table[n,i]==cost_table[n,i+1]+lostcost):
						if DUP not in path_table[n,i]: path_table[n,i]+=DUP 
					elif (cost_table[n,i]>cost_table[n,i+1]+lostcost):
						cost_table[n,i]=min(cost_table[n,i+1]+lostcost, cost_table[n,i])
						path_table[n,i]=DUP
					i-=1

				i=pos+1
				while(i<max_y):
					if(cost_table[n,i]==cost_table[n,i-1]+lostcost):
						if LOST not in path_table[n,i]: path_table[n,i]+=LOST 
					elif (cost_table[n,i]>cost_table[n,i-1]+lostcost):
						cost_table[n,i]=min(cost_table[n,i-1]+lostcost, cost_table[n,i])
						path_table[n,i]=LOST
					i+=1

			# Should also consider a special case for filling the table here


	# find the shape of the cost_table
	xsize,ysize=cost_table.shape

	# Use the list of current leave in the genetree to find all the possible path
	# in the tree construction
	curr_leaves=[x for x in genetree.get_children()] 

	#t=treeConstruct(cost_table, path_table, row_node_corr, xsize-1, 0, genetree, curr_leaves)['topo']
	#for i in t:
	#	print i.get_ascii(show_internal=True, attributes=['name','type'])

	# find the path from the position (x_size-1, 0) to the leaves
	paths= findPathFromTable(path_table, row_node_corr, count, xsize-1, 0);
	solution=[]
	mapGene= getReversedMap(genetree, specietree)

	if(verbose):
		print "Matrix M: \n"
		print cost_table
		print
		print "Path table for tree construction: \n"
		print path_table
		print
		print "Correspondance: \n"
		pprint(row_node_corr)
		print 
		print "Gene Tree:\n"
		print genetree.get_ascii(attributes=['species', 'name'])
		print 
		print "Specie Tree:\n"
		print specietree.get_ascii()
		print "\nNumber of Tree found : ", len(paths), "\n"
	i=1
	for path in paths:
		if(limit>0 and i>limit):
			break
		#pprint(path)
		#print "**Tree ", i
		solution.append(constructFromPath(path, genetree, specietree, mapGene, numpy.copy(gene_matrix), node_order[:], 1e305))
		i+=1

	for t in solution:
		t.add_features(cost=cost_table[xsize-1,0])
	return solution

def findPathFromTable(path_table, row_node_corr, count, xpos, ypos):
	""" Find all the possible path from the lower left case to the leaves"""

	chemin=[]
	#Case 1: current position correspond to a leaf
	if(row_node_corr[xpos].is_leaf() and (ypos<0 or path_table[xpos, ypos] is None)):
		case=row_node_corr[xpos].name+':%i'%(ypos+1)
		chemin.append(case)
	
	#Case 2 : this a internal node
	else:
		#each case can have multiple path
		for c in path_table[xpos, ypos]:

			# we found a speciation
			if c=='s':
				# this is bad for perfomance
				spec_pos_1 = [x for x in row_node_corr.keys() if row_node_corr[xpos].get_child_at(0)==row_node_corr[x]][0]
				spec_pos_2 = [x for x in row_node_corr.keys() if row_node_corr[xpos].get_child_at(1)==row_node_corr[x]][0]
				snode=row_node_corr[xpos].name
				nb_node = count[snode]
				spec_1= findPathFromTable(path_table, row_node_corr, count, spec_pos_1, ypos-nb_node)
				spec_2= findPathFromTable(path_table, row_node_corr, count, spec_pos_2, ypos-nb_node)
				# add all possible path from the children
				for path1 in spec_1:
					for path2 in spec_2:
						chemin.append(",".join([row_node_corr[xpos].name+':%i'%(ypos+1),path1, path2]))

			# we found a duplication
			elif c=='d':
				dup=findPathFromTable(path_table, row_node_corr,count, xpos, ypos+1)
				# add possible path of the case that lead to this duplication
				for path1 in dup:
					chemin.append(",".join([row_node_corr[xpos].name+':%i'%(ypos+1), path1]))

			# instead we found a lost
			elif c=='l':
				lost=findPathFromTable(path_table, row_node_corr, count, xpos, ypos-1)
				# add possible path of the case that lead to this lost
				for path1 in lost:
					chemin.append(",".join([row_node_corr[xpos].name+':%i'%(ypos+1), path1]))
				
	# return a list of path, the path are specially formated string
	return chemin


def constructFromPath(chemin, genetree, specietree, sptree_mapping, gene_matrix, node_order, maxVal, verbose=0):
	"""Construct tree from a path using the upgma method"""
	# get the node order in the path
	node_list= list(reversed(chemin.split(','))) 
	# find the node to show in the tree construction
	node_in_tree = genetree.get_children()
	leaf_list = [x.name for x in specietree.get_leaves()]
	gene_tree_desc_species= genetree.get_descendant_species()
	species_list=[x.name for x in specietree.traverse()]
	gene_root_species= specietree.get_common_ancestor([x for x in gene_tree_desc_species if(x in species_list)])
	#leaf_list = genetree.get_children_species()
	#keep this in case a lost node is in the path
	lost_nodes=[]
		# total number of node
	tot_node = len(node_list)
	deja_vu=[]
	# traverse the list of node in the path and construct the tree
	for indice in range(tot_node): #0 - len(node_list)-1
		
		#find the current node and its position
		[node, pos]= node_list[indice].split(":")
		# find the next node and its position
		[n_node, n_pos]=[None, -1]
		if(indice+1 != tot_node):
			[n_node, n_pos]= node_list[indice+1].split(":")
		pos=int(pos)
		n_pos=int(n_pos)
		
		#list of node which specie is the same as the current node
		node_structs=[n for n in node_in_tree if n.species==node]
		#index to keep in the node order after for upgma
		ind_to_keep=[]

		# simple cas, the node is a leaf
		if(node in leaf_list):
			# we have a duplication here,
			if(n_node==node and n_pos<pos): 
				# ind to keep is empty on purpose, we have to find the good index
				# and remove the row we don't want for the upgma algorithme to work with
				upgma=Upgma.UPGMA_cluster(getMatrix(node_structs, gene_matrix, node_order, ind_to_keep), node_structs, maxVal, 1)
				# find the resulting duplication tree
				dup_tree= upgma[0]
				# set the name of the duplication, i'll try noName next time
				dup_tree.name="-".join([dup_tree.get_child_at(0).name, dup_tree.get_child_at(1).name])
				dup_tree.add_features(type="DUP")
				#find the merged index and remove them from the distance matrice
				merged_index=map(lambda x: ind_to_keep[x], upgma[2])
				dup_tree.add_features(score=gene_matrix[merged_index[1], merged_index[0]])
				dup_tree.add_features(species=node)
				# Update the gene matrix and the node order
				gene_matrix=Upgma.del_row_column(gene_matrix,merged_index, maxVal)
				node_order[merged_index[0]]=dup_tree.name
				node_order.pop(merged_index[1])

				for n in dup_tree.get_children():
					node_in_tree.remove(n)
				node_in_tree.append(dup_tree)

			# the leaf is a lost leaf, add it to the list of lost_nodes
			# which will be checked when we want to construct a internal node with a lost child
			elif(n_node==node and n_pos>pos):
				lost=TreeClass()
				lost.name=node+'_%i'%(n_pos)
				lost.species=node
				lost.add_features(type='LOST')
				lost_nodes.append(lost)


		else: #cas de noeud interne
			if(n_node==node and pos<n_pos):
				"Lost internal node"
				lost=TreeClass()
				lost.name=node+'_%i'%(n_pos)
				lost.species=node
				lost.add_features(type='LOST')
				lost_nodes.append(lost)

			elif(pos>=1 and node not in deja_vu):
				i=pos
				deja_vu.append(node)
				# Case where we have this node n times, n>=1. We should then construct the node all the n times
				i=len([x for x in node_in_tree if x.name==node])
				while(i<pos):
					#Speciation case:
					s_node = specietree.search_nodes(name=node)[0]
					right_child= s_node.get_child_at(0)
					left_child= s_node.get_child_at(1)

					# Case where the child we're looking for is not even in the geneTree
					r_child_not_there, l_child_not_there= (False, False)

					if(right_child.name not in gene_root_species.get_descendant_name()):
						r_child_not_there=True

					elif(left_child.name not in gene_root_species.get_descendant_name()):
							l_child_not_there=True	

					s_node.add_features(species=node)
					child_name = s_node.get_children_name()
					#Now we should find the index in the matrix of the best node to join
					ind_to_keep=findSpeciationBestJoint(gene_matrix, node_order, s_node, node_in_tree, maxVal)
					# the two node to join are found
					if(ind_to_keep and len(ind_to_keep)==2):
						#node_structs= [ node_s for node_s in node_in_tree if (node_s.name in map(lambda x:node_order[x],ind_to_keep ) and ((node_s.is_leaf() and node_s.species == genetree.search_nodes(name=node_s.name)[0].species) or not(node_s.is_leaf())))]
						# carefully choose the true node, not the added(in case of lost)
						node_structs= [ node_s for node_s in node_in_tree if (node_s.name in map(lambda x:node_order[x],ind_to_keep ) and (not node_s.has_feature('lostnode') or not node_s.lostnode==1))]
						upgma=Upgma.UPGMA_cluster(getMatrix(node_structs, gene_matrix, node_order, ind_to_keep, got_ind=True), node_structs, maxVal)
						spec_tree= upgma[0]
						spec_tree.name="-".join([spec_tree.get_child_at(0).name, spec_tree.get_child_at(1).name])
						spec_tree.add_features(type="SPEC")
						merged_index=map(lambda x: ind_to_keep[x], upgma[2])
						spec_tree.add_features(species=node)
						spec_tree.add_features(score=gene_matrix[merged_index[0], merged_index[1]])
						gene_matrix=Upgma.del_row_column(gene_matrix,merged_index, maxVal)
						node_order[merged_index[0]]=spec_tree.name
						node_order.pop(merged_index[1])
						#Remove child from path_table and add parent
						possible_path_name= [x.name for x in node_in_tree]
						for n in spec_tree.get_children_name():
							while(n in possible_path_name):
								node_in_tree.pop(possible_path_name.index(n))
								possible_path_name.pop(possible_path_name.index(n))
						node_in_tree.append(spec_tree)

						#too much repetition, should make a function to do this part and the duplication part

					# we certainly have a case of a lost node here
					else:
						#avoir la liste des enfants probable et la liste des enfants perdus dans la table
						r_child_list = [x for x in node_in_tree if x.species==right_child.name]
						r_child_lst_list = [x for x in lost_nodes if x.species==right_child.name]
						
						l_child_list = [x for x in node_in_tree if x.species==left_child.name]
						l_child_lst_list = [x for x in lost_nodes if x.species==left_child.name]
						

						# Control 2, avoir la liste des enfants non-presents dans geneTree
						# check which node is lost and make the correct choice
						if(not r_child_list and r_child_lst_list):
							"Cas d'un noeud avec une perte droite"
							for n in l_child_list:
								n_copy= n.copy()
								n_copy.add_features(species=s_node.species)
								n_copy.add_features(lostnode=1)
								already_in= [nd for nd in node_in_tree if(nd.name==n_copy.name and nd.species==n_copy.species and nd.has_feature('lostnode'))]
								if(not already_in):
									node_in_tree.append(n_copy)

						elif(not l_child_list and l_child_lst_list):
							"Cas d'un noeud avec une perte gauche"
							for n in r_child_list:
								n_copy= n.copy()
								n_copy.add_features(species=s_node.species)
								n_copy.add_features(lostnode=1)
								already_in= [nd for nd in node_in_tree if(nd.name==n_copy.name and nd.species==n_copy.species and nd.has_feature('lostnode'))]
								if(not already_in):
									node_in_tree.append(n_copy)


						if(r_child_not_there or l_child_not_there):
							"Cas de l'absence de l'espece dans le geneTree"
							try:
								n= right_child if r_child_not_there else left_child
								match_node=[nd for nd in node_in_tree if (nd.species==s_node.name and nd.has_feature('lostnode') and nd.lostnode==1)]
								for nd in match_node:
									node_in_tree.remove(nd)
							except:
								pass
					i += 1

			# we have a duplication here, construct the duplicated node
			if(n_node==node and pos>n_pos):
				node_structs= [x for x in node_in_tree if x.species== node]
				node_struct_name = [ x.name for x in node_structs]
				node_structs.extend([x for x in node_in_tree if(x.name in node_struct_name and not x in node_structs)])
				node_structs= [ node_s for node_s in node_structs if (not node_s.has_feature('lostnode') or not node_s.lostnode==1)]
				ind_to_keep=[]
				upgma=Upgma.UPGMA_cluster(getMatrix(node_structs, gene_matrix, node_order, ind_to_keep), node_structs, maxVal, 1)
				dup_tree= upgma[0]
				dup_tree.name="-".join([dup_tree.get_child_at(0).name, dup_tree.get_child_at(1).name])
				dup_tree.add_features(type="DUP")

				merged_index=map(lambda x: ind_to_keep[x], upgma[2])
				dup_tree.add_features(species=node)
				dup_tree.add_features(score=gene_matrix[merged_index[0], merged_index[1]])
				gene_matrix=Upgma.del_row_column(gene_matrix,merged_index, maxVal)
				node_order[merged_index[0]]=dup_tree.name
				node_order.pop(merged_index[1])

				possible_path_name= [x.name for x in node_in_tree]
				for n in dup_tree.get_children_name():
					while(n in possible_path_name):
						node_in_tree.pop(possible_path_name.index(n))
						possible_path_name.pop(possible_path_name.index(n))
				node_in_tree.append(dup_tree)


	# the tree is constructed, we should only have one tree in the node_in_tree list
	if(node_in_tree and len(node_in_tree)==1):
		if(verbose):
			print "Total number of node = ", len(node_in_tree[0]), " root specie = ", node_in_tree[0].species, " and score = ",node_in_tree[0].score, "\n\n"
		return node_in_tree[0]
	else:
		print "The tree construction from you path went wrong!, node still not used : ",len(node_in_tree), "\n"
		# for node in node_in_tree:
		# 	print node.get_ascii(show_internal=True, attributes=['species' ])
		# 	print len(node)
		return None

def findSpeciationBestJoint(matrice, node_order, parent_node, node_in_tree, maxVal):
	""" findSpeciationBestJoint find the best node to joint in case of speciation"""
	
	child_0 = parent_node.get_child_at(0).name # find the left child specie
	child_1 = parent_node.get_child_at(1).name
	
	# list of node with the same specie as child_0/child_1
	child_1_list= []
	child_0_list= []

	for x in node_in_tree:
		if x.species==child_1:
			if x.name in node_order:
				child_1_list.append(node_order.index(x.name))
			else:
				child_1_list.append(-1)

		elif x.species==child_0:

			if(x.name in node_order):
				child_0_list.append(node_order.index(x.name))
			else: 
				child_0_list.append(-1)

	# find the best node to join (minimal cost)
	min_val = maxVal
	join_index=[]
	for x_0 in child_0_list:
		for x_1 in child_1_list:
			if(x_0>=0 and x_1>=0 and matrice[x_0, x_1]<min_val):
					min_val=matrice[x_0, x_1]
					join_index=[x_0, x_1]
	
	return join_index



def getMatrix(node_struct, gene_matrix, node_order, ind_to_keep, got_ind=False):
	"""Extract the correct position of the gene_matrix for the next upgma join"""

	if(not got_ind):
		for node in node_struct:
			ind_to_keep.append(node_order.index(node.name))
	matrix= gene_matrix #same reference here
	del_array=[]
	for x in range(len(node_order)):
		if(x not in ind_to_keep):
			del_array.append(x)
	matrix=numpy.delete(matrix, del_array,0)
	matrix=numpy.delete(matrix,del_array,1)
	return matrix


def getReversedMap(genetree, specietree):
	"""Find the reversed map (map between specietree node and genetree node)"""
	mapGene={}
	for node in specietree.traverse():
		mapGene[node]=genetree.search_nodes(species=node.name)

	return mapGene


def getSpecieCount(tree):
	"""Species distribution in the genetree"""
	count= collections.defaultdict(int)
	for node in tree.get_children():
		count[node.species]+=1
	return count


def polytomy_preprocessing(polytomy, specietree, gene_matrix, node_order, maxVal):
	
	for node in polytomy.traverse("postorder"):
		if(node.is_binary()):
			children_list= node.get_children()
			ind_order=[]
			if (not node.has_feature('species')) or node.species=="Unknown":
				species = set([i.species for i in children_list])
				if(len(species)>1):
					node.add_features(species=specietree.get_common_ancestor(species).name)
				else:
					node.add_features(species=species.pop())

			if(	node.name==TreeClass.DEFAULT_NAME):
				global PARTIAL_RESOLUTION_ITERATOR
				node.name="%s_%i"%(node.species, PARTIAL_RESOLUTION_ITERATOR)
				PARTIAL_RESOLUTION_ITERATOR+=1

			upgma=Upgma.UPGMA_cluster(getMatrix(children_list, gene_matrix, node_order, ind_order, got_ind=False), children_list, maxVal)
			merged_index=map(lambda x: ind_order[x], upgma[2])
			node.add_features(score=gene_matrix[merged_index[0], merged_index[1]])
			gene_matrix=Upgma.del_row_column(gene_matrix,merged_index, maxVal)
			node_order[merged_index[0]]=node.name
			node_order.pop(merged_index[1])
	return gene_matrix, node_order


def findMaxX(polytomy, specietree):
	
	polytomy_name_set=set(polytomy.get_children_species())
	if(len(polytomy_name_set)==1):
		polytomy_specie_ancestor=specietree.search_nodes(name=list(polytomy_name_set)[0])[0]
	else:
		polytomy_specie_ancestor= specietree.get_common_ancestor(polytomy_name_set)

	for leaf in polytomy_specie_ancestor.traverse("postorder"):
		parent= leaf.up
		if( not leaf.is_leaf()):
			if(len(set(leaf.get_descendant_name()).intersection(polytomy_name_set))==0):
				parent.remove_child(leaf)

			else:
				polytomy_name_set.add(leaf.name)

		else:
			if(parent is not None) and (len(set(parent.get_children_name()).intersection(polytomy_name_set))==0):
				parent.remove_child(leaf)
			else:
				polytomy_name_set.add(leaf.name)
	

	row_node_corr={}
	n_row=len(polytomy_name_set)-1

	for node in specietree.traverse("levelorder"):
		if(node.name in polytomy_name_set):
			row_node_corr[n_row]=node
			n_row-=1

	return polytomy_name_set, row_node_corr


def solvePolytomy(genetree_nw, specietree_nw, distance_file=None, verbose=0, poly_limit=-1):

	maxVal=1e305
	"""specietree=TreeClass("(((a,b)e,g)h,(c,d)f)i;", format=1)
	genetree= TreeClass("((a_1, b_1, a_2)a_b_1, c_2, b_2, c_1, (d_1,d_2)d_1-d_2, (c_3, d_3)f_1);", format=1)
	s_map={'c_1':'c','c_2':'c', 'd_1':'d', 'a_1':'a', 'b_1':'b', 'a_2':'a', 'b_2':'b','d_2':'d', 'd_3':'d', 'c_3':'c'}
	genetree.set_species(speciesMap=s_map)
	gene_matrix = Utils.makeFakeDstMatrice(10,1,4,1e305)
	node_order =[ 'c_1', 'd_1', 'd_2','c_2', 'c_3', 'd_3' ,'a_1', 'b_1','a_2', 'b_2']
	"""
	#specietree input
	speciename, specie_sep=TreeUtils.newick_preprocessing(specietree_nw)
	specietree=TreeClass(speciename)
	TreeUtils.label_internal_node(specietree)
	#genetree input
	genename, gene_sep=TreeUtils.newick_preprocessing(genetree_nw)
	genetree =TreeClass(genename)
	genetree.set_species(sep=gene_sep, capitalize=True, pos="postfix")

	TreeUtils.reset_node_name(genetree, gene_sep)

	#distance matrice input
	if(distance_file):
		gene_matrix, node_order= Utils.distMatProcessor(distance_file, maxVal)
	else:
		node_order= genetree.get_leaf_names()
		gene_matrix= Utils.makeFakeDstMatrice(len(node_order), 0, 1, maxVal)

	#Find list of species not in genetree
	specieGeneList= set(genetree.get_leaf_species())
	specieList= set([x.name for x in specietree.get_leaves()])
	print "Specie not in genetree = ", specieGeneList-specieList
	#"""
	#Start with only one polytomy
	nb_polytomy=0
	polysolution = [genetree]
	while True:
		next_tree_solution=[] #next list of partially resolved polytomies
		for tree in polysolution:
			for polytomy in tree.iter_polytomies(strategy="postorder"):
				nb_polytomy+=1
				#copying the input for each step, necessary in order to not modify by reference
				matrice = numpy.copy(gene_matrix)
				sptree= specietree.copy("newick") 
				ptree= polytomy.copy()
				order= node_order[:]
				poly_parent= polytomy.up
				node_to_replace=polytomy
				matrice, order=polytomy_preprocessing(ptree, sptree, matrice, order, maxVal)
				#print ptree.get_ascii(show_internal=False, attributes=['species', 'name'])
				solution=polySolver(ptree,sptree, matrice, order,poly_limit)

				if(poly_parent== None):
					#print "root"
					next_tree_solution.extend(solution)

				else:
					#print "internal polytomy"
					for sol in solution:
						poly_parent.replace_node(node_to_replace, sol)
						node_to_replace=sol
						next_tree_solution.append(tree.copy())
				# only solve one polytomy per iteration		
				break
		if not next_tree_solution:
			break
		polysolution = next_tree_solution

	if(nb_polytomy<1):
		raise ValueError("Polytomy not found in your gene tree")

	return [t.copy("newick-extended") for t in polysolution]

if __name__ == '__main__':
	
	
	corrected_tree= TreeClass("famille_1.corrected")
	c_tree= TreeClass("famille_1.corrected")
	corrected_tree.set_species(sep='__', capitalize=True, pos="postfix")
	TreeUtils.reset_node_name(corrected_tree, '__')
	TreeUtils.reset_node_name(c_tree, '__')

	taille=len(corrected_tree.get_descendants())
	arbre_species= TreeClass('speciestree.newick')
	mapping=TreeUtils.lcaMapping(corrected_tree, arbre_species)
	#print corrected_tree.get_ascii(attributes=['name', 'species'], show_internal=False)

	TreeUtils.reconcile(corrected_tree, mapping, "yes")
	c_score=TreeUtils.ComputeDupLostScore(corrected_tree)
	#"""
	polysolution = solvePolytomy('nostar_genetree.tree', 'speciestree.newick', 'nostar_famille_1.dist', poly_limit=-1)
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
		Utils.customTreeCompare(TreeClass('nostar_genetree.tree'), c_tree, tree)
		print '**tree has', len(tree.get_descendants())+1, "nodes and corrected tree has", taille+1 , 'nodes'
		#tree.set_species(sep='_', capitalize=False, pos="prefix")
		mapping=TreeUtils.lcaMapping(tree, arbre_species)
		TreeUtils.reconcile(tree, mapping, "yes")
		score=TreeUtils.ComputeDupLostScore(tree)
		#TreeUtils.CleanFeatures(tree, ['type'])
		print "**m_score= ", m_score_sum, 'r_score= ', score, "c_score= ",c_score, 'reconcile tree node= ', len(tree.get_descendants())+1, " leaves: ", len(tree), " reconcile corrected tree node= ", len(corrected_tree.get_descendants())+1, " ", len(corrected_tree)
		print "**rf=", rf, "max_rf=", max_rf
		print "**tree still has polytomies= ", tree.has_polytomies()


