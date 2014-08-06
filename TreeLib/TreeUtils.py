"""
Author: Emmanuel Noutahi
Date: 02/2014
TreeUtils is a python class that offer function related to phylogeny tree, using TreeClass
"""
from TreeClass import TreeClass
import ClusterUtils as clu
from ete2 import Phyloxml
from ete2.parser.newick import NewickError
import httplib2
import hashlib, re
import os
import collections, string


#TreeUtils:

def fetch_ensembl_genetree_by_id(treeID=None,aligned=0, sequence="none", output="nh", nh_format="full"):
	"""Fetch genetree from ensembl tree ID

	:argument treeID: the ensembl tree ID, this is mandatory
	:argument aligned: boolean (0/1), used with sequence to retrieve aligned sequence
	:argument sequence: none / protein /cdna /gene, should we retrieve sequence also?, work only with phyloxml nh_format
	:argument output: nh / phyloxml, type of output we are looking for!
	:argument nh_format: full / display_label_composite / simple / species / species_short_name / ncbi_taxon / ncbi_name / njtree / phylip, The format of the nh output, only useful when the output is set to nh

	"""
	if not treeID:
		raise valueError('Please provide a genetree id')
	else:
		http = httplib2.Http(".cache")
		server = "http://beta.rest.ensembl.org"
		ext = "/genetree/id/%s?sequence=%s;aligned=%i" %(treeID,sequence,aligned)
		if(output=="nh"):
			ext = ext+";nh_format=%s" %nh_format
		output="text/x-"+output
		resp, content = http.request(server+ext, method="GET", headers={"Content-Type":output})
		if not resp.status == 200:
			print "Invalid response: ", resp.status
			raise ValueError('Failled to process request!')

		if(output.lower()!="text/x-phyloxml"):
			return TreeClass(content)
		else:
			return getTreeFromPhyloxml(content)




def fetch_ensembl_genetree_by_member(memberID=None, species=None, id_type=None, output="nh", nh_format="full"):

	"""Fetch genetree from a member ID

	:argument memberID: the ensembl gene ID member of the tree to fetch, this is mandatory! EX: ENSG00000157764
	:argument species: Registry name/aliases used to restrict searches by. Only required if a stable ID is not unique to a species (not the case with Ensembl databases) EX: human, homo_sapiens
	:argument id_type: Object type to restrict searches to. Used when a stable ID is not unique to a single class. EX: gene, transcript
	:argument output: nh / phyloxml, type of output we are looking for!
	:argument nh_format: full / display_label_composite / simple / species / species_short_name / ncbi_taxon / ncbi_name / njtree / phylip, The format of the nh output, only useful when the output is set to nh

	"""
	if not memberID:
		raise valueError('Please provide a genetree id')
	else:
		http = httplib2.Http(".cache")
		server = "http://beta.rest.ensembl.org"
		ext = "/genetree/member/id/%s?" %(memberID)
		if species:
			ext=ext+"species="+species+";"
		if id_type:
			ext= ext+"object_type="+id_type+";"
		if(output=="nh"):
			ext = ext+"nh_format=%s;" %nh_format
		output="text/x-"+output
		resp, content = http.request(server+ext, method="GET", headers={"Content-Type":output})
		if not resp.status == 200:
			print "Invalid response: ", resp.status
			raise ValueError('Failled to process request!')
		if(output.lower()!="text/x-phyloxml"):
			return TreeClass(content)
		else:
			return getTreeFromPhyloxml(content)


def lcaMapping(geneTree, specieTree):
	"""Map the geneTree to a specieTree"""
	mapping ={}
	try:
		for node in geneTree.traverse(strategy="postorder"):
			if not node.is_leaf():
				leaf_under_node= node.get_leaves()
				species = set([i.species for i in leaf_under_node])
				if(len(species)>1):
					mapping[node]= specieTree.get_common_ancestor(species)
				else:
					mapping[node]=specieTree.get_leaves_by_name(list(species)[0])[0]
				#node.add_features(species=mapping[node].name.replace("/",",")) ######Change will depends
				node.add_features(species=",".join(species))
			else:
				mapping[node]=specieTree.search_nodes(name=node.species)[0]
	except Exception as e:
		print type(e)
		print("Leaves without species")
	else :
		return mapping


def reconcile(geneTree=None, lcaMap=None, lost="no"):
	"""Reconcile genetree topology to a specieTree, using an adequate mapping obtained with lcaMapping.
	'reconcile' will infer evolutionary events like gene lost, gene speciation and gene duplication with distinction between AD and NAD
	"""
	if(map is None or geneTree is None):
		raise Exception("lcaMapping or geneTree not found")
	else :
		for node in geneTree.traverse("levelorder"):
			node.add_features(type=TreeClass.SPEC)
			#print node.name , node.species, " and children name ", node.get_children_name()," and children species ", node.get_children_species()
			if(not node.is_leaf() and (lcaMap[node]==lcaMap[node.get_child_at(0)] or lcaMap[node]==lcaMap[node.get_child_at(1)])):
				node.type=TreeClass.AD
				#print "\n\nnode = ", node, "\n\nand children : ", node.children
				if not (set(node.get_child_at(0).get_species()).intersection(set(node.get_child_at(1).get_species()))):
					node.type=TreeClass.NAD

		if(lost.upper()=="YES"):
			for node in geneTree.traverse("postorder"):
				children_list=node.get_children()
				for child_c in children_list:
					if((lcaMap[child_c].up != lcaMap[node] and lcaMap[child_c] != lcaMap[node]) or (node.type==TreeClass.AD and lcaMap[node]!=lcaMap[child_c])):

						while((lcaMap[child_c].up!=lcaMap[node] and node.type==TreeClass.SPEC) or (lcaMap[child_c]!=lcaMap[node] and node.type!=TreeClass.SPEC)):
							lostnode=TreeClass()
							intern_lost=TreeClass()
							intern_lost.type=TreeClass.SPEC
							if lcaMap[child_c].is_root():
								intern_lost.species=",".join(lcaMap[child_c].get_leaf_names())
								lcaMap.update({intern_lost:lcaMap[child_c]})

							else:
								intern_lost.species=",".join(lcaMap[child_c].up.get_leaf_names())
								lcaMap.update({intern_lost:lcaMap[child_c].up})


							#change here to display a subtree and not a leaf with a lot of specie
							lostnode.species=",".join(set(lcaMap[intern_lost].get_leaf_names())-set(child_c.species.split(",")))
							lostnode.type=TreeClass.LOST
							child_c.detach()
							#print "***********************\n\n** node : ", node, "\n\n** child_c: ", child_c, "\n\n** child parent", child_c.up
							#node.remove_child(child_c)
							intern_lost.add_child(child=lostnode)
							intern_lost.add_child(child=child_c)
							child_c=intern_lost
						node.add_child(child_c)
						children_list.append(child_c)

				#Case of polytomie in species tree....
				if not node.is_leaf():
					specie_list = ",".join([",".join(lcaMap[child_c].get_leaf_names()) for child_c in node.get_children()])
					child_specie_set=set(specie_list.split(","))
					real_specie_list=set(lcaMap[node].get_leaf_names())
					unadded_specie=real_specie_list-child_specie_set
					#print unadded_specie, child_specie_set, real_specie_list
					#print node.species
					if(unadded_specie):
						lostnode=TreeClass()
						lostnode.type=TreeClass.LOST
						lostnode.species=",".join(unadded_specie)
						node.add_child(lostnode)


def ComputeDupLostScore(genetree=None):
	"""
	Compute the duplication and lost cost
	"""
	if(genetree is None or 'type' not in genetree.get_all_features()):
		raise Exception("Your Genetree didn't undergoes reconciliation yet")
	count=0
	for node in genetree.traverse("levelorder"):
		if(node.has_feature('type') and node.type==TreeClass.NAD or node.type==TreeClass.LOST or node.type==TreeClass.AD):
			count+=1
	return count


def CleanFeatures(tree=None, features=[]):
	cleaned=False
	if(tree):
		for node in tree.traverse():
			for f in features:
				if(node.has_feature(f)):
					node.del_feature(f)
					cleaned=True
	return cleaned


def getTreeFromPhyloxml(xml, saveToFile="default.xml", delFile=True):
	"""
	Read a phylogeny tree from a phyloxml string and return a TreeClass object
	or a list of TreeClass object
	"""
	project = Phyloxml()
	fo=open(saveToFile, "w+")
	fo.write(xml)
	fo.close()
	project.build_from_file(saveToFile)
	treeList=[]
	for tree in project.get_phylogeny():
		treeList.append(TreeClass.import_from_PhyloxmlTree(tree))

	if(delFile):
		os.remove(saveToFile)
	if len(treeList)==1:
		return treeList[0]
	return treeList


def reset_node_name(tree,sep):
	for x in tree.traverse():
		x.name=x.name.split(sep)[0]
	return tree


def make_random_tree(names=list(string.lowercase), contract_seuil=0, feature_to_contract='support'):
	"""Make a random Gene Tree"""
	tree= TreeClass()
	tree.populate(len(names), names_library=names, random_branches=True)
	tree.contract_tree(seuil=contract_seuil, feature=feature_to_contract)
	return tree


def getSpecieCount(tree):
	"""Species distribution in the genetree"""
	count= collections.defaultdict(int)
	for node in tree.get_children():
		count[node.species]+=1
	return count


def getSpecieGeneMap(genetree, specietree):
	"""Find the reversed map (map between specietree node and genetree node)"""
	mapGene={}
	for node in specietree.traverse():
		mapGene[node]=genetree.search_nodes(species=node.name)

	return mapGene


def treeHash(tree):
	"""Hashing the tree based on the sorted node name"""
	newick_str= re.sub("(?<=\()([^()]+?)(?=\))",lambda m: ",".join(sorted(m.group(1).split(","))), tree.write(format=9))
	#print "newick: ", tree.write(format=9), "parsing: ", newick_str
	return hashlib.sha384(newick_str).hexdigest()


def newick_preprocessing(newick, gene_sep=None):
	"""Newick format pre-processing in order to assure its correctness"""
	DEF_SEP_LIST= [';;', ';', '-', '|', '%', '+','/']

	if isinstance(newick, basestring):
		if os.path.exists(newick):
			nw = open(newick, 'rU').read()
		else:
			nw = newick
		nw = nw.strip()
		if nw.endswith(';'):
			nw = nw[:-1]

		if gene_sep is None:
			i=0
			while i< len(DEF_SEP_LIST) and DEF_SEP_LIST[i] not in nw:
				i+=1
			if i<len(DEF_SEP_LIST):
				gene_sep= '%%'
				nw=nw.replace(DEF_SEP_LIST[i], gene_sep)

			elif i>=len(DEF_SEP_LIST) or ';' in nw:
				raise NewickError, \
				'Unable to format your newick file, Bad gene-specie separator or too much special chars'
		nw+=';'
		return nw, gene_sep
	else:
		raise NewickError, \
		"'newick' argument must be either a filename or a newick string."


def polySolverPreprocessing(genetree, specietree, distance_file, capitalize=False, gene_sep = None, specie_pos="postfix", dist_diagonal=1e305, nFlag=False):
	#################################################################
	#TODO :
	#	1) Correct newick
	#	2) Sequence retrieve
	#	3) PhyML to align sequence and make a distance matrice
	#
	#################################################################

	#genetree input
	if isinstance(genetree, basestring):
		genetree, gene_sep=newick_preprocessing(genetree, gene_sep)
		genetree= TreeClass(genetree)
	genetree.set_species(sep=gene_sep, capitalize=capitalize, pos=specie_pos)

	#specietree input
	if isinstance(specietree, basestring):
		specietree, sep=newick_preprocessing(specietree, '')
		specietree= TreeClass(specietree)
		specietree.label_internal_node()

	#distance matrice input
	if(distance_file):
		gene_matrix, node_order= clu.distMatProcessor(distance_file, dist_diagonal, nFlag)
		#Difference check 1
		if set(node_order).difference(set(genetree.get_leaf_names())):
			reset_node_name(genetree, gene_sep)
	else:
		node_order= genetree.get_leaf_names()
		gene_matrix= clu.makeFakeDstMatrice(len(node_order), 0, 1, dist_diagonal) #Alternative, retrieve aligned sequence and run phyML

	#Find list of species not in genetree
	specieGeneList= set(genetree.get_leaf_species())
	specieList= set([x.name for x in specietree.get_leaves()])
	if(specieGeneList-specieList):
		raise Exception("Species in genetree but not in specietree : %s" %(", ".join(specieGeneList-specieList)))

	return genetree, specietree, gene_matrix, node_order


def customTreeCompare(original_t, corrected_t, t):

	#Leaves remaining test and original binary node test
	ct_leaves=[]
	t_leaves=[]
	t_binary=[]
	success=[]
	ct_binary=[]
	for node in original_t.traverse("levelorder"):
		desc_name=set(node.get_leaf_names())
		ct_parent= corrected_t.get_common_ancestor(desc_name)
		t_parent= t.get_common_ancestor(desc_name)
		ctl=set(ct_parent.get_leaf_names())
		tl=set(t_parent.get_leaf_names())
		ct_leaves.append(ctl.difference(desc_name))
		t_leaves.append(tl.difference(desc_name))
		if(node.is_binary() and not node.has_polytomies()):
			ct_binary.append(ct_parent.robinson_foulds(node)[0:3])
			t_binary.append(t_parent.robinson_foulds(node)[0:3])
		success.append(len(tl.difference(ctl))<1)



	ct_success=filter(None, map(lambda x:len(x)<1,ct_leaves))
	t_success=filter(None, map(lambda x:len(x)<1,t_leaves))

	print "\nCorrected Tree binary list rf_fould\n"
	print "\n".join(map(lambda x: "\t".join([str(v) for v in x]), ct_binary))
	print "\nTree binary list rf_fould\n"
	print "\n".join(map(lambda x: "\t".join([str(v) for v in x]), t_binary))

	if(len(ct_success)==len(ct_leaves)):
		print "**Leave remaining success for corrected tree"
		#print "\n".join([str(h) for h in t_success])
	else:
		print "**Corrected tree doesn't follow patern"
		#print "\n".join(map(lambda x: "\t".join([str(v) for v in x]), ct_leaves))


	if(len(t_success)==len(t_leaves)):
		print "**Leave remaining success for tree"
		#print "\n".join([str(h) for h in t_success])
	else:
		print "**Tree doesn't follow patern"
		#print "\n".join(map(lambda x: "\t".join([str(v) for v in x]), t_leaves))

	print "**Compatibility test between tree: ", all(success)
