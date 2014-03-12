"""
Author: Emmanuel Noutahi
Date: 02/2014
TreeUtils is a python class that offer function related to phylogeny tree, using TreeClass
"""
from TreeClass import TreeClass
from ete2 import Phyloxml
import httplib2
import os

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
			if(not node.is_leaf() and (lcaMap[node]==lcaMap[node.get_child_at(0)] or lcaMap[node]==lcaMap[node.get_child_at(1)])):
				node.type=TreeClass.AD
				if not (set(node.get_child_at(0).get_species()).intersection(set(node.get_child_at(1).get_species()))):
					node.type=TreeClass.NAD

		if(lost.upper()=="YES"):
			for node in geneTree.traverse("postorder"):
				children_list=node.get_children()
				for child_c in children_list:
					if((lcaMap[child_c].up != lcaMap[node] and lcaMap[child_c] != lcaMap[node]) or (node.type==TreeClass.AD and lcaMap[node]!=lcaMap[child_c])):

						while((lcaMap[child_c].up!=lcaMap[node] and node.type==TreeClass.SPEC) or (lcaMap[child_c]!=lcaMap[node] and node.type!=TreeClass.SPEC)):
							lost=TreeClass()
							intern_lost=TreeClass()
							intern_lost.type=TreeClass.SPEC
							if lcaMap[child_c].is_root():
								intern_lost.species=",".join(lcaMap[child_c].get_leaf_names())
								lcaMap.update({intern_lost:lcaMap[child_c]})

							else:
								intern_lost.species=",".join(lcaMap[child_c].up.get_leaf_names())
								lcaMap.update({intern_lost:lcaMap[child_c].up})


							#change here to display a subtree and not a leaf with a lot of specie
							lost.species=",".join(set(lcaMap[intern_lost].get_leaf_names())-set(child_c.species.split(",")))
							lost.type=TreeClass.LOST
							temp_node=child_c.detach()
							intern_lost.add_child(child=lost)
							intern_lost.add_child(child=temp_node)
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
						lost=TreeClass()
						lost.type=TreeClass.LOST
						lost.species=",".join(unadded_specie)
						node.add_child(lost)

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