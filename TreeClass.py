"""
Author: Emmanuel Noutahi
Date: 02/2014
TreeClass is a python class derived from TreeNode class from the ete2 package.
TreeClass add additional specific function
"""
from ete2 import TreeNode
from ete2.phylo import spoverlap
import types
import collections

class TreeClass(TreeNode):

	DEFAULT_SPECIE="Unknown"
	DEFAULT_GENE="Unknown"
	AD=1
	LOST=-1
	SPEC=0
	NAD=2

	def __init__(self, newick=None, format=0, dist=None, support=None,name=None):
		"""	Default init for the TreeClass. This works better than wrapping the entire class"""
		TreeNode.__init__(self, newick=newick, format=format, dist=dist, support=support, name=name)


	def get_child_at(self, i=0):
		"""Return child at a specific position in the children list of a node"""
		children_list=self.get_children()
		if(i<len(children_list)):
			return children_list[i]
		return None


	def has_ancestor(self, ancestor):
		"""Check if "ancestor" is an ancestor of the current TreeNode"""
		ancestor = self.translate_nodes(ancestor)
		ancestors=self.get_ancestors()
		return True if len(list(filter(lambda x : x in ancestors, ancestor)))==len(ancestor) else False


	def has_descendant(self, descendant):
		"""Check if "descendant" is a descendant of the current TreeNode"""
		descendant=self.translate_nodes(descendant)
		descendants=self.get_descendants()
		return True if len(list(filter(lambda x : x in descendants, descendant)))==len(descendant) else False


	def translate_nodes(self, *target_nodes):
		"""Translate list of node name into Node"""

		if len(target_nodes) == 1 and type(target_nodes[0]) in set([set, tuple, list, frozenset]):
			target_nodes = target_nodes[0]

		try:
			target_nodes = [n if isinstance(n, TreeNode) else self.search_nodes(name=n)[0] for n in target_nodes]
			return target_nodes

		except (ValueError, IndexError) as e:
			print "You may have name which cannot be associated with a node"
			raise e


	def insert_child_at(self, i, newNode, replace_if_exist=True):
		"""Insertion of a node in a specific position"""
		removing_child=None
		try:
			if i<len(self.get_children()):
				removing_child=self.children[i]
				if(replace_if_exist):
					self.children[i]=newNode
			else:
				self.add_child(newNode)
		except ValueError, e:
			raise e
		else:
			return removing_child


	def remove_childAt(self, i):
		"""Remove a child at a specific position"""

		try:
			child=self.children[i]
			self.children.remove(child)
		except (IndexError, ValueError), e:
			raise e
		else:
			child.up = None
			return child


	def delete_leaf(self, leafList):
		"""Delete a list of leaf"""
		if len(leafList) == 1 and type(leafList[0]) in set([set, tuple, list, frozenset]):
			leafList = leafList[0]

		for leaf in  leafList:
			if(leaf.is_leaf()):
				leaf.delete()


	def get_degree(self):
		"""Return the degree of the current Node"""
		child_number=len(self.get_children())
		return child_number +1 if self.is_leaf() else child_number


	def set_species(self, speciesMap=None, sep="_", pos="postfix", use_fn=None):

		"""Set species feature for each leaf in the tree.

		:argument speciesMap: Default=None. speciesMap is a Map of species for the geneTree. Each key is a leaf name from the genetree and the value is the corresponding specie name
		:argument sep: Default ="_" , the separator for the default species extraction using the leaf name
		:argument pos: Default="postfix", the species position in the leaf name for the default extraction. Should be used with sep. Can take for value, "prefix", which
                means "specie-sep-gene" or "postfix" for "gene-sep-specie"
		argument fn: Pointer to a parsing python function that receives a node as first argument and returns the species name.

		"""

		for leaf in self:
			if speciesMap is not None:
				leaf.add_features(species=speciesMap.get(leaf.name, TreeClass.DEFAULT_SPECIE))
			elif use_fn is not None :
				leaf.add_features(species=use_fn(leaf))
			else:
				leaf.add_features(species=leaf._extractFeatureName(separator=sep, order=pos))


	def set_genes(self, genesMap=None, sep="_", pos="postfix", use_fn=None):

		"""Set gene feature for each leaf in the tree.

		:argument genesMap: Default=None. genesMap is a Map of genes for the geneTree. Each key is a leaf name from the genetree and the value is the corresponding genes name
		:argument sep: Default ="_" , the separator for the default genes extraction using the leaf name
		:argument pos: Default="postfix", the gene position in the leaf name for the default extraction. Should be used with sep. Can take for value, "postfix", which
                means "specie-sep-gene" or "prefix" for "gene-sep-specie"
		argument fn: Pointer to a parsing python function that receives a node as first argument and returns the genes name.

		"""

		for leaf in self:
			if genesMap is not None:
				leaf.add_features(genes=genesMap.get(leaf.name, TreeClass.DEFAULT_GENE))
			elif use_fn is not None :
				leaf.add_features(genes=use_fn(leaf))
			else:
				leaf.add_features(genes=leaf._extractFeatureName(separator=sep, order=pos))


	def get_species(self, sep=","):
		"""Return the a list of species for the current node"""
		return self.species.split(sep)

	def get_genes(self, sep=","):
		"""Return the a list of genes for the current node"""
		return self.genes.split(sep)

	def _extractFeatureName(self, separator=None, order=None):
		"""Private function, extract feature name (e.g. genes, species) based on the node name"""
		l=self.name.split(separator)
		if len(l)>1 and order=="postfix":
			feature=l[-1]
		elif len(l)>1 and order=="prefix":
			feature=l[0]
		else:
			feature=self.name
		return feature


	def treeContraction(self, seuil=0):
		"""Contract tree based on the dist between node, using a threshold. Any branches shorter than "seuil" will be removed """
		for node in self.traverse():
			if(node.support<seuil):
				node.delete()


	def restrictToSpecies(self, species=[]):
		"""Restrict the current genetree to the list of species passed in argument"""
		hits=[]
		try:
			for value in species:
				hits.extend(self.get_leaves_by_feature(species=value))
			self.prune(hits)
		except Exception as e:
			print e
			print "Check if this tree have species as feature"


	def toPolytomy(self):
		"""Move every leaves to the node by deleting all the internals nodes"""
		for node in self.traverse():
			if(not node.is_leaf() and not node.is_root()):
				node.delete()


	def get_leaves_by_feature(self, **condition):
		"""Return leaves that match the features passed as argument"""
		condition['children']=[]
		return self.search_nodes(**condition)


	def get_all_features(self):
		"""Return all the features of a node in a list object"""
		features_list=[]
		for node in self.traverse():
			features_list.extend(list(node.features))
		return set(features_list)

	def __repr__(self):
		return "Tree Class '%s' (%s)" %(self.name, hex(self.__hash__()))


	def get_my_evol_events(self, sos_thr=0.0):
		""" Returns a list of duplication and speciation events in
		which the current node has been involved. Scanned nodes are
		also labeled internally as dup=True|False. You can access this
		labels using the 'node.dup' sintaxis.

		Method: the algorithm scans all nodes from the given leafName to
		the root. Nodes are assumed to be duplications when a species
		overlap is found between its child linages. Method is described
		more detail in:

		"The Human Phylome." Huerta-Cepas J, Dopazo H, Dopazo J, Gabaldon
		T. Genome Biol. 2007;8(6):R109.
		"""
		return spoverlap.get_evol_events_from_leaf(self, sos_thr=sos_thr)

	def get_descendant_evol_events(self, sos_thr=0.0):
		""" Returns a list of **all** duplication and speciation
		events detected after this node. Nodes are assumed to be
		duplications when a species overlap is found between its child
		linages. Method is described more detail in:

		"The Human Phylome." Huerta-Cepas J, Dopazo H, Dopazo J, Gabaldon
		T. Genome Biol. 2007;8(6):R109.
		"""
		return spoverlap.get_evol_events_from_root(self, sos_thr=sos_thr)

	def is_monophyletic(specieSet):
		""" Returns True id species names under this node are all
		included in a given list or set of species names."""
		if type(specieSet) != set:
			specieSet = set(specieSet)
		return set(self.get_species()).issubset(species)


	#"""Should I implement this???"""
	def delete_single_child_descendant(self,*args, **kargs): pass


	#"""Too Lazy to do that, let just ignore it for the moment"""
	def insertNodeBetween(): pass

	@classmethod
	def import_from_PhyloxmlTree(cls,phyloxml):
		"""import Tree structure and useful Tree features from a _Phyloxml.PhyloxmlTree to a TreeClass
		This is really dirty but it does the job!!.
		****For each attribut, only the first xml value in the clade is transferred to the TreeClass instance***
		accessible feature: Most of the features are list() or dict(). Because phyloxml format support multi Tree and we can have muliple infos per node!
		-code : the taxon code at the node
		-sc_name: scientific_name for the current node, if it was retrieved from emsembl
		-c_name:common_name for the current node
		-taxon_id: taxon_id from the node
		-'type' : types of the sequence (cdna/protein...
		-'symbol': symbol of the gene
		-'name': name of the gene
		-'accession': accession number of the gene
		-'mol_seq': aa sequence or nuc sequence
		"""
		if(phyloxml.__class__.__name__!="PhyloxmlTree"):
			raise ValueError("Please provide a phyloxml class")

		for node in phyloxml:
			clade = node.phyloxml_clade
			sequence=collections.defaultdict(list)
			taxa=collections.defaultdict(list)
			for seq in clade.get_sequence():
				sequence['name'].append(seq.get_name())
				sequence['symbol'].append(seq.get_symbol())
				if(seq.get_mol_seq() is not None):
					sequence['mol_seq'].append(seq.get_mol_seq().get_valueOf_())
				sequence['type'].append(seq.get_type())
				if(seq.get_accession() is not None):
					sequence['accession'].append(seq.get_accession().get_valueOf_())

			for taxon in clade.get_taxonomy():
				taxa['common_name'].append(taxon.common_name)
				taxa['sc_name'].append(taxon.scientific_name)
				taxa['code'].append(taxon.code)
				if(taxon.id is not None):
					taxa['taxon_id'].append(taxon.id.get_valueOf_())

			if(len(taxa['code'])>=1):
				node.add_features(code=taxa['code'][0])
			if(len(taxa['sc_name'])>=1):
				node.add_features(sc_name=taxa['sc_name'][0])
			if(len(taxa['common_name'])>=1):
				node.add_features(c_name=taxa['common_name'][0])
			if(len(taxa['taxon_id'])>=1):
				node.add_features(tax_id=taxa['taxon_id'][0])

			if(len(sequence['accession'])>=1):
				node.add_features(accession=sequence['accession'][0])
			if(len(sequence['mol_seq']) >=1):
				node.add_features(sequence=sequence['mol_seq'][0])

			if(len(sequence['name'])>=1):
				node.add_features(seqname=sequence['name'][0])
			if(len(sequence['symbol'])>=1):
				node.add_features(symbol=sequence['symbol'][0])

		return TreeClass(phyloxml.write(features=[],format_root_node=True))


	def writeSeqToFasta(self, out='seq.fasta', comment=1):
		if("sequence" in self.get_all_features()):
			with open(out, 'w') as outfile:
				for leaf in self:
					if("sequence" in leaf.features):
						id=">%s" % leaf.name
						if(comment and "accession" in leaf.features and "seqname" in leaf.features):
							id= id+ " %s;%s" %(leaf.accession, leaf.seqname)
							if("sc_name" in leaf.features):
								id= id+ (";%s"%leaf.sc_name)
						seq=leaf.sequence+"\n"
						id=id+"\n"
						outfile.write(id)
						outfile.write(seq)
