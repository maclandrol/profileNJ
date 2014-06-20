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
import copy

try:
	import cPickle as pickle
except:
	import pickle

class TreeClass(TreeNode):

	DEFAULT_SPECIE="Unknown"
	DEFAULT_NAME= "NoName"
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
		else:
			raise IndexError("Index out of bound, Can't access the child at index: %i"%i)


	def copy(self, method="cpickle", nw_features=[], nw_format_root_node=True, binary_correct=False):
		""" .. override of ete TreeNode original copy

		Returns a copy of the current node.

		:var cpickle method: Protocol used to copy the node
		structure. The following values are accepted:

			- "newick": Tree topology, node names, branch lengths and
			branch support values will be copied by as represented in
			the newick string (copy by newick string serialisation).

			- "newick-extended": Tree topology and all node features
			will be copied based on the extended newick format
			representation. Only node features will be copied, thus
			excluding other node attributes. As this method is also
			based on newick serialisation, features will be converted
			into text strings when making the copy.

			- "cpickle": The whole node structure and its content is
			cloned based on cPickle object serialisation (slower, but
			recommended for full tree copying)
	
			- "deepcopy": The whole node structure and its content is
			copied based on the standard "copy" Python functionality
			(this is the slowest method but it allows to copy complex
			objects even if attributes point to lambda functions,
			etc.)

		"""
		if method=="newick":
			new_node = self.__class__(self.write(features=["name"], format_root_node=True))
		elif method=="newick-extended":
			new_node = self.__class__(self.write(features=nw_features, format_root_node=nw_format_root_node))
		elif method == "deepcopy":
			parent = self.up
			self.up = None
			new_node = copy.deepcopy(self)
			self.up = parent
		elif method == "cpickle":
			parent = self.up
			self.up = None
			new_node = pickle.loads(pickle.dumps(self, 2))
			self.up = parent

		elif method == 'simplecopy':
			parent = self.up
			self.up = None
			new_node = self._recur_simple_copy(nw_features)
			self.up = parent

		else:
			raise ValuerError("Invalid copy method")

		return self._correct_copy(new_node) if binary_correct else new_node


	def _recur_simple_copy(self, features=[]):
		"""Simple copy of a node by a recursive call"""
		root= self._copy_node(features)
		for node in self.get_children():
			root.add_child(node._recur_simple_copy(features))
		return root


	def _iter_simple_copy(self, features=[]):
		"""Iteratif simple copy, this can be optimized"""
		ori_parents=[self]
		copy_parents=[self._copy_node(features)]

		while ori_parents:
			next_parent=[]
			next_copy_parent=[]
			for i in xrange(len(ori_parents)):
				parent= ori_parents[i]
				root=copy_parents[i]
				for node in parent.get_children():
					copy_node=node._copy_node(features)
					root.add_child(copy_node)
					next_parent.append(node)
					next_copy_parent.append(copy_node)
			ori_parents=next_parent
			copy_parents=next_copy_parent
		return root.get_tree_root()


	def _copy_node(self, features=[]):
		"Copy a node and its features to a new node"
		copy=TreeClass()
		if not features:
			features= self.features
		for feature in features:
			if(self.has_feature(feature)):
				copy.add_feature(feature, getattr(self,feature))
		return copy


	def _correct_copy(self, copy):
		"""Correct the structure of new node copied using newick method"""
		for node in copy.traverse("postorder"):
			if not node.is_root() and (node.name not in list(self.get_descendant_name())):
				node.detach()
			if node.is_internal() and len(node.get_children())<2:
				child=node.get_child_at(0)
				if(child.is_leaf()):
					ori_node=(self&(child.name)).up
					print ori_node
				else:
					ori_node=self.get_common_ancestor(child.get_leaf_name()).up
				node.up.replace_child(node, ori_node.copy("newick-extended"))
		return copy

	def has_ancestor(self, ancestor):
		"""Check if "ancestor" is an ancestor of the current Node"""
		ancestor = self.translate_nodes(ancestor)
		ancestors=self.get_ancestors()
		return True if len(list(filter(lambda x : x in ancestors, ancestor)))==len(ancestor) else False


	def has_descendant(self, descendant):
		"""Check if "descendant" is a descendant of the current Node"""
		descendant=self.translate_nodes(descendant)
		descendants=self.get_descendants()
		return True if len(list(filter(lambda x : x in descendants, descendant)))==len(descendant) else False


	def translate_nodes(self, *target_nodes):
		"""Translate list of node name into Node"""

		if len(target_nodes) == 1 and type(target_nodes[0]) in set([set, tuple, list, frozenset]):
			target_nodes = target_nodes[0]

		try:
			target_nodes = [n if isinstance(n, self.__class__) else self&(n) for n in target_nodes]
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
					newNode.up=self
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


	def set_species(self, speciesMap=None, sep="_", capitalize=False, pos="postfix", use_fn=None):

		"""Set species feature for each leaf in the tree.

		:argument speciesMap: Default=None. speciesMap is a Map of species for the geneTree. Each key is a leaf name from the genetree and the value is the corresponding specie name
		:argument sep: Default ="_" , the separator for the default species extraction using the leaf name
		:argument pos: Default="postfix", the species position in the leaf name for the default extraction. Should be used with sep. Can take for value, "prefix", which
                means "specie-sep-gene" or "postfix" for "gene-sep-specie"
		argument fn: Pointer to a parsing python function that receives a node as first argument and returns the species name.

		"""
		if speciesMap is not None:
			for node in self.traverse():
				node.add_features(species=speciesMap.get(node.name, TreeClass.DEFAULT_SPECIE))
		else:
			for leaf in self:
				if use_fn is not None :
					leaf.add_features(species=use_fn(leaf))
				else:
					leaf.add_features(species=leaf._extractFeatureName(separator=sep, order=pos, cap=capitalize))

	
	def set_genes(self, genesMap=None, sep="_", capitalize=False, pos="postfix", use_fn=None):
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
				leaf.add_features(genes=leaf._extractFeatureName(separator=sep, order=pos,cap=capitalize))


	def get_species(self, sep=","):
		"""Return the list of species for the current node 
		(under the current node after reconciliation)"""
		return self.species.split(sep)


	def get_genes(self, sep=","):
		"""Return the list of genes for the current node
		(under the current node after reconciliation)"""
		return self.genes.split(sep)


	def _extractFeatureName(self, separator=None, order=None, cap=False):
		"""Private function, extract feature name (e.g. genes, species) based on the node name"""
		l=self.name.split(separator)
		if len(l)>1 and order=="postfix":
			feature=l[-1]
		elif len(l)>1 and order=="prefix":
			feature=l[0]
		else:
			feature=self.name
		if cap:
			feature=self.__class__._capitalize(feature)
		return feature


	def contract_tree(self, seuil=0, feature='support', break_tree_topo=False):
		"""Contract tree based on the dist between node, using a threshold. `contract_tree` proceed bottom-up. Any branches with a support less than "seuil" will be removed 
		if `break_tree_topo` is set to True, all the branch under this node will be recursively removed"""
		for node in self.traverse("postorder"):
			if(node.has_feature(feature) and node.is_internal() and getattr(node,feature)<seuil):
				node.toPolytomy(break_tree_topo)
		

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


	def toPolytomy(self, break_tree_topo=False):
		"""Move every leaves to the node by deleting all the internals nodes"""
		if(break_tree_topo):
			for node in self.traverse():
				if(not node.is_leaf() and not node.is_root()):
					node.delete()
		else:
			self.delete()


	def get_leaves_by_feature(self, **condition):
		"""Return leaves that match the features passed as argument"""
		match= self.search_nodes(**condition)
		return [node for node in match if node.is_leaf()]


	def is_polytomy(self):
		"""
        Return True if current node is a polytomy.
    	"""
		return len(self.children)>2


	def is_binary(self):
		"""
        Return True if current node is a binary node.
    	"""
		return len(self.children)==2


	def is_internal(self):
		"""
        Return True if current node is an internal node.
    	"""
		return (not self.is_root() and not self.is_leaf())


	def get_all_features(self):
		"""Return all the features of all nodes under self in a set"""
		features_list=[]
		for node in self.traverse():
			features_list.extend(list(node.features))
		return set(features_list)


	def has_feature(self, feature):
		"""Return weither or not this node has feature in its list of features"""
		return (feature in self.features)


	def __repr__(self):
		return "Tree Class '%s' (%s)" %(self.name, hex(self.__hash__()))


	def reroot(self, root_node=True):
		"""reroot tree at each node"""
		#self.label_internal_node()
		for node in self.iter_descendants():
			c_tree= self.copy("newick-extended", nw_format_root_node=True)
			#c_node =c_tree&node.name
			c_node = c_tree.get_common_ancestor(node.get_leaf_name()) if node.is_internal() else c_tree&node.name
			c_tree.set_outgroup(c_node)
			#case where we root at the node and not at the branch
			if(root_node and not node.is_leaf()):
				root= c_tree.get_tree_root()
				#new_child= [child for child in root.get_children() if child !=node.name][0]
				#rooting_node = [child for child in root.get_children() if child.name ==node.name][0]
				new_child= [child for child in root.get_children() if set(child.get_leaf_name()).symmetric_difference(set(node.get_leaf_name()))][0]
				rooting_node = [child for child in root.get_children() if child!=new_child][0]
				c_tree= rooting_node.detach()
				new_child.detach()
				#new_child.label_internal_node()
				c_tree.add_child(new_child)
				yield c_tree

			elif not root_node:
				yield c_tree



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
	

	def is_monophyletic(self, specieSet):
		""" Returns True if species names under this node are all
		included in a given list or set of species names."""

		if type(specieSet) != set:
			specieSet = set(specieSet)
		return self.get_leaf_species().issubset(specieSet)


	def has_polytomies(self):
		"""Return whether or not this tree has polytomies
		"""
		return len(self.get_polytomies())>0


	def get_children_species(self):
		""" Return the species list of the children under this particular node
		"""
		c_species= set([])
		for node in self.get_children():
			c_species.add(node.species)
		return c_species


	def get_children_name(self):
		""" Return the names of the children under this particular node
		"""
		c_names= set([])
		for node in self.get_children():
			c_names.add(node.name)
		return c_names


	def get_descendant_species(self):
		""" Return the species list of the descendants under this particular node
		"""
		c_species= set([])
		for node in self.get_descendants():
			c_species.add(node.species)
		return c_species


	def get_leaf_species(self, is_leaf_fn=None):
		""" Return the species list of the leave under this node
		"""
		return set([leaf.species for leaf in self.iter_leaves(is_leaf_fn=is_leaf_fn)])


	def get_descendant_name(self):
		""" Return the names of the descendants under this particular node
		"""
		c_names= set([])
		for node in self.get_descendants():
			c_names.add(node.name)
		return c_names


	def get_ancestor_name(self):
		"""Return the names of all the ancestor of this node in a set
		"""
		c_names=set([])
		parent=self.up
		while parent is not None:
			c_names.add(parent.name)
			parent=parent.up
		return c_names


	def get_leaf_name(self, is_leaf_fn=None):
		return self.get_leaf_names(is_leaf_fn)


	def delete_single_child_internal(self):
		for node in self.traverse("postorder"):
			if(node.is_internal() and len(node.children)<2):
				node.delete()


	def iter_polytomies(self, is_polytomy_fn=None, strategy="postorder"):
		"""
		Returns an iterator over the polytomies starting from the curent node
		:argument None is_polytomy_fn: See :func:`TreeNode.traverse` for
		documentation.
		"""
		for n in self.traverse(strategy=strategy):
			if not is_polytomy_fn:
				if n.is_polytomy():
					yield n
			else:
				if is_polytomy_fn(n):
					yield n


	def get_polytomies(self, ind=-1 ,is_polytomy_fn=None):
		"""
		Return a list of polytomies under this node
		"""
		polytomies= [pol for pol in self.iter_polytomies(is_polytomy_fn=is_polytomy_fn)]
		if(ind>-1 and ind<len(polytomies)):
			return polytomies[ind]
		else:
			return polytomies

	
	def label_internal_node(self):
		"""Label the internal node of a specietree for the polysolver algorithm"""
		count=1
		for node in self.traverse(strategy='levelorder'):
			if not node.is_leaf() and node.name in [TreeClass.DEFAULT_NAME, TreeClass.DEFAULT_GENE]:
				node.name="%i"%(count)
				count+=1
		return self


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


	def replace_child(self, child_to_replace, new_child):
		if (self is None) or (child_to_replace not in self.get_children()):
			raise ValueError("Node is None or child_to_replace not valid")
		else :
			self.remove_child(child_to_replace)
			self.add_child(new_child)
			#print "***node swap"
			#print "to replace parent ", child_to_replace.up
			#print "new child parent ", new_child.up
			#print "real parent ", self
			return self


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


	@staticmethod
	def _capitalize(line):
		return "".join([line[0].upper(), line[1:]])
