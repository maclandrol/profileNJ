from Multipolysolver import *
from TreeLib import *
from TreeLib.SupportUtils import *

# mltree, smap, specietree, alignfile, correctPhylo, seuil, mltree_ext, r_option, slimit, plimit):
smap="exp/fungi.smap"
specietree="exp/fungi.stree"
alignfile="exp/0.align"
mltree="exp/0.align.root.bootstrap.tree"
intree="exp/0.align.bootstrap.tree"
gtree="exp/0.align.bootstrap.nw"

mltree_ext=".align.root.bootstrap.tree"
r_option="none"
slimit=-1
plimit=-1
correctPhylo="exp/0.tree"
seuil=95

maxLTree= TreeClass(intree)
maxLTree.contract_tree(seuil=seuil)
maxLTree.write(outfile=gtree, format=0)

methodCompare(mltree, smap, specietree, alignfile, gtree, seuil, mltree_ext, r_option, slimit, plimit, correctPhylo)

# show help
#treefix_compute --type likelihood -m treefix.models.raxmlmodel.RAxMLModel --show-help
#treefix_compute --type cost -m treefix.models.duplossmodel.DupLossModel --show-help




#=============================
# test dup/loss module

# show help

# compute cost for RAxML tree
#treefix_compute --type cost -m treefix.models.duplossmodel.DupLossModel -r -s config/fungi.stree -S config/fungi.smap -o .nt.raxml.tree  sim-fungi/0/0.nt.raxml.tree

# compute cost for true tree
#treefix_compute --type cost -m treefix.models.duplossmodel.DupLossModel -r -s config/fungi.stree -S config/fungi.smap -o .tree sim-fungi/0/0.tree






