"""
Author: Emmanuel Noutahi
Date: 02/2014
TreeUtils is a python class that offer function related to phylogeny tree, using TreeClass
"""


dupcost, losscost=1,1
def set(dup, loss):
	global dupcost
	global losscost
	dupcost, losscost= dup, loss
	