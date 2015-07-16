"""
Author: Emmanuel Noutahi
Date: 02/2014
params set duplication and losses cost for profileNJ and reconciliation
"""

dupcost, losscost = 1, 1

def set(dup, loss):
	global dupcost
	global losscost
	dupcost, losscost = dup, loss

def getdup(specie):
	global dupcost
	if isinstance(dupcost, dict):
		return dupcost.get(specie, 1)
	else:
		return dupcost

def getloss(species):
	global losscost
	if isinstance(losscost, dict):
		return losscost.get(specie, 1)
	else:
		return losscost

