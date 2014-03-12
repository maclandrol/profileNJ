def logging(func):
	"""Wrapper for any method returning a TreeNode This will cast it in TreeClass"""
	def wrapper(*args, **kwargs):
		res = func(*args, **kwargs)
		if type(res) in set([set, tuple, list, frozenset]):
			res= map(lambda x:castToMyNode(x),res)

		elif isinstance(res, TreeNode):
			res.__class__=TreeClass
		return res
	return wrapper

def castToMyNode(x):
	if(isinstance(x, TreeNode)):
		x.__class__=TreeClass
	return x

def hasmethod(obj, name):
	"""Return true is name is a method of obj"""
	return hasattr(obj, name) and type(getattr(obj, name)) == types.MethodType

def loggify(myclass):
	"""Run logging on all the method of myclass except name Mangling method"""
	for x in filter(lambda x:"__" not in x, dir(myclass)):
		if hasmethod(myclass,x):
			# print(x)
			setattr(myclass,x,logging(getattr(myclass,x)))
	return myclass