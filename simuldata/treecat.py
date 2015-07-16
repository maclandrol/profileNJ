#!/usr/bin/env python
import sys

tree_list = []
i = 0;
for tree in sys.argv[1:]:
	with open(tree, 'r') as INPUT: 
		header = ""
		for line in INPUT:
			line = line.strip()
			if line.startswith('>'):
				try:
					header = line.split(';')[1]
				except:
					header = ""
					
			else:
				tree_list.append(">" + tree + " %d "%i + header)
				i += 1
				tree_list.append(line)

for l in tree_list:
	print l