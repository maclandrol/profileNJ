from Multipolysolver import *
from TreeLib import *
from TreeLib.SupportUtils import *

file = convert_to_phylip('0.align', 'fasta', 'align')
construct_phyml_tree(file)


