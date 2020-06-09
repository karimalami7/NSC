#/bin/bash
###############################
#
#   argv[1] = ANTI, INDE, CORR, PERS, other
#
#	argv[2] : if  argv[1]==other : dataset path
# 			  if  argv[1]!=other : k = distinct values 
#	
#	argv[3] : n = number of tuples 
#
#	argv[4] : d = space or number of attributes
#
#	argv[5] : NBTHREARDS = number of parallel threads to be run.
#
#	argv[6-15] = list of wanted methods {NSC, NSCwM, TREE, NAIF, CSC}
# 
#

distribution=INDE
dataset_size=100000
distint_values=$(($dataset_size/100))
space=8
nb_threads=8

./exec_NSC $distribution $distint_values $dataset_size $space $nb_threads NSC