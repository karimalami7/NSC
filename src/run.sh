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
#	argv[6-15] = list of wanted methods
# 
#

./main INDE 100 10000 10 4 NSCwM