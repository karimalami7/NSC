#/bin/bash
###############################
#
#    argv[1] = ANTI, INDE, CORR, PERS, autre
#
#	 argv[2] : si autre : path (chemin du fichier des donnees)
# 			   si !autre : k 
#	argv[3] : n : number of tuples 
#
#	argv[4] : d : space or number of attributes
#
#	argv[5] : NBTHREARDS : number of parallel threads to be run.
#
#	argv[6-15] : list of wanted methods
#
#

make clean

make

./main INDE 100 10000 14 4 NSCwM