CFLAGS=-O3 -lm -fopenmp -finline -std=c++11  

NSC: main.o DataGenerator.o Gaussian.o ValueGenerator.o VectorCommon.o lp.o sphere.o operation.o data_utility.o
	g++ $(CFLAGS) main.o DataGenerator.o Gaussian.o ValueGenerator.o VectorCommon.o lp.o sphere.o operation.o data_utility.o -lglpk -o exec_NSC

main.o: 
	g++ $(CFLAGS) -c main.cpp

DataGenerator.o:
	g++ $(CFLAGS) -c generator/DataGenerator.cpp

Gaussian.o:
	g++ $(CFLAGS) -c generator/Gaussian.cpp

ValueGenerator.o:
	g++ $(CFLAGS) -c generator/ValueGenerator.cpp

VectorCommon.o:
	g++ $(CFLAGS) -c generator/VectorCommon.cpp

lp.o:
	g++ $(CFLAGS) -c sphere/lp.cpp

sphere.o:
	g++ $(CFLAGS) -c sphere/sphere.cpp

operation.o:
	g++ $(CFLAGS) -c sphere/operation.cpp

data_utility.o:
	g++ $(CFLAGS) -c sphere/data_utility.cpp


clean:
	rm *.o exec_NSC