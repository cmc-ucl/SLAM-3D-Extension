#CC=/usr/local/Cellar/gcc/12.1.0/bin/g++-12 -std=c++20 -I /usr/local/Cellar/eigen/3.4.0_1/include/eigen3
#ICC=icc -std=c++20 -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -fopenmp -I/usr/local/Cellar/libomp/14.0.6/include -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3/Eigen
#ICC=icc -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -fopenmp -I/usr/local/Cellar/libomp/14.0.6/include -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3/Eigen
CC=icpc -std=c++17 -I /home/uccawkj/eigen/eigen-3.4.0/


#slam_ext.x : ext_main.o Cell.o Manager.o Atom.hpp LonePairMatrix.o Integral_lib.o
out.x : ext_main.o Cell.o Manager.o Atom.hpp LonePairMatrix.o Integral_lib.o
	$(CC) -O3 -o out.x ext_main.o Cell.o Manager.o LonePairMatrix.o Integral_lib.o -qopenmp

Cell.o : Cell.cpp Cell.hpp Atom.hpp Manager.hpp
	$(CC) -O3 -c Cell.cpp

Manager.o : Manager.cpp Manager.hpp
	$(CC) -O3 -c Manager.cpp

ext_main.o : ext_main.cpp
	$(CC) -O3 -c ext_main.cpp

LonePairMatrix.o : LonePairMatrix.cpp LonePairMatrix.hpp Integral_lib.hpp
	$(CC) -O3 -c LonePairMatrix.cpp -qopenmp

Integral_lib.o : Integral_lib.cpp Integral_lib.hpp
	$(CC) -O3 -c Integral_lib.cpp

clean : 
	rm -rf *.o

