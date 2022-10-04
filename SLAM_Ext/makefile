#CC=/usr/local/Cellar/gcc/12.1.0/bin/g++-12 -std=c++20 -I /usr/local/Cellar/eigen/3.4.0_1/include/eigen3
#ICC=icc -std=c++20 -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -fopenmp -I/usr/local/Cellar/libomp/14.0.6/include -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3/Eigen
#ICC=icc -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -fopenmp -I/usr/local/Cellar/libomp/14.0.6/include -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3/Eigen

CC=/usr/local/Cellar/gcc/12.2.0/bin/g++-12 -std=c++20 -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3 -fopenmp -I/usr/local/Cellar/libomp/14.0.6/include

GSL_I=-I/usr/local/include/
GSL_L=-L/usr/local/lib/ -lgsl

slam_ext.x : test_main.o Cell.o Manager.o Atom.hpp LonePairMatrix.o Integral_lib.o
	#$(CC) -O3 -o slam_ext.x test_main.o Cell.o Manager.o LonePairMatrix.o Integral_lib.o $(GSL_L) 	// MacBookAir
	$(CC) -O3 -o slam_ext.x test_main.o Cell.o Manager.o LonePairMatrix.o Integral_lib.o 				# iMac

Cell.o : Cell.cpp Cell.hpp Atom.hpp Manager.hpp
	$(CC) -c Cell.cpp

Manager.o : Manager.cpp Manager.hpp
	$(CC) -c Manager.cpp

test_main.o : test_main.cpp
	$(CC) -c test_main.cpp

LonePairMatrix.o : LonePairMatrix.cpp LonePairMatrix.hpp Integral_lib.hpp
	$(CC) -O3 -c LonePairMatrix.cpp 

Integral_lib.o : Integral_lib.cpp Integral_lib.hpp
	$(CC) -O3 -c Integral_lib.cpp $(GSL_I)

clean : 
	rm -rf *.o

