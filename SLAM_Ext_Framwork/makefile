#CC=/usr/local/Cellar/gcc/10.2.0/bin/g++-10 -std=c++20 -I /usr/local/Cellar/eigen/3.4.0_1/include/eigen3
CC=/usr/local/Cellar/gcc/12.2.0/bin/g++-12 -std=c++20 -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3 -fopenmp -I/usr/local/Cellar/libomp/14.0.6/include

run : test.o Cell.o Manager.o Atom.hpp
	$(CC) -O3 -o run test.o Cell.o Manager.o

Cell.o : Cell.cpp Cell.hpp Atom.hpp Manager.hpp
	$(CC) -O3 -c Cell.cpp

Manager.o : Manager.cpp Manager.hpp
	$(CC) -O3 -c Manager.cpp

test.o : test.cpp
	$(CC) -O3 -c test.cpp


clean : 
	rm -rf *.o
