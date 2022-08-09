#CC=/usr/local/Cellar/gcc/12.1.0/bin/g++-12 -std=c++20 -I /usr/local/Cellar/eigen/3.4.0_1/include/eigen3
CC=/usr/local/Cellar/gcc/12.1.0/bin/g++-12 -std=c++20 -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3 -fopenmp -I/usr/local/Cellar/libomp/14.0.6/include

run : test_main.o Cell.o Manager.o Atom.hpp
	$(CC) -o run test_main.o Cell.o Manager.o

Cell.o : Cell.cpp Cell.hpp Atom.hpp Manager.hpp
	$(CC) -c Cell.cpp

Manager.o : Manager.cpp Manager.hpp
	$(CC) -c Manager.cpp

test_main.o : test_main.cpp
	$(CC) -c test_main.cpp


clean : 
	rm -rf *.o
