CC=/usr/local/Cellar/gcc/10.2.0/bin/g++-10 -std=c++20 -I /usr/local/Cellar/eigen/3.4.0_1/include/eigen3

run : test.o Cell.o Manager.o Atom.hpp
	$(CC) -o run test.o Cell.o Manager.o

Cell.o : Cell.cpp Cell.hpp Atom.hpp Manager.hpp
	$(CC) -c Cell.cpp

Manager.o : Manager.cpp Manager.hpp
	$(CC) -c Manager.cpp

test.o : test.cpp
	$(CC) -c test.cpp


clean : 
	rm -rf *.o
