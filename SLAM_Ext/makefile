CC=icpc -std=c++11 -I/opt/include/eigen3/
MCC=mpiicpc -std=c++11 -I/opt/include/eigen3/

out.x : ext_main.o Cell.o Manager.o Atom.hpp LonePairMatrix.o Integral_lib.o
	$(CC) -O3 -o out.x ext_main.o Cell.o Manager.o LonePairMatrix.o Integral_lib.o -qopenmp

Cell.o : Cell.cpp Cell.hpp Atom.hpp Manager.hpp
	$(CC) -O3 -static -c Cell.cpp

Manager.o : Manager.cpp Manager.hpp
	$(CC) -O3 -static -c Manager.cpp

ext_main.o : ext_main.cpp
	$(CC) -O3 -static -c ext_main.cpp

LonePairMatrix.o : LonePairMatrix.cpp LonePairMatrix.hpp Integral_lib.hpp
	$(CC) -O3 -static -c LonePairMatrix.cpp -qopenmp

Integral_lib.o : Integral_lib.cpp Integral_lib.hpp
	$(CC) -O3 -static -c Integral_lib.cpp

clean : 
	rm -rf *.o

