# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= c++ -O3 -std=c++11
# Here we define the name of the executable
PROG= testcode
${PROG} :	   	tests-main.o  test-functions.o ising.o
			${CPPflags} tests-main.o test-functions.o ising.o -o ${PROG}

tests-main.o :			tests-main.cpp
		        	${CPPflags} -c tests-main.cpp

test-functions.o :		test-functions.cpp
		        	${CPPflags} -c test-functions.cpp

ising.o :			ising.cpp
		        	${CPPflags} -c ising.cpp
