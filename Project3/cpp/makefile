# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= c++ -O3
# Here we define the name of the executable
PROG= testcode
${PROG} :	   	tests-main.o  test-functions.o montecarlo.o lib.o
			${CPPflags} tests-main.o test-functions.o montecarlo.o lib.o -o ${PROG}

tests-main.o :			tests-main.cpp
		        	${CPPflags} -c tests-main.cpp

test-functitons.o :			test-functions.cpp
		        	${CPPflags} -c test-functions.cpp

montecarlo.o :		montecarlo.cpp
		        	${CPPflags} -c montecarlo.cpp

lib.o:		lib.cpp
					${CPPflags} -c lib.cpp
