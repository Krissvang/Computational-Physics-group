# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= c++ -std=c++11
# Here we define the name of the executable
PROG= testcode
${PROG} :	   	tests-main.o  test-functions.o varmontecarlo.o
			${CPPflags} tests-main.o test-functions.o varmontecarlo.o -o ${PROG}

tests-main.o :			tests-main.cpp
		        	${CPPflags} -c tests-main.cpp

test-functions.o :		test-functions.cpp
		        	${CPPflags} -c test-functions.cpp

varmontecarlo.o :			varmontecarlo.cpp
		        	${CPPflags} -c varmontecarlo.cpp
