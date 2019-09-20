# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= c++ -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG= testcode
${PROG} :	   	tests-main.o  test-functions.o jacobi.o main.o beam.o qho_int.o qho_no_int.o
			${CPPflags} tests-main.o test-functions.o jacobi.o main.o beam.o qho_int.o qho_no_int.o ${LIB} -o ${PROG}

tests-main.o :			tests-main.cpp 
		        	${CPPflags} -c tests-main.cpp

test-functitons.o :			test-functions.cpp
		        	${CPPflags} -c test-functions.cpp

jacobi.o :		jacobi.cpp 
		        	${CPPflags} -c jacobi.cpp

main.o:		main.cpp
					${CPPflags} -c main.cpp

beam.o:		beam.cpp
					${CPPflags} -c beam.cpp

qho_int.o:		qho_int.cpp
					${CPPflags} -c qho_int.cpp

qho_no_int.o:		qho_no_int.cpp
					${CPPflags} -c qho_no_int.cpp

