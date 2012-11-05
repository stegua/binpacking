## Stefano Gualandi, stefano.gualandi@gmail.com
## August, 2012

# Change this file to use a different set of local variable (machine dependant_
include config.topsy

gepack: ${LIB}/dag_pack.o ${SRC}/gepack.cc
	${COMPILER} -c ${SRC}/gepack.cc -o ${LIB}/gepack.o -I${GECODE_INCLUDE} -I${BOOST_INCLUDE} -I${INCLUDE}
	${LINKER} -o ${BIN}/gepack ${LIB}/gepack.o ${LIB}/dag_pack.o ${GECODE_LIB}

binpacking: ${LIB}/dag_pack.o ${SRC}/binpacking.cc
	${COMPILER} -c ${SRC}/binpacking.cc -o ${LIB}/binpacking.o -I${BOOST_INCLUDE} -I${INCLUDE}
	${LINKER} -o ${BIN}/binpacking ${LIB}/binpacking.o ${LIB}/dag_pack.o

# Build everything
all: binpacking

#--------------
# MY LIBRARIES
#--------------
${LIB}/dag_pack.o: ${SRC}/dag_pack.cc
	${COMPILER} -c ${SRC}/dag_pack.cc -o ${LIB}/dag_pack.o -I${BOOST_INCLUDE} -I${INCLUDE}

clean::
	rm -f *.o
	rm -f ${LIB}/*.o
	rm -f *~
	rm -f ${SRC}/*~ ${INCLUDE}/*~
