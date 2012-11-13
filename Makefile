## Stefano Gualandi, stefano.gualandi@gmail.com
## August, 2012

# Change this file to use a different set of local variable (machine dependant)
include config.topsy

gegap: ${LIB}/dag_pack.o ${LIB}/cost_binpacking.o ${SRC}/gegap.cc
	${COMPILER} -c ${SRC}/gegap.cc -o ${LIB}/gegap.o -I${GECODE_INCLUDE} -I${BOOST_INCLUDE} -I${INCLUDE}
	${LINKER} -o ${BIN}/gegap ${LIB}/gegap.o ${LIB}/dag_pack.o ${LIB}/cost_binpacking.o ${LIB}/propagate.o ${GECODE_LIB}

gepack: ${LIB}/dag_pack.o ${SRC}/gepack.cc
	${COMPILER} -c ${SRC}/gepack.cc -o ${LIB}/gepack.o -I${GECODE_INCLUDE} -I${BOOST_INCLUDE} -I${INCLUDE}
	${LINKER} -o ${BIN}/gepack ${LIB}/gepack.o ${LIB}/dag_pack.o ${GECODE_LIB}

binpacking: ${LIB}/dag_pack.o ${SRC}/binpacking.cc
	${COMPILER} -c ${SRC}/binpacking.cc -o ${LIB}/binpacking.o -I${BOOST_INCLUDE} -I${INCLUDE}
	${LINKER} -o ${BIN}/binpacking ${LIB}/binpacking.o ${LIB}/dag_pack.o

# Build everything
all: binpacking

# Propagator
${LIB}/cost_binpacking.o: ${LIB}/propagate.o ${LIB}/dag_pack.o ${SRC}/cost_binpacking.cc
	${COMPILER} -c ${SRC}/cost_binpacking.cc -o ${LIB}/cost_binpacking.o -I${GECODE_INCLUDE} -I${BOOST_INCLUDE} -I${INCLUDE}

${LIB}/propagate.o: ${LIB}/dag_pack.o ${SRC}/propagate.cc
	${COMPILER} -c ${SRC}/propagate.cc -o ${LIB}/propagate.o -I${GECODE_INCLUDE} -I${BOOST_INCLUDE} -I${INCLUDE}

#--------------
# MY LIBRARIES
#--------------
${LIB}/dag_pack.o: ${SRC}/dag_pack.cc
	${COMPILER} -c ${SRC}/dag_pack.cc -o ${LIB}/dag_pack.o -I${BOOST_INCLUDE} -I${INCLUDE} -I${GECODE_INCLUDE}

clean::
	rm -f *.o
	rm -f ${LIB}/*.o
	rm -f *~
	rm -f ${SRC}/*~ ${INCLUDE}/*~
