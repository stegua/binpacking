## Stefano Gualandi, stefano.gualandi@gmail.com
## March, 2013

# Change this file to use a different set of local variable (machine dependant)
include config.mac

example: ${LIB}/multibin-packing.o ${SRC}/example.cc
	${COMPILER} -c ${SRC}/example.cc -o ${LIB}/example.o \
		-I${GECODE_INCLUDE} -I${BOOST_INCLUDE} -I${INCLUDE} -I${CLIQUER_INC}
	${LINKER} -o ${BIN}/example ${LIB}/example.o ${LIB}/multibin-packing.o ${GECODE_LIB} ${CLIQUER_LIB}

gecode_multibin: ${LIB}/propagate_multi.o ${LIB}/dag_pack.o ${LIB}/path.o ${LIB}/cost_multibin.o ${SRC}/gecode_multibin.cc
	${COMPILER} -c ${SRC}/gecode_multibin.cc -o ${LIB}/gecode_multibin.o \
		-I${GECODE_INCLUDE} -I${BOOST_INCLUDE} -I${INCLUDE} -I${GUROBI_INC} -I${CLIQUER_INC}
	${LINKER} -o ${BIN}/gecode_multibin ${LIB}/gecode_multibin.o \
		${LIB}/path.o ${LIB}/dag_pack.o ${LIB}/propagate_multi.o ${LIB}/cost_multibin.o \
		${GUROBI_LIB} ${GECODE_LIB} ${QSOPT}/qsopt.a ${CLIQUER_LIB}

multibinpacking: ${LIB}/path.o ${LIB}/dag_pack.o ${SRC}/multibinpacking.cc
	${COMPILER} -c ${SRC}/multibinpacking.cc -o ${LIB}/multibinpacking.o \
		-I${QSOPT} -I${BOOST_INCLUDE} -I${INCLUDE} -I${GECODE_INCLUDE}
	${LINKER} -o ${BIN}/multibinpacking ${LIB}/multibinpacking.o ${LIB}/dag_pack.o ${LIB}/path.o \
		${GECODE_LIB} ${QSOPT}/qsopt.a

gegap: ${LIB}/path.o ${LIB}/dag_pack.o ${LIB}/cost_binpacking.o ${SRC}/gegap.cc
	${COMPILER} -c ${SRC}/gegap.cc -o ${LIB}/gegap.o -I${GECODE_INCLUDE} -I${BOOST_INCLUDE} -I${INCLUDE}
	${LINKER} -o ${BIN}/gegap ${LIB}/gegap.o ${LIB}/dag_pack.o ${LIB}/cost_binpacking.o ${LIB}/propagate.o \
		${LIB}/path.o ${GECODE_LIB} ${QSOPT}/qsopt.a

gepack: ${LIB}/dag_pack.o ${SRC}/gepack.cc
	${COMPILER} -c ${SRC}/gepack.cc -o ${LIB}/gepack.o -I${GECODE_INCLUDE} -I${BOOST_INCLUDE} -I${INCLUDE}
	${LINKER} -o ${BIN}/gepack ${LIB}/gepack.o ${LIB}/dag_pack.o ${GECODE_LIB} ${QSOPT}/qsopt.a

binpacking: ${LIB}/dag_pack.o ${SRC}/binpacking.cc
	${COMPILER} -c ${SRC}/binpacking.cc -o ${LIB}/binpacking.o -I${BOOST_INCLUDE} -I${INCLUDE}
	${LINKER} -o ${BIN}/binpacking ${LIB}/binpacking.o ${LIB}/dag_pack.o

# Build everything
all: binpacking

#--------------
# MY LIBRARIES
#--------------
${LIB}/multibin-packing.o: ${LIB}/dag_pack.o ${SRC}/multibin-packing.cpp
	${COMPILER} -c ${SRC}/multibin-packing.cpp -o ${LIB}/multibin-packing.o \
		-I${GECODE_INCLUDE} -I${BOOST_INCLUDE} -I${INCLUDE} -I${CLIQUER_INC}

${LIB}/cost_binpacking.o: ${LIB}/propagate.o ${LIB}/dag_pack.o ${SRC}/cost_binpacking.cc
	${COMPILER} -c ${SRC}/cost_binpacking.cc -o ${LIB}/cost_binpacking.o \
		-I${GECODE_INCLUDE} -I${BOOST_INCLUDE} -I${INCLUDE}

${LIB}/propagate.o: ${LIB}/dag_pack.o ${SRC}/propagate.cc
	${COMPILER} -c ${SRC}/propagate.cc -o ${LIB}/propagate.o \
		-I${GECODE_INCLUDE} -I${BOOST_INCLUDE} -I${INCLUDE} 

${LIB}/cost_multibin.o: ${LIB}/propagate_multi.o ${LIB}/dag_pack.o ${SRC}/cost_multibin.cc
	${COMPILER} -c ${SRC}/cost_multibin.cc -o ${LIB}/cost_multibin.o \
		-I${GECODE_INCLUDE} -I${BOOST_INCLUDE} -I${INCLUDE} -I${GUROBI_INC} -I${CLIQUER_INC}

${LIB}/propagate_multi.o: ${LIB}/dag_pack.o ${SRC}/propagate_multi.cc
	${COMPILER} -c ${SRC}/propagate_multi.cc -o ${LIB}/propagate_multi.o \
		-I${QSOPT} -I${GECODE_INCLUDE} -I${BOOST_INCLUDE} -I${INCLUDE} -I${GUROBI_INC} -I${CLIQUER_INC}

${LIB}/path.o: ${LIB}/dag_pack.o ${SRC}/path.cc
	${COMPILER} -c ${SRC}/path.cc -o ${LIB}/path.o \
		-I${QSOPT} -I${BOOST_INCLUDE} -I${INCLUDE} -I${GECODE_INCLUDE} -I${GUROBI_INC}

${LIB}/dag_pack.o: ${SRC}/dag_pack.cc
	${COMPILER} -c ${SRC}/dag_pack.cc -o ${LIB}/dag_pack.o \
		-I${QSOPT} -I${BOOST_INCLUDE} -I${INCLUDE} -I${GECODE_INCLUDE} -I${GUROBI_INC}

clean::
	rm -f *.o
	rm -f ${LIB}/*.o
	rm -f *~
	rm -f ${SRC}/*~ ${INCLUDE}/*~
