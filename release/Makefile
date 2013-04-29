## Stefano Gualandi, stefano.gualandi@gmail.com
## April, 2013

# README:
# -------
# In order to compile this example you need the following libraries:
#  - Gecode v4.0.0  (www.gecode.org) 
#  - Cliquer v1.21  (http://users.tkk.fi/pat/cliquer.html)

# Gecode directory
GECODE_INCLUDE  = /Users/gualandi/LIBS/gecode/include
GECODE_LIB      = -L/Users/gualandi/LIBS/gecode/lib \
   -lgecodedriver -lgecodeint -lgecodekernel -lgecodeminimodel \
   -lgecodesearch -lgecodeset -lgecodesupport 

# Cliquer
CLIQUER_INC = /Users/gualandi/MyLIB/cliquer-1.21/
CLIQUER_LIB = /Users/gualandi/MyLIB/cliquer-1.21/cliquer.o \
   /Users/gualandi/MyLIB/cliquer-1.21/graph.o \
   /Users/gualandi/MyLIB/cliquer-1.21/reorder.o

COMPILER    = g++ -O2 -DNDEBUG -Wall -g -fPIC
LINKER      = g++ -O2 -DNDEBUG -lm -lpthread -g

example: multibin-packing.o example.cc
	${COMPILER} -c example.cc -o example.o -I. -I${GECODE_INCLUDE} -I${CLIQUER_INC}
	${LINKER} -o example example.o multibin-packing.o ${GECODE_LIB} ${CLIQUER_LIB}

multibin-packing.o: multibin-packing.cpp
	${COMPILER} -c multibin-packing.cpp -o multibin-packing.o -I. -I${GECODE_INCLUDE} -I${CLIQUER_INC}

clean::
	rm -f *.o *~
