# Makefile for all programs relating to scoreassoc

# Destination for executables, change this if you want
DCBIN = ../bin

C = gcc
CC = g++

MAX_LOCI_MAK = 12000
MAX_ALL_MAK = 40
MAX_SUB_MAK = 10000

OURFLAGS = $(CFLAGS) -DMAX_LOCI=$(MAX_LOCI_MAK) -DMAX_ALL=$(MAX_ALL_MAK) -DMAX_SUB=$(MAX_SUB_MAK) 

HEADERS = cdflib.h  dcerror.hpp  dcexpr.hpp  fisher.h  sagcutils.h  safilterfuncs.hpp  scoreassoc.hpp
# cheat and just assume all code dependent on all of these

ifdef INOBJ
all: scoreassoc pathwayAssoc permPathwayAssoc 
else
all:
	if [ ! -e ../obj ] ; then mkdir ../obj ; fi ; \
	if [ ! -e ${DCBIN} ] ; then mkdir ${DCBIN} ; fi ; \
	cd ../obj; \
	make -f ../scoreassocCode/scoreassoc.mak INOBJ=INOBJ ; \
	cp scoreassoc pathwayAssoc permPathwayAssoc ${DCBIN} ; \
	echo copied executables to ${DCBIN} ; \
	cd ../src
endif

clean:
	rm ../obj/*.o

VPATH=../src
	
%.o: ../src/%.cpp $(HEADERS)
	$(CC) $(OURFLAGS) -c $< -o ../obj/$@
	
%.o: ../src/%.c $(HEADERS)
	$(C) $(OURFLAGS) -c $< -o ../obj/$@

scoreassoc: scoreassoc.o saglobals.o scoreassocfuncs.o sarecfuncs.o sahaprecfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o 
	$(CC) -o scoreassoc scoreassoc.o saglobals.o scoreassocfuncs.o sarecfuncs.o sahaprecfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o -lm

pathwayAssoc: pathwayAssoc.o dcdflib.o ipmpar.o dcerror.o 
	$(CC) -o pathwayAssoc pathwayAssoc.o dcdflib.o ipmpar.o dcerror.o -lm

permPathwayAssoc: permPathwayAssoc.o dcdflib.o ipmpar.o dcerror.o 
	$(CC) -o permPathwayAssoc permPathwayAssoc.o dcdflib.o ipmpar.o dcerror.o -lm

