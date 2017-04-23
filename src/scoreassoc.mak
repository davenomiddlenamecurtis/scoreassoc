# Makefile for all programs relating to scoreassoc

# Destination for executables, change this if you want
DCBIN = ../bin

C = gcc
CC = g++

MAX_LOCI_MAK = 12000
MAX_ALL_MAK = 40
MAX_SUB_MAK = 15000

OURFLAGS = $(CFLAGS) -DMAX_LOCI=$(MAX_LOCI_MAK) -DMAX_ALL=$(MAX_ALL_MAK) -DMAX_SUB=$(MAX_SUB_MAK) 
# DEBUGFLAG = -g

HEADERS = cdflib.h  dcerror.hpp  dcexpr.hpp  fisher.h  sagcutils.h  safilterfuncs.hpp  scoreassoc.hpp
# cheat and just assume all code dependent on all of these

EXES = combineGeneScores getVarScores scoreassoc pathwayAssoc permPathwayAssoc combineCohorts makeScoreTable 

ifdef INOBJ
all: ${EXES}
else
all:
	if [ ! -e ../obj ] ; then mkdir ../obj ; fi ; \
	if [ ! -e ${DCBIN} ] ; then mkdir ${DCBIN} ; fi ; \
	cd ../obj; \
	make -f ../src/scoreassoc.mak INOBJ=INOBJ ; \
	cp ${EXES} ${DCBIN} ; \
	echo copied executables to ${DCBIN} ; \
	cd ../src
endif
 
clean:
	rm ../obj/*

VPATH=../src
	
%.o: ../src/%.cpp $(HEADERS)
	$(CC) $(OURFLAGS) ${DEBUGFLAG} -c $< -o ../obj/$@
	
%.o: ../src/%.c $(HEADERS)
	$(C) $(OURFLAGS) ${DEBUGFLAG} -c $< -o ../obj/$@
	
combineGeneScores: combineGeneScores.o dcerror.o
	$(CC) ${DEBUGFLAG} -o combineGeneScores combineGeneScores.o dcerror.o -lm

getVarScores: getVarScores.o saglobals.o scoreassocfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o 
	$(CC) ${DEBUGFLAG} -o getVarScores getVarScores.o saglobals.o scoreassocfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o -lm

scoreassoc: scoreassoc.o saglobals.o scoreassocfuncs.o sarecfuncs.o sahaprecfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o 
	$(CC) ${DEBUGFLAG} -o scoreassoc scoreassoc.o saglobals.o scoreassocfuncs.o sarecfuncs.o sahaprecfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o -lm

combineCohorts: combineCohorts.o dcdflib.o ipmpar.o dcerror.o 
	$(CC) ${DEBUGFLAG} -o combineCohorts combineCohorts.o dcdflib.o ipmpar.o dcerror.o -lm

makeScoreTable: makeScoreTable.o dcerror.o 
	$(CC) ${DEBUGFLAG} -o makeScoreTable makeScoreTable.o dcerror.o -lm

pathwayAssoc: pathwayAssoc.o dcdflib.o ipmpar.o dcerror.o 
	$(CC) ${DEBUGFLAG} -o pathwayAssoc pathwayAssoc.o dcdflib.o ipmpar.o dcerror.o -lm

permPathwayAssoc: permPathwayAssoc.o dcdflib.o ipmpar.o dcerror.o 
	$(CC) ${DEBUGFLAG} -o permPathwayAssoc permPathwayAssoc.o dcdflib.o ipmpar.o dcerror.o -lm

