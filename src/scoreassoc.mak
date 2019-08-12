# Makefile for all programs relating to scoreassoc
# You must edit the line below so that it points to the righht folder for your installation of DLIB (obtainable from dlib.net)
DLIB = /home/rejudcu/dlib-19.4

# Destination for executables, change this if you want
DCBIN = ../bin
# on cluster must enter this to get right compiler versiion:
# scl enable devtoolset-3 bash
# DLIB will only compile with later versions of C/C++

C = gcc
CC = g++

MAX_LOCI_MAK = 12000
MAX_ALL_MAK = 40
MAX_SUB_MAK = 15000

MYFLAGS = $(CFLAGS) -DMAX_LOCI=$(MAX_LOCI_MAK) -DMAX_ALL=$(MAX_ALL_MAK) -DMAX_SUB=$(MAX_SUB_MAK)  # -std=gnu++0x
CPPFLAGS = -std=c++14 
OURFLAGS = $(MYFLAGS) $(EXTRAFLAGS)

# so to compile for debugging use make -f scoreassoc.mak DEBUGFLAG=-g

HEADERS = cdflib.h  dcerror.hpp  dcexpr.hpp  fisher.h  sagcutils.h  safilterfuncs.hpp  scoreassoc.hpp
# cheat and just assume all code dependent on all of these

EXES = scoreassoc pathwayAssoc permPathwayAssoc makeScoreTable # fitScores combineGeneScores getVarScores combineCohorts 

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
	$(CC) $(OURFLAGS) ${DEBUGFLAG} ${CPPFLAGS} -c $< -o ../obj/$@
	
%.o: ../src/%.c $(HEADERS)
	$(C) $(OURFLAGS) ${DEBUGFLAG} -c $< -o ../obj/$@
	
glModel.o: ../src/glModel.cpp ${DLIB}/dlib/optimization.h
	$(CC) $(OURFLAGS) $(CPPFLAGS) ${DEBUGFLAG} -c ../src/glModel.cpp  -o ../obj/glModel.o -I ${DLIB}
	
lrModel.o: ../src/lrModel.cpp ${DLIB}/dlib/optimization.h
	$(CC) $(OURFLAGS) $(CPPFLAGS) ${DEBUGFLAG} -c ../src/lrModel.cpp  -o ../obj/lrModel.o -I ${DLIB}
	
combineGeneScores: combineGeneScores.o dcerror.o
	$(CC) ${DEBUGFLAG} -o combineGeneScores combineGeneScores.o dcerror.o -lm

fitScoresWithDlib.o: ../src/fitScores.cpp ${DLIB}/dlib/optimization.h
	$(CC) $(OURFLAGS) $(CPPFLAGS) ${DEBUGFLAG} -c ../src/fitScores.cpp  -o ../obj/fitScoresWithDlib.o -I ${DLIB}

fitScores: fitScoresWithDlib.o dcerror.o
	$(CC) ${DEBUGFLAG} -o fitScores ../obj/fitScoresWithDlib.o ../obj/dcerror.o -lm

getVarScores: getVarScores.o saglobals.o scoreassocfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o 
	$(CC) ${DEBUGFLAG} -o getVarScores getVarScores.o saglobals.o scoreassocfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o -lm

scoreassoc: scoreassoc.o sainit.o salrtests.o glModel.o saglobals.o scoreassocfuncs.o sarecfuncs.o sahaprecfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o 
	$(CC) ${DEBUGFLAG} -o scoreassoc scoreassoc.o sainit.o salrtests.o glModel.o saglobals.o scoreassocfuncs.o sarecfuncs.o sahaprecfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o -lm

combineCohorts: combineCohorts.o dcdflib.o ipmpar.o dcerror.o 
	$(CC) ${DEBUGFLAG} -o combineCohorts combineCohorts.o dcdflib.o ipmpar.o dcerror.o -lm

makeScoreTable: makeScoreTable.o dcerror.o 
	$(CC) ${DEBUGFLAG} -o makeScoreTable makeScoreTable.o dcerror.o -lm

pathwayAssoc: pathwayAssoc.o dcdflib.o ipmpar.o dcerror.o salrtests.o glModel.o 
	$(CC) ${DEBUGFLAG} -o pathwayAssoc pathwayAssoc.o dcdflib.o ipmpar.o dcerror.o salrtests.o glModel.o -lm

permPathwayAssoc: permPathwayAssoc.o dcdflib.o ipmpar.o dcerror.o 
	$(CC) ${DEBUGFLAG} -o permPathwayAssoc permPathwayAssoc.o dcdflib.o ipmpar.o dcerror.o -lm

