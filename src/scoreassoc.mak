# Makefile for all programs relating to scoreassoc

# Please note that this makefile call itself. If you save a local copy for yourself then you must change the below accordingly.
# For example, if you change this file and call it mymakefile.mak, you also need to change the command below to read:
#    make -f ../src/mymakefile.mak INOBJ=INOBJ ; \  
# The line needs to finish with a backspace character with no spaces after it.

# You must edit the line below so that it points to the right folder for your installation of DLIB (obtainable from dlib.net)
DLIB = /home/dcurtis/dlib/dlib-19.19

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
MAXLRVARIABLES_MAK = 100

MYFLAGS = $(CFLAGS) -DMAX_LOCI=$(MAX_LOCI_MAK) -DMAX_ALL=$(MAX_ALL_MAK) -DMAX_SUB=$(MAX_SUB_MAK) -DMAXLRVARIABLES=$(MAXLRVARIABLES_MAK)
DLIBFLAGS = -DDLIB_NO_GUI_SUPPORT=1 -I ${DLIB} 
# BLAS and LAPACK are not used by default unless cmake is used for installation

CPPFLAGS = -std=c++14 
OURFLAGS = $(MYFLAGS) $(EXTRAFLAGS)

# to compile for debugging use make -f scoreassoc.mak DEBUGFLAG=-g 

HEADERS = cdflib.h  dcerror.hpp  dcexpr.hpp  fisher.h  sagcutils.h  safilterfuncs.hpp  scoreassoc.hpp
# cheat and just assume all code dependent on all of these

EXES = scoreassoc pathwayAssoc permPathwayAssoc makeScoreTable testGeneSets # fitScores combineGeneScores getVarScores combineCohorts 

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
# need to change the line above if the name of this makefile is changed
 
clean:
	rm ../obj/*

VPATH=../src
	
%.o: ../src/%.cpp $(HEADERS)
	$(CC) $(OURFLAGS) ${DEBUGFLAG} ${CPPFLAGS} -c $< -o ../obj/$@
	
%.o: ../src/%.c $(HEADERS)
	$(C) $(OURFLAGS) ${DEBUGFLAG} -c $< -o ../obj/$@

source.o: ${DLIB}/dlib/all/source.cpp
	$(CC) $(OURFLAGS) $(CPPFLAGS) ${DEBUGFLAG} ${DLIBFLAGS} -c ${DLIB}/dlib/all/source.cpp  -o ../obj/source.o
	
glModel.o: ../src/glModel.cpp ${DLIB}/dlib/optimization.h
	$(CC) $(OURFLAGS) $(CPPFLAGS) ${DEBUGFLAG} ${DLIBFLAGS} -c ../src/glModel.cpp  -o ../obj/glModel.o
	
lrModel.o: ../src/lrModel.cpp ${DLIB}/dlib/optimization.h
	$(CC) $(OURFLAGS) $(CPPFLAGS) ${DEBUGFLAG} ${DLIBFLAGS} -c ../src/lrModel.cpp  -o ../obj/lrModel.o
	
combineGeneScores: combineGeneScores.o dcerror.o
	$(CC) ${DEBUGFLAG} -o combineGeneScores combineGeneScores.o dcerror.o -lm

fitScoresWithDlib.o: ../src/fitScores.cpp ${DLIB}/dlib/optimization.h
	$(CC) $(OURFLAGS) $(CPPFLAGS) ${DEBUGFLAG} ${DLIBFLAGS} -c ../src/fitScores.cpp  -o ../obj/fitScoresWithDlib.o

testGeneSets: testGeneSets.o dcerror.o dcdflib.o ipmpar.o
	$(CC) ${DEBUGFLAG} -o testGeneSets testGeneSets.o dcerror.o dcdflib.o ipmpar.o -lm

fitScores: fitScoresWithDlib.o dcerror.o source.o
	$(CC) ${DEBUGFLAG} -o fitScores ../obj/fitScoresWithDlib.o ../obj/dcerror.o source.o -lm -lpthread

getVarScores: getVarScores.o saglobals.o scoreassocfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o 
	$(CC) ${DEBUGFLAG} -o getVarScores getVarScores.o saglobals.o scoreassocfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o -lm

scoreassoc: source.o scoreassoc.o sainit.o salrtests.o glModel.o saglobals.o scoreassocfuncs.o sarecfuncs.o sahaprecfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o 
	$(CC) ${DEBUGFLAG} -o scoreassoc source.o scoreassoc.o sainit.o salrtests.o glModel.o saglobals.o scoreassocfuncs.o sarecfuncs.o sahaprecfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o -lm -lpthread

combineCohorts: combineCohorts.o dcdflib.o ipmpar.o dcerror.o 
	$(CC) ${DEBUGFLAG} -o combineCohorts combineCohorts.o dcdflib.o ipmpar.o dcerror.o -lm

makeScoreTable: makeScoreTable.o dcerror.o 
	$(CC) ${DEBUGFLAG} -o makeScoreTable makeScoreTable.o dcerror.o -lm

pathwayAssoc: pathwayAssoc.o dcdflib.o ipmpar.o dcerror.o salrtests.o glModel.o source.o
	$(CC) ${DEBUGFLAG} -o pathwayAssoc pathwayAssoc.o dcdflib.o ipmpar.o dcerror.o salrtests.o glModel.o source.o -lm -lpthread

permPathwayAssoc: permPathwayAssoc.o dcdflib.o ipmpar.o dcerror.o 
	$(CC) ${DEBUGFLAG} -o permPathwayAssoc permPathwayAssoc.o dcdflib.o ipmpar.o dcerror.o -lm

