# Makefile for all programs relating to scoreassoc

# Destination for executables, change this if you want
DCBIN = ../bin

# on cluster must enter this to get right compiler version:
#  scl enable devtoolset-3 bash

C = gcc
CC = g++

MAX_LOCI_MAK = 12000
MAX_ALL_MAK = 40
MAX_SUB_MAK = 15000

MYFLAGS = $(CFLAGS) -DMAX_LOCI=$(MAX_LOCI_MAK) -DMAX_ALL=$(MAX_ALL_MAK) -DMAX_SUB=$(MAX_SUB_MAK) -std=c++14 
OURFLAGS = $(MYFLAGS) $(EXTRAFLAGS)

# so to compile for debugging use make -f scoreassoc.mak DEBUGFLAG=-g

HEADERS = cdflib.h  dcerror.hpp  dcexpr.hpp  fisher.h  sagcutils.h  safilterfuncs.hpp  scoreassoc.hpp
# cheat and just assume all code dependent on all of these

EXES = lrBurdenAssoc
DLIB = /home/rejudcu/dlib-19.4
# needed only for fitScores

ifdef INOBJ
all: ${EXES}
else
all:
	if [ ! -e ../obj ] ; then mkdir ../obj ; fi ; \
	if [ ! -e ${DCBIN} ] ; then mkdir ${DCBIN} ; fi ; \
	cd ../obj; \
	make -f ../src/lrBurdenAssoc.mak INOBJ=INOBJ ; \
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
	
lrModel.o: ../src/lrModel.cpp ${DLIB}/dlib/optimization.h
	$(CC) $(OURFLAGS) ${DEBUGFLAG} -c ../src/lrModel.cpp  -o ../obj/lrModel.o -I ${DLIB}
	
lrBurdenAssoc : lrBurdenAssoc.o lrBurdenAssocGlobals.o lrBurdenAssocFuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o lrBAFilterFuncs.o lrBAInit.o lrModel.o
	$(CC) ${DEBUGFLAG} -o lrBurdenAssoc lrBurdenAssoc.o lrBurdenAssocGlobals.o lrBurdenAssocFuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o lrBAFilterFuncs.o lrBAInit.o lrModel.o -lm

