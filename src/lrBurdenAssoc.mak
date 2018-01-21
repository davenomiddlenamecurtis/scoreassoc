# Makefile for all programs relating to scoreassoc

# Destination for executables, change this if you want
DCBIN = ../bin
C = gcc
CC = g++

MAX_LOCI_MAK = 12000
MAX_ALL_MAK = 40
MAX_SUB_MAK = 15000

MYFLAGS = -DMAX_LOCI=$(MAX_LOCI_MAK) -DMAX_ALL=$(MAX_ALL_MAK) -DMAX_SUB=$(MAX_SUB_MAK) -std=gnu++0x
OURFLAGS = $(MYFLAGS) $(EXTRAFLAGS)
CFLAGS = $(OURFLAGS) ${DEBUGFLAG}

# so to compile for debugging use make -f scoreassoc.mak DEBUGFLAG=-g

HEADERS = cdflib.h  dcerror.hpp  dcexpr.hpp  fisher.h  sagcutils.h  safilterfuncs.hpp  lrBAFilterFuncs.hpp lrBurdenAsssoc.hpp lrModel.hpp
# cheat and just assume all code dependent on all of these

EXES = lrBurdenAssoc # combineGeneScores getVarScores combineCohorts 
DLIB = /home/rejudcu/dlib-19.4

ifdef INOBJ
all: ${EXES}
else
all:
	if [ ! -e ../obj ] ; then mkdir ../obj ; fi ; \
	if [ ! -e ${DCBIN} ] ; then mkdir ${DCBIN} ; fi ; \
	cd ../obj; \
	make -f ../src/lrBurdenAssoc.mak INOBJ=INOBJ CFLAGS=$(CFLAGS) ; \
	cp ${EXES} ${DCBIN} ; \
	echo copied executables to ${DCBIN} ; \
	cd ../src
endif
 
clean:
	rm ../obj/*

VPATH=../src
	
%.o: ../src/%.c $(HEADERS)
	$(C) $(CFLAGS) -I ${DLIB} -c $< -o ../obj/$@ 
	
lrBurdenAssoc: lrBurdenAssoc.o lrBurdenAssocGlobals.o lrBurdenAssocFuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o lrBAFilterFuncs.o lrBAInit.o lrModel.o
	$(CC) ${DEBUGFLAG} -o lrBurdenAssoc lrBurdenAssoc.o lrBurdenAssocGlobals.o lrBurdenAssocFuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o lrBAFilterFuncs.o lrBAInit.o lrModel.o -lm

