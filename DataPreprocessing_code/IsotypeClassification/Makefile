###Makefile for LocalAlign
###By Feng 03/08/2014 @bu
#######v1.0

###Define variables###### 
#GCC= gcc  <------
GCC = ${CC}
#GXX = ${GXX}  
GXX = g++ 

#CFLAG= ${CPPFLAGS} -Wall -g -Werror -O -std=c++11
CFLAG=  -Wall -g -Werror -O -std=c++11  

#LOADFLAG=  -s -lm -lz ${LDFLAGS} -L/home/ffeng/miniconda3/lib
LOADFLAG= -s -lm -lz 
CXXFLAG=${CFLAG}

ACCESSDIR=Accessory/


SRCS_4=	NGSMapping_Isotype_main.cpp ${ACCESSDIR}string_ext.cpp score.cpp ${ACCESSDIR}SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp ${ACCESSDIR}FastaHandler.cpp SequenceHandlerIsotype.cpp TracebackTable.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp LocalAlignment.cpp SequenceHandlerCommon.cpp ${ACCESSDIR}GzTools.cpp ${ACCESSDIR}FastqHandler.cpp ${ACCESSDIR}FileHandler.cpp ${ACCESSDIR}FASTQ.cpp Path.cpp


OBJS_4=${SRCS_4:.cpp=.o}

PROG_4=ngsmapping_Isotype
DEPEND=$(GXX) $(CFLAG) -MM

######Rules######

all: $(PROG_4) 

.PHONY: clean all depend

clean:
	rm -fr *.o *~ core $(PROG_4) 
	cd $(ACCESSDIR) && rm -fr *.o *~ core; #ls ;   #pwd ; # note makefile is a script and each line (of commands) is doing its own sub-precess (&& is trying to stop if the "cd" command fails.

.cpp.o:   #old fasion suffix rule, double suffix rule
	$(GXX) $(CXXFLAG) -c $< -o $(addsuffix .o, $(basename $<))


$(PROG_4): $(OBJS_4)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
	@echo "";	
	@echo "******Make complete"


depend: .depend

.depend: Makefile $(SRCS_4) 
	$(GXX) -MM *.cpp >.depend
	@echo " "
	@echo "****Dependencies generated successfully."


include .depend 

