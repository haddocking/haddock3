#####################################
#                                   #
#     HADDOCK tools makefile        #
#                                   #
#####################################

CPP=g++
CPFLAGS=-O2

EXEC=contact_fcc

all: 
	make $(EXEC)

contact_fcc: contact_fcc.cpp
	$(CPP) $(CPFLAGS) -o contact_fcc contact_fcc.cpp
	chmod u+x contact_fcc

clean :
	/bin/rm $(EXEC)

