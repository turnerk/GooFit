include ../../Makefile.config

#------------------------------------------------------------------------------
GOOFITDIR = ${PWD}/../../

INCLUDES += -I$(GOOFITDIR) -I$(GOOFITDIR)/rootstuff -I$(GOOFITDIR)/PDFs/
LIBS += -L$(GOOFITDIR)/rootstuff -lRootUtils

# These are for user-level programs that want access to the ROOT plotting stuff, 
# not just the fitting stuff included in the GooFit-local ripped library. 
ROOT_INCLUDES = -I$(shell root-config --incdir)
ROOT_LIBS     = $(shell root-config --libs)
WRKDIR=$(GOOFITDIR)/wrkdir/

THRUSTO		= $(WRKDIR)/Variable.o $(WRKDIR)/FitManager.o $(WRKDIR)/GooPdfCUDA.o $(WRKDIR)/Faddeeva.o $(WRKDIR)/FitControl.o $(WRKDIR)/PdfBase.o $(WRKDIR)/DataSet.o $(WRKDIR)/BinnedDataSet.o $(WRKDIR)/UnbinnedDataSet.o $(WRKDIR)/FunctorWriter.o 
ROOTRIPDIR	= $(GOOFITDIR)/rootstuff
ROOTUTILLIB	= $(ROOTRIPDIR)/libRootUtils.so 

.SUFFIXES: 

all:	zach

%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES) $(DEFINEFLAGS) $(CXXFLAGS) -c -o $@ $<

zach:		ZachFit.o $(THRUSTO) $(ROOTUTILLIB) 
		$(LD) $(LDFLAGS) $^ $(LIBS) $(ROOT_LIBS) -o $@
		@echo "$@ done"

clean:
		@rm -f *.o core zach zach_mchist.png 
