## Makefile ##

#
# notes:
#  WCSIM_INCDIR  points to WCSim .hh files
#  WCSIM_LIBDIR  points to WCSim .o files
#

CXX           = g++
CXXDEPEND     = -MM
CXXFLAGS      = -g -Wall -fPIC -O0
LD            = g++
LDFLAGS       = -g 

UNAME := $(shell uname)

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLDFLAGS  := $(shell root-config --ldflags)
ROOTLIBS     := $(shell root-config --evelibs) 
# $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXXFLAGS     += $(ROOTCFLAGS)
LDFLAGS      += $(ROOTLDFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

INCDIR = ./include
SRCDIR = ./src
TMPDIR = ./tmp
LIBDIR = ./lib
BINDIR = ./bin

INCLUDES = -I$(INCDIR)

#WCSIM_INCDIR = $(WCSIM)/wcsim/include
#WCSIM_LIBDIR = $(WCSIM)/x86_64-slc5-gcc43-dbg
WCSIM_INCDIR = $(WCSIM)/include
WCSIM_LIBDIR = $(WCSIM)/tmp

WCSIM_INCLUDES = -I$(WCSIM_INCDIR)
WCSIM_LDFLAGS += -L$(WCSIM_LIBDIR)
WCSIM_LIBS += -lWCSim

.PHONY: 
	all

all: clean rootcint shared

ROOTSO := $(LIBDIR)/libWCLAnalysis.so

ROOTDICT := $(SRCDIR)/WCLRootDict.cc

ROOTSRC := $(SRCDIR)/WCSimRecoObjectTable.cc $(INCDIR)/WCSimRecoObjectTable.hh $(SRCDIR)/WCSimRecoDigit.cc $(INCDIR)/WCSimRecoDigit.hh $(SRCDIR)/WCSimRecoCluster.cc $(INCDIR)/WCSimRecoCluster.hh $(SRCDIR)/SandBoxPMTcoverage.cc $(INCDIR)/SandBoxPMTcoverage.hh $(SRCDIR)/WCSimTrueLight.cc $(INCDIR)/WCSimTrueLight.hh $(SRCDIR)/WCSimTruePart.cc $(INCDIR)/WCSimTruePart.hh $(SRCDIR)/WCSimTrueCapture.cc $(INCDIR)/WCSimTrueCapture.hh $(SRCDIR)/WCSimTrueLightCluster.cc $(INCDIR)/WCSimTrueLightCluster.hh $(SRCDIR)/WCLTreeReader.cc $(INCDIR)/WCLTreeReader.hh $(SRCDIR)/WCLTreeWriter.cc $(INCDIR)/WCLTreeWriter.hh $(INCDIR)/WCLRootLinkDef.hh

ROOTOBJS := $(TMPDIR)/WCSimRecoObjectTable.o $(TMPDIR)/WCSimRecoDigit.o $(TMPDIR)/WCSimRecoCluster.o $(TMPDIR)/SandBoxPMTcoverage.o $(TMPDIR)/WCSimTrueLight.o $(TMPDIR)/WCSimTruePart.o $(TMPDIR)/WCSimTrueCapture.o $(TMPDIR)/WCSimTrueLightCluster.o $(TMPDIR)/WCLTreeReader.o $(TMPDIR)/WCLTreeWriter.o  $(TMPDIR)/WCLRootDict.o

$(TMPDIR)/%.o : $(SRCDIR)/%.cc
	@echo "<**Compiling $@**>"
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(WCSIM_INCLUDES) -c $< -o $@

$(TMPDIR)/%.d: $(SRCDIR)/%.cc
	@echo "<**Depend $@**>"
	set -e; $(CXX) $(CXXDEPEND) $(CXXFLAGS) $(INCLUDES) $(WCSIM_INCLUDES) $< \
	| sed 's!$*\.o!& $@!' >$@;\
	[ -s $@ ] || rm -f $@

$(ROOTDICT) : $(ROOTSRC)

rootcint : $(ROOTSRC)
	@echo "<**Rootcint**>"
	rootcint -f $(ROOTDICT) -c -I$(INCDIR) -I$(WCSIM_INCDIR) -I$(shell root-config --incdir) WCSimRecoObjectTable.hh WCSimRecoDigit.hh WCSimRecoCluster.hh SandBoxPMTcoverage.hh WCSimTrueLight.hh WCSimTruePart.hh WCSimTrueCapture.hh WCSimTrueLightCluster.hh WCLTreeReader.hh WCLTreeWriter.hh WCLRootLinkDef.hh

shared: $(ROOTDICT) $(ROOTSRC) $(ROOTOBJS)
	@echo "<**Shared**>"
ifeq ($(UNAME), Darwin) 
	g++ -dynamiclib $(ROOTLIBS) $(ROOTGLIBS) -O $(ROOTOBJS) -lMinuit $(WCSIM)/libWCSimRoot.so -o $(ROOTSO)
endif
ifeq ($(UNAME), Linux) 
	g++ -shared $(ROOTLIBS) $(ROOTGLIBS) -O $(ROOTOBJS) -o $(ROOTSO)
endif

clean :
	@echo "<**Clean**>"
	rm -f $(SRCDIR)/*~ $(INCDIR)/*~ $(TMPDIR)/*.o $(TMPDIR)/*.d $(TMPDIR)/*.a $(LIBDIR)/*.so $(SRCDIR)/WCSimAnalysisRootDict.*

reco.o : reco.cc 
	@echo "making $@ from $<"
	@echo $(CXX) $(CXXFLAGS) $(INCLUDES) $(WCSIM_INCLUDES) -c $< -o $@
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(WCSIM_INCLUDES) -c $< -o $@

reco : reco.o $(ROOTSO)
	@echo "Making reco executable for WCSimAnalysis..."
	g++ -o reco $(CXXFLAGS) reco.o $(WCSIM)/libWCSimRoot.so $(ROOTSO) $(ROOTLIBS) $(ROOTGLIBS) -lMinuit 

DEPS = $(ROOTOBJS:$(TMPDIR)/%.o=$(TMPDIR)/%.d)

ifeq ($(MAKECMDGOALS),all)
 include $(DEPS)
endif

ifeq ($(MAKECMDGOALS),shared)
 include $(DEPS)
endif
