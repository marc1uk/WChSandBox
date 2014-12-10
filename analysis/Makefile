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

ROOTSRC := $(SRCDIR)/WCSimRecoObjectTable.cc $(INCDIR)/WCSimRecoObjectTable.hh $(SRCDIR)/WCSimRecoDigit.cc $(INCDIR)/WCSimParameters.hh $(SRCDIR)/WCSimParameters.cc $(SRCDIR)/EventDisplay2D.cc $(INCDIR)/EventDisplay2D.hh $(SRCDIR)/EventDisplay3D.cc $(INCDIR)/EventDisplay3D.hh $(INCDIR)/WCSimWaterModel.hh $(SRCDIR)/WCSimWaterModel.cc $(INCDIR)/WCSimRecoDigit.hh $(SRCDIR)/WCSimRecoCluster.cc $(INCDIR)/WCSimRecoCluster.hh $(SRCDIR)/SandBoxPMTcoverageBox.cc $(INCDIR)/SandBoxPMTcoverageBox.hh $(SRCDIR)/SandBoxPMTcoverage.cc $(INCDIR)/SandBoxPMTcoverage.hh $(SRCDIR)/WCSimIsochron.cc $(INCDIR)/WCSimIsochron.hh $(SRCDIR)/WCSimVertexGeometry.cc $(INCDIR)/WCSimVertexGeometry.hh $(SRCDIR)/WCSimTrueLight.cc $(INCDIR)/WCSimTrueLight.hh $(SRCDIR)/WCSimTruePart.cc $(INCDIR)/WCSimTruePart.hh $(SRCDIR)/WCSimTrueCapture.cc $(INCDIR)/WCSimTrueCapture.hh $(SRCDIR)/WCSimTrueLightCluster.cc $(INCDIR)/WCSimTrueLightCluster.hh $(SRCDIR)/WCLTreeReader.cc $(INCDIR)/WCLTreeReader.hh $(SRCDIR)/WCLTreeWriter.cc $(INCDIR)/WCLTreeWriter.hh $(SRCDIR)/RasterizeCircle.cc $(INCDIR)/RasterizeCircle.hh $(SRCDIR)/WChRecoLite.cc $(INCDIR)/WChRecoLite.hh $(SRCDIR)/WCSimIsoChronTransform.cc $(INCDIR)/WCSimIsoChronTransform.hh $(INCDIR)/WCLRootLinkDef.hh

ROOTOBJS := $(TMPDIR)/WCSimRecoObjectTable.o $(TMPDIR)/WCSimRecoDigit.o $(TMPDIR)/WCSimParameters.o $(TMPDIR)/EventDisplay2D.o $(TMPDIR)/EventDisplay3D.o $(TMPDIR)/WCSimWaterModel.o $(TMPDIR)/WCSimRecoCluster.o $(TMPDIR)/SandBoxPMTcoverageBox.o $(TMPDIR)/SandBoxPMTcoverage.o $(TMPDIR)/WCSimIsochron.o $(TMPDIR)/WCSimVertexGeometry.o $(TMPDIR)/WCSimTrueLight.o $(TMPDIR)/WCSimTruePart.o $(TMPDIR)/WCSimTrueCapture.o $(TMPDIR)/WCSimTrueLightCluster.o $(TMPDIR)/WCLTreeReader.o $(TMPDIR)/WCLTreeWriter.o $(TMPDIR)/RasterizeCircle.o $(TMPDIR)/WChRecoLite.o $(TMPDIR)/WCSimIsoChronTransform.o $(TMPDIR)/WCLRootDict.o

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
	rootcint -f $(ROOTDICT) -c -I$(INCDIR) -I$(WCSIM_INCDIR) -I$(shell root-config --incdir) WCSimRecoObjectTable.hh WCSimRecoDigit.hh WCSimParameters.hh EventDisplay2D.hh  EventDisplay3D.hh WCSimWaterModel.hh WCSimIsochron.hh WCSimRecoCluster.hh SandBoxPMTcoverageBox.hh  SandBoxPMTcoverage.hh WCSimVertexGeometry.hh WCSimTrueLight.hh WCSimTruePart.hh WCSimTrueCapture.hh WCSimTrueLightCluster.hh WCLTreeReader.hh WCLTreeWriter.hh RasterizeCircle.hh  WChRecoLite.hh WCSimIsoChronTransform.hh WCLRootLinkDef.hh

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
