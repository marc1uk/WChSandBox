# $Id: GNUmakefile,v 1.2 2003/01/23 15:31:39 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := WChSandBox
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

CPPFLAGS += $(shell root-config --cflags)
LDFLAGS  += $(shell root-config --libs)

EXTRALIBS = $(shell $(ROOTSYS)/bin/root-config --glibs)

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
LDFLAGS := $(LDFLAGS) -L$(G4LIB)/$(G4SYSTEM)
visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
