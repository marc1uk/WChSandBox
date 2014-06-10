#!/bin/env python
import ROOT  # @UnresolvedImport
ROOT.PyConfig.IgnoreCommandLineOptions = True
import pkg_resources
import argparse

###############################################################################

def chktrack2(inputfile, outputfile):
#   Int_t i, j;
#   Int_t nevents;
#   double cosx,cosy,cosz;
#   double p;
#   Int_t stat;
    kPi = ROOT.TMath.Pi();
  
#   TTree  *tn;
#   NeutVtx *nvtx;
#   NeutVect *nvect;

    ran = ROOT.TRandom3(0);

    f = ROOT.TFile(inputfile);
    tn = f.Get("nRooTracker"); # try to get NEUT tree
    if not tn:
        tn = f.Get("gRooTracker"); # try to get GENIE tree
    nevents = tn.GetEntries();
    fcomp = open(outputfile,"w");
    for j, entry in enumerate(tn):
        m_to_cm = 100.0
        GeV_to_MeV = 1000.0
        try:
            mode = entry.NEneutmode
        except AttributeError:
            mode = entry.G2NeutEvtCode
        vtx = entry.EvtVtx
    
        x = vtx[0] * m_to_cm
        y = vtx[1] * m_to_cm
        z = vtx[2] * m_to_cm
        t = vtx[3] * m_to_cm
    
        fcomp.write(" $begin %d\n" % j);
        fcomp.write(" $nuance %d\n" % mode);
        fcomp.write(" $vertex %5.4f %5.4f %5.4f %5.4f\n" % (x, y, z, t) );

        npart = entry.StdHepN
        stdhepp4 = ROOT.getRooTrackerHepP4(entry.StdHepN, entry.StdHepP4)
        tracklist = []
        infoline_index = 0
        for i in xrange(npart):
            #get the information that we need
            pvec = stdhepp4.at(i)
            pid = entry.StdHepPdg[i]
            stdhepstatus = entry.StdHepStatus[i]
            #p = ROOT.TMath.Sqrt(
            #                    pvec.E()*(pvec.E()-part.fMass*part.fMass)
            #                    )
            pvec = ROOT.TLorentzVector(pvec[0] * GeV_to_MeV, 
                                       pvec[1] * GeV_to_MeV, 
                                       pvec[2] * GeV_to_MeV, 
                                       pvec[3] * GeV_to_MeV
                                       )
            p = pvec.P()
            if p == 0: 
                p = 0.000001;
            cosx = pvec.Px() / p;
            cosy = pvec.Py() / p;
            cosz = pvec.Pz() / p;

            if stdhepstatus==0:
                #initial state particle
                stat = -1;
            elif stdhepstatus == 1:
                #final state particle
                stat = 0
            else:
                #unknown (intermediate?)
                stat = -2;

            if (stat == -2):
                continue; # get rid of intermediate particles

            if stat == -1:
                infoline_index = i

            if pid > 10000:
                continue;
            tracklist.append(" $track %d %5.4f %5.4f %5.4f %5.4f %d\n" % (pid, pvec.E(), cosx, cosy, cosz, stat))
        #insert the info line into the tracklist
        #infoline = " $info 0 0 %d\n" % j
        #tracklist.insert(infoline_index, infoline)
        for trk in tracklist:
            fcomp.write(trk)
        #fcomp.write(" $headerend %d\n" % j);
        fcomp.write(" $end %d\n" % j);
    fcomp.close()
    return

###############################################################################

def loadlib():
    srcfile = pkg_resources.resource_filename("vectorgen", "readRooTracker.C+")
    ROOT.gROOT.ProcessLine(".L " + srcfile)
    return

###############################################################################

def parsecml():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=str, help="Input genev file name.")
    parser.add_argument("-o", "--outfile", dest="outfile", type=str, help="Output file name genev file name.", default=None)
    return parser.parse_args()

###############################################################################

def getfilenames(opt):
    inname = opt.infile #"genev_nu_nd280_r200_z400_Z_10000.root"
    outname = opt.outfile #"nuance_nu_nd280_r200_z400_Z_10000.txt"
    if outname is None:
        #automatically determine output file name
        if not ("genev" in inname and ".root" in inname):
            raise Exception("Cannot automatically determine outputfile name from input. You must specify it manually with the -")
        outname = inname.replace("genev", "nuance").replace(".root", ".dat")
    return inname, outname

###############################################################################

def main():
    loadlib()
    opt = parsecml()
    inname, outname = getfilenames(opt)
    chktrack2(inname, outname)
    return

###############################################################################

if __name__ == "__main__":
    main()

