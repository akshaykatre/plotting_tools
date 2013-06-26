#!/usr/bin/env python 

import ROOT
import os 
import math 
from ROOT import kRed, kBlue
import argparse
from Tools import CanvasMaker
import glob
import code
#print os.getcwd()

parser = argparse.ArgumentParser(prog="mc_vs_data", description='')

parser.add_argument('-d', '--data', dest='DATAPATH', help="Destination of the data files; along with a pattern by which to recognize them eg. $HOME/*MonFile* or just MonFile if files belong to this directory")
parser.add_argument('-mc', '--mc', dest='MCPATH', help="Destination of MC files; along with pattern by which to recognize them eg. $HOME/*MonFile* or just MonFile if the files belong to this directory")
parser.add_argument('-mcl', '--mclabel', dest='MCLABEL', default='MC12', help="Label for MC plots")
parser.add_argument('-dl', '--datalabel', dest='DATALABEL', default='Data', help="Label for data plots")
parser.add_argument('-s', '--save', dest='SAVE', action='store_true', help="Saves the histograms in .png format")
parser.add_argument('-n', '--norm', dest='NORM', action='store_true', help="Normalize histograms")
args = parser.parse_args()
#### Color Scheme 
##Red - Data plots
##Blue - MC plots
####

#print args.DATAPATH
#print glob.glob(args.DATAPATH)

files_data = []
files_mc = []
if args.DATAPATH:
    files_data = glob.glob(args.DATAPATH)
if args.MCPATH:
    files_mc = glob.glob(args.MCPATH)

# for i in range(0, 13):
#     if os.path.exists('user.akatre.000064.MonFile._%05d.root' %i):
#         files_data.append('user.akatre.000064.MonFile._%05d.root' %i)

# files_mc.append("MC12_SingleParticle.root")


chain_mc = ROOT.TChain("MonAnaTree")
for mcfiles in files_mc:
    chain_mc.Add(mcfiles)
#chain_mc.Add(files_mc[0])


chain_data = ROOT.TChain('MonAnaTree')
ROOT.gStyle.SetOptStat(0)
for files in files_data:
    chain_data.Add(files)

#for data 
clusterE = ROOT.TH1F("clusterEnergy", "clusterEnergy", 100, 0, 150)
clusterE.SetLineColor(kRed)
clusterE.SetLineWidth(4)

clusterEt = ROOT.TH1F("clusterEt", "clusterEt", 100, 0, 150)
clusterEt.SetLineColor(kRed)
clusterEt.SetLineWidth(4)

clusterEta = ROOT.TH1F("clusterEta", "ClusterEta", 100,-3, 3)
clusterEta.SetLineColor(kRed)
clusterEta.SetLineWidth(4)

clusterPhi = ROOT.TH1F("clusterPhi", "ClusterPhi", 100, -3, 3)
clusterPhi.SetLineColor(kRed)
clusterPhi.SetLineWidth(4)

nHTTRT_data = ROOT.TH1F("nHTTRT_Data", "nHTTRT_Data", 100, 0, 100)
nHTTRT_data.SetLineColor(kRed)
nHTTRT_data.SetLineWidth(4)

nTRT_data = ROOT.TH1F("nTRT_Data", "nTRT_Data",  100, 0, 200)
nTRT_data.SetLineColor(kRed)
nTRT_data.SetLineWidth(4)

frac_data = ROOT.TH1F("frac_Data", "frac_Data", 101, -0.005, 1.005)
frac_data.SetLineColor(kRed)
frac_data.SetLineWidth(4)

deltaR = ROOT.TH1F("deltaR", "deltaR", 100, 0,0.1)
deltaR.SetLineColor(kRed)
deltaR.SetLineWidth(4)

clustetasize = ROOT.TH1F("clustetasize", "clustetasize", 100, 0, .1)
clustetasize.SetLineColor(kRed)
clustetasize.SetLineWidth(4)

clustphisize = ROOT.TH1F("clustphisize", "clustphisize", 100, 0, .1)
clustphisize.SetLineColor(kRed)
clustphisize.SetLineWidth(4)

deltaRvsfraction = ROOT.TH2F("deltaRvsFrac", "deltaRvsFrac", 100, 0, 0.1, 101, -0.005, 1.005)
deltaRvsfraction.SetLineColor(kRed)
deltaRvsfraction.SetLineWidth(4)

clustEmaxinEM1 = ROOT.TH1F("clustEmaxinEM1", "clustEmaxinEM1", 100, 0, 100)
clustEmaxinEM1.SetLineColor(kRed)
clustEmaxinEM1.SetLineWidth(4)

clustEmaxinEM2 = ROOT.TH1F("clustEmaxinEM2", "clustEmaxinEM2", 100, 0, 100)
clustEmaxinEM2.SetLineColor(kRed)
clustEmaxinEM2.SetLineWidth(4)

clustEtaSizeB0  = ROOT.TH1F("clustetasizeinpre", "clustetasizeinpre", 100,0, 0.01)
clustEtaSizeB0.SetLineColor(kRed)
clustEtaSizeB0.SetLineWidth(4)

clustTime_data = ROOT.TH1F("clusterTiming", "clusterTiming", 100, -10, 10)
clustTime_data.SetLineColor(kRed)
clustTime_data.SetLineWidth(4)

PreEnergy = ROOT.TH1F("EnergyinPre", "EnergyinPre", 100, -0.5, 99.5)
PreEnergy.SetLineColor(kRed)
PreEnergy.SetLineWidth(4)

EM1Energy = ROOT.TH1F("EnergyinEM1", "EnergyinEM1", 100, -0.5, 99.5)
EM1Energy.SetLineColor(kRed)
EM1Energy.SetLineWidth(4)

EM2Energy = ROOT.TH1F("EnergyinEM2", "EnergyinEM2", 100, -0.5, 99.5)
EM2Energy.SetLineColor(kRed)
EM2Energy.SetLineWidth(4)

BEMEnergy = ROOT.TH1F("BEMEnergy", "BEMEnergy", 100, -0.5, 99.5)
BEMEnergy.SetLineColor(kRed)
BEMEnergy.SetLineWidth(4)

MECinEM1 = ROOT.TH1F("MECinEM1", "MECinEM1", 100, -0.5, 99.5)
MECinEM1.SetLineColor(kRed)
MECinEM1.SetLineWidth(4)

MEC2inEM1 = ROOT.TH1F("Sec&FirstMECinEM1", "Sec&FirstMECinEM1", 100, -0.5, 99.5)
MEC2inEM1.SetLineColor(kRed)
MEC2inEM1.SetLineWidth(4)

FracEinMECinEM1 = ROOT.TH1F("FracofEinMECinEM1", "FracofEinMECinEM1", 101, -0.005, 1.005)
FracEinMECinEM1.SetLineColor(kRed)
FracEinMECinEM1.SetLineWidth(4)

FracEinMEC2inEM1 = ROOT.TH1F("FracofEinSec&FirstMECinEM1", "FracofEinSec&FirstMECinEM1", 101, -0.005, 1.005)
FracEinMEC2inEM1.SetLineColor(kRed)
FracEinMEC2inEM1.SetLineWidth(4)

nClusters = ROOT.TH1F("Numberofclusters", "Numberofclusters", 7,-0.5, 6.5)
nClusters.SetLineColor(kRed)
nClusters.SetLineWidth(4)

MEClustEnergy = ROOT.TH1F("MostEnergeticCluster", "MostEnergeticCluster", 100, 0, 100)
MEClustEnergy.SetLineColor(kRed)
MEClustEnergy.SetLineWidth(4)

MEClustEta = ROOT.TH1F("MEClustEta", "MEClustEta", 300, -3.01, 2.99)
MEClustEta.SetLineColor(kRed)
MEClustEta.SetLineWidth(4)

MEClustPhi = ROOT.TH1F("MEClustPhi", "MEClustPhi", 300, -3.01, 2.99)
MEClustPhi.SetLineColor(kRed)
MEClustPhi.SetLineWidth(4)

Eprofile = ROOT.TH2F("Eofclust1vsEofothers", "Eofclust1vsEofothers", 100, 0, 500, 100, 0, 100)
Eprofile.GetXaxis().SetTitle("Energy of first cluster")
Eprofile.GetYaxis().SetTitle("Energy of secondary clusters")
Eprofile.SetLineColor(kRed)
Eprofile.SetLineWidth(4)


Eprofile_max = ROOT.TH2F("Eprofile_max", "Eprofile_max", 100, 0, 500, 100, 0, 100)
Eprofile_max.GetXaxis().SetTitle("Energy of most energetic cluster")
Eprofile_max.GetYaxis().SetTitle("Energy of secondary clusters")
Eprofile_max.SetLineColor(kRed)
Eprofile_max.SetLineWidth(4)

eta_difference = ROOT.TH1F("eta_closenes", "eta_closeness", 100, -1.01, 0.99)
eta_difference.GetXaxis().SetTitle("Eta difference between two cluster candidates") 
eta_difference.GetYaxis().SetTitle("Number of candidates")
eta_difference.SetLineColor(kRed)
eta_difference.SetLineWidth(4)

OtherEnergyclust = ROOT.TH1F("otherenergyclust", "otherenergyclut", 100, 0, 500)
OtherEnergyclust.SetLineColor(kRed)
OtherEnergyclust.SetLineWidth(4)

eta_diff_max = ROOT.TH1F("eta_diff_max", "eta_diff_max", 100, -1.01, 0.99)
eta_diff_max.SetLineColor(kRed)
eta_diff_max.SetLineWidth(4)


nHTTRTCl_data = ROOT.TH1F("nHTTRTCl_Data", "nHTTRTCl_Data", 100, 0, 500)
nHTTRTCl_data.SetLineColor(kRed)
nHTTRTCl_data.SetLineWidth(4)

nTRTCl_data = ROOT.TH1F("nTRTCl_Data", "nTRTCl_Data",  100, 0, 1200)
nTRTCl_data.SetLineColor(kRed)
nTRTCl_data.SetLineWidth(4)

frac_dataCl = ROOT.TH1F("frac_DataCl", "frac_DataCl", 100, 0, 1)
frac_dataCl.SetLineColor(kRed)
frac_dataCl.SetLineWidth(4)

trig_nTRT = ROOT.TH1F("Trig_nTRT", "Trig_nTRT", 100, 0, 200)
trig_nTRT.SetLineColor(kRed)
trig_nTRT.SetLineWidth(4)

trig_nHTTRT = ROOT.TH1F("Trig_nHTTRT", "Trig_nHTTRT", 100,0, 200)
trig_nHTTRT.SetLineColor(kRed)
trig_nHTTRT.SetLineWidth(4)

trig_fraction = ROOT.TH1F("Trig_fraction", "Trig_fraction", 101, -0.005, 1.005)
trig_fraction.SetLineColor(kRed)
trig_fraction.SetLineWidth(4)

MEC = ROOT.TH1F("MEC", "MEC", 100, -0.5, 99.5)
MEC.SetLineColor(kRed)
MEC.SetLineWidth(4)

MEC2 = ROOT.TH1F("SecMEC", "SecMEC", 100, -0.5, 99.5)
MEC2.SetLineColor(kRed)
MEC2.SetLineWidth(4)

frac_vs_eta = ROOT.TH2F("frac_vs_eta", "frac_vs_eta", 101, -0.005, 1.005, 100, -3, 3)
frac_vs_eta.SetMarkerStyle(20)
frac_vs_eta.SetMarkerColor(kRed)

frac_vs_phi = ROOT.TH2F("frac_vs_phi", "frac_vs_phi", 101, -0.005, 1.005, 100, -3, 3)
frac_vs_phi.SetMarkerStyle(20)
frac_vs_phi.SetMarkerColor(kRed)

MEC2_and_MEC1 = ROOT.TH1F("MEC2+MEC1", "MEC2+MEC1",100, -0.5, 99.5)
MEC2_and_MEC1.SetLineWidth(4)
MEC2_and_MEC1.SetLineColor(kRed)

Sum_MEC_123_inEM1 = ROOT.TH1F("Sum_MEC123", "Sum_MEC123", 100, 0.05, 99.5)
Sum_MEC_123_inEM1.SetLineWidth(4)
Sum_MEC_123_inEM1.SetLineColor(kRed)

Frac_MEC_123_inEM1 = ROOT.TH1F("Frac_MEC123", "Frac_MEC123", 101, -0.005, 1.005)
Frac_MEC_123_inEM1.SetLineColor(kRed)
Frac_MEC_123_inEM1.SetLineWidth(4)

#data hitogram filling
for data_ev in range(0, chain_data.GetEntries()):
    tmp = chain_data.GetEntry(data_ev)

    if chain_data.nHTTRT2 > 0 and chain_data.EF_g_nocut_hiptrtL2[0] == 1 :
        nClusters.Fill(chain_data.nClust[0])
        for hits in range(0, chain_data.nClust[0]):
            if abs(chain_data.clustEta[hits]) < 3:
                nHTTRT_data.Fill(chain_data.nHTTRT2[hits])
                nTRT_data.Fill(chain_data.nTRT2[hits])
                nTRTCl_data.Fill(chain_data.nTRT[hits])
                nHTTRTCl_data.Fill(chain_data.nHTTRT[hits])
                if chain_data.nTRT2[hits] > 0:
                    frac_data.Fill(chain_data.nHTTRT2[hits]/float(chain_data.nTRT2[hits]))
                    frac_vs_eta.Fill(chain_data.nHTTRT2[hits]/float(chain_data.nTRT2[hits]), chain_data.clustEta[hits])
                    frac_vs_phi.Fill(chain_data.nHTTRT2[hits]/float(chain_data.nTRT2[hits]), chain_data.clustPhi[hits])
                if chain_data.nTRT[hits] > 0:
                    frac_dataCl.Fill(chain_data.nHTTRT[hits]/float(chain_data.nTRT[hits]))
                if chain_data.clustEta237[hits] > 0:
                    size_eta = chain_data.clustEta237[hits]
                    size_phi = chain_data.clustPhi237[hits]
                    clustetasize.Fill(size_eta)
                    clustphisize.Fill(size_phi)
                    delta_r = math.sqrt(size_eta*size_eta + size_phi*size_phi)
                    deltaR.Fill(delta_r)
                    if chain_data.nTRT2[hits] > 0:
                        deltaRvsfraction.Fill(delta_r, chain_data.nHTTRT2[hits]/float(chain_data.nTRT2[hits]))

                clustEmaxinEM1.Fill(chain_data.clustEmaxBar1[hits]/1000.)
                clustEmaxinEM2.Fill(chain_data.clustEmaxBar2[hits]/1000.)
            #print chain_data.clustEmaxBar1[hits]
                clustEtaSizeB0.Fill(chain_data.clustEtaSizeB0[hits])
                clusterE.Fill(chain_data.clustE[hits]/1000.)
                clusterEt.Fill(chain_data.clustEt[hits]/1000.)
                clusterEta.Fill(chain_data.clustEta[hits])
                clusterPhi.Fill(chain_data.clustPhi[hits])
                clustTime_data.Fill(chain_data.clustTime[hits])
                PreEnergy.Fill(chain_data.TotalEnergyinPre[hits]/1000.)
                EM1Energy.Fill(chain_data.TotalEnergyinEM1[hits]/1000.)
                EM2Energy.Fill(chain_data.TotalEnergyinEM2[hits]/1000.)
                BEMEnergy.Fill(chain_data.TotalEnergyBeyondEM[hits]/1000.)
                MECinEM1.Fill(chain_data.hotCellinEM1Energy[hits][0]/1000.)
                MEC2inEM1.Fill((chain_data.hotCellinEM1Energy[hits][1]+chain_data.hotCellinEM1Energy[hits][0])/1000.)
                MEC.Fill(chain_data.hotCellEnergy[hits][0]/1000.)
                MEC2.Fill(chain_data.hotCellEnergy[hits][1]/1000.)
                MEC2_and_MEC1.Fill((chain_data.hotCellEnergy[hits][0]+chain_data.hotCellEnergy[hits][1])/1000.)
                Sum_MEC_123_inEM1.Fill((chain_data.hotCellinEM1Energy[hits][0]+chain_data.hotCellinEM1Energy[hits][1]+chain_data.hotCellinEM1Energy[hits][2])/1000.)

                if chain_data.TotalEnergyinEM1[hits] > 0:
                    FracEinMECinEM1.Fill(chain_data.hotCellinEM1Energy[hits][0]/chain_data.TotalEnergyinEM1[hits])
                    FracEinMEC2inEM1.Fill((chain_data.hotCellinEM1Energy[hits][1]+chain_data.hotCellinEM1Energy[hits][0])/chain_data.TotalEnergyinEM1[hits])
                    Frac_MEC_123_inEM1.Fill((chain_data.hotCellinEM1Energy[hits][0]+chain_data.hotCellinEM1Energy[hits][1]+chain_data.hotCellinEM1Energy[hits][2])/chain_data.TotalEnergyinEM1[hits])
                
                for trig_ev in range(0, len(chain_data.HT_hits_in_phi0)):
                    trig_nTRT.Fill(chain_data.Total_number_of_hits_in_phi0[trig_ev])
                    trig_nHTTRT.Fill(chain_data.HT_hits_in_phi0[trig_ev])
                    trig_fraction.Fill(chain_data.HT_hits_in_phi0[trig_ev]/float(chain_data.Total_number_of_hits_in_phi0[trig_ev]))
                    
##The following was to make histograms checking for the most energetic cluster in an multiple-cluster event
                # max_E= 0
                # index = 0
                # if chain_data.nClust[0] > 1 and chain_data.EF_g_nocut_hiptrtL2[0] == 1:
                #     for clust in range(0, chain_data.nClust[0]):

                #         if max_E < chain_data.clustE[clust]:
                #             max_E = chain_data.clustE[clust]
                #             index = clust
                #     MEClustEnergy.Fill(chain_data.clustE[index]/1000.)
                #     MEClustEta.Fill(chain_data.clustEta[index])
                #     MEClustPhi.Fill(chain_data.clustPhi[index])
                #     for clustnew in range(0, chain_data.nClust[0]):
                #         if clustnew ==0:
                #             for clustnew2 in range(1, chain_data.nClust[0]):
                #                 Eprofile.Fill(chain_data.clustE[clustnew]/1000., chain_data.clustE[clustnew2]/1000.)
                #         if clustnew == index :
                #             continue
                #         eta_diff_max.Fill(chain_data.clustEta[index]-chain_data.clustEta[clustnew])
                #         OtherEnergyclust.Fill(chain_data.clustE[clustnew]/1000.)
                #         Eprofile_max.Fill(chain_data.clustE[index]/1000. , chain_data.clustE[clustnew]/1000.)




#for MC histogram filling 
clusterE_mc = ROOT.TH1F("clustE_mc", "clustE_mc", 100,0 ,150)
clusterE_mc.SetLineWidth(4)
clusterE_mc.SetLineColor(kBlue)

clusterEt_mc = ROOT.TH1F("clustEt_mc", "clustEt_mc" ,100,0,150)
clusterEt_mc.SetLineColor(kBlue)
clusterEt_mc.SetLineWidth(4)

clusterEta_mc = ROOT.TH1F("clusterEta_mc", "clusterEta_mc", 100, -3 ,3)
clusterEta_mc.SetLineColor(kBlue)
clusterEta_mc.SetLineWidth(4)

clusterPhi_mc = ROOT.TH1F("clusterPhi_mc", "clusterPhi_mc", 100, -3 ,3)
clusterPhi_mc.SetLineColor(kBlue)
clusterPhi_mc.SetLineWidth(4)

nHTTRT_mc = ROOT.TH1F("nHTTRT_MC", "nHTTRT_MC", 100, 0, 100)
nHTTRT_mc.SetLineColor(kBlue)
nHTTRT_mc.SetLineWidth(4)

nTRT_mc = ROOT.TH1F("nTRT_MC", "nTRT_MC", 100, 0, 200)
nTRT_mc.SetLineColor(kBlue)
nTRT_mc.SetLineWidth(4)

frac_mc = ROOT.TH1F("frac_MC", "frac_MC", 101, -0.005, 1.005)
frac_mc.SetLineColor(kBlue)
frac_mc.SetLineWidth(4)

delta_R_mc = ROOT.TH1F("delta_R_mc", "delta_R_mc", 100, 0, 0.1)
delta_R_mc.SetLineColor(kBlue)
delta_R_mc.SetLineWidth(4)

deltaR_vs_frac_mc = ROOT.TH2F("deltaRvsfrcmc", "deltaRvsfrcmc", 100,0,0.1, 01, -0.005, 1.005)
deltaR_vs_frac_mc.SetLineColor(kBlue)
deltaR_vs_frac_mc.SetLineWidth(4)

clustTime_mc = ROOT.TH1F("clusterTiming_mc", "clusterTiming_mc", 100, -10, 10)
clustTime_mc.SetLineColor(kBlue)
clustTime_mc.SetLineWidth(4)

PreEnergy_mc = ROOT.TH1F("EnergyinPre_mc", "EnergyinPre_mc", 100, -0.5, 99.5)
PreEnergy_mc.SetLineColor(kBlue)
PreEnergy_mc.SetLineWidth(4)

EM1Energy_mc = ROOT.TH1F("EM1Energy_mc", "EM1Energy_mc", 100, -0.5, 99.5)
EM1Energy_mc.SetLineColor(kBlue)
EM1Energy_mc.SetLineWidth(4)

EM2Energy_mc = ROOT.TH1F("EM2Energy_mc", "EM2Energy_mc", 100, -0.5, 99.5)
EM2Energy_mc.SetLineColor(kBlue)
EM2Energy_mc.SetLineWidth(4)

BEMEnergy_mc = ROOT.TH1F("BEMEnergy_mc", "BEMEnergy_mc", 100, -0.5, 99.5)
BEMEnergy_mc.SetLineColor(kBlue)
BEMEnergy_mc.SetLineWidth(4)

MECinEM1_mc = ROOT.TH1F("MECinEM1_mc", "MECinEM1_mc", 100, -0.5, 99.5)
MECinEM1_mc.SetLineColor(kBlue)
MECinEM1_mc.SetLineWidth(4)

MEC2inEM1_mc = ROOT.TH1F("Sec&FirstMECinEM1_mc", "Sec&FirstMECinEM1_mc", 100, -0.5, 99.5)
MEC2inEM1_mc.SetLineColor(kBlue)
MEC2inEM1_mc.SetLineWidth(4)

FracEinMECinEM1_mc = ROOT.TH1F("FracofEinMECinEM1_mc", "FracofEinMECinEM1_mc", 101, -0.005, 1.005)
FracEinMECinEM1_mc.SetLineColor(kBlue)
FracEinMECinEM1_mc.SetLineWidth(4)

FracEinMEC2inEM1_mc = ROOT.TH1F("FracofEinSec&FirstMECinEM1_mc", "FracofEinSec&FirstMECinEM1_mc", 101, -0.005, 1.005)
FracEinMEC2inEM1_mc.SetLineColor(kBlue)
FracEinMEC2inEM1_mc.SetLineWidth(4)

nClusters_mc = ROOT.TH1F("numberofclustersMC", "numberofclustersMC", 7, -0.5, 6.5)
nClusters_mc.SetLineColor(kBlue)
nClusters_mc.SetLineWidth(4)

MEClustEnergy_mc = ROOT.TH1F("MostEnergeticCluster_mc", "MostEnergeticCluster_mc", 100, 0, 100)
MEClustEnergy_mc.SetLineColor(kBlue)
MEClustEnergy_mc.SetLineWidth(4)

MEClustEta_mc = ROOT.TH1F("MEClustEta_mc", "MEClustEta_mc", 300, -3.01, 2.99)
MEClustEta_mc.SetLineColor(kBlue)
MEClustEta_mc.SetLineWidth(4)

MEClustPhi_mc = ROOT.TH1F("MEClustPhi_mc", "MEClustPhi_mc", 300, -3.01, 2.99)
MEClustPhi_mc.SetLineColor(kBlue)
MEClustPhi_mc.SetLineWidth(4)

Eprofile_mc = ROOT.TH2F("Eofclust1vsEofothers_mc", "Eofclust1vsEofothers_mc", 100, 0, 500, 100, 0, 100)
Eprofile_mc.GetXaxis().SetTitle("Energy of first cluster")
Eprofile_mc.GetYaxis().SetTitle("Energy of secondary clusters")
Eprofile_mc.SetLineColor(kBlue)
Eprofile_mc.SetLineWidth(4)


Eprofile_max_mc = ROOT.TH2F("Eprofile_max_mc", "Eprofile_max_mc", 100, 0, 500, 100, 0, 100)
Eprofile_max_mc.GetXaxis().SetTitle("Energy of most energetic cluster")
Eprofile_max_mc.GetYaxis().SetTitle("Energy of secondary clusters")
Eprofile_max_mc.SetLineColor(kBlue)
Eprofile_max_mc.SetLineWidth(4)

# eta_difference = ROOT.TH1F("eta_closenes", "eta_closeness", 100, -1.01, 0.99)
# eta_difference.GetXaxis().SetTitle("Eta difference between two cluster candidates") 
# eta_difference.GetYaxis().SetTitle("Number of candidates")
# eta_difference.SetLineColor(kBlue)
# eta_difference.SetLineWidth(4)

OtherEnergyclust_mc = ROOT.TH1F("otherenergyclust_mc", "otherenergyclut_mc", 100, 0, 500)
OtherEnergyclust_mc.SetLineColor(kBlue)
OtherEnergyclust_mc.SetLineWidth(4)

eta_diff_max_mc = ROOT.TH1F("eta_diff_max_mc", "eta_diff_max_mc", 100, -1.01, 0.99)
eta_diff_max_mc.SetLineColor(kBlue)
eta_diff_max_mc.SetLineWidth(4)

nHTTRTCl_mc = ROOT.TH1F("nHTTRTCl_MC", "nHTTRTCl_MC", 100, 0, 500)
nHTTRTCl_mc.SetLineColor(kBlue)
nHTTRTCl_mc.SetLineWidth(4)

nTRTCl_mc = ROOT.TH1F("nTRTCl_MC", "nTRTCl_MC", 100, 0, 1200)
nTRTCl_mc.SetLineColor(kBlue)
nTRTCl_mc.SetLineWidth(4)

frac_mcCl = ROOT.TH1F("frac_MCCl", "frac_MCCl", 101, -0.005, 1.005)
frac_mcCl.SetLineColor(kBlue)
frac_mcCl.SetLineWidth(4)

trig_nTRT_mc = ROOT.TH1F("Trig_nTRT_mc", "Trig_nTRT_mc", 100, 0, 200)
trig_nTRT_mc.SetLineColor(kBlue)
trig_nTRT_mc.SetLineWidth(4)

trig_nHTTRT_mc = ROOT.TH1F("Trig_nHTTR_mcT", "Trig_nHTTRT_mc", 100,0, 200)
trig_nHTTRT_mc.SetLineColor(kBlue)
trig_nHTTRT_mc.SetLineWidth(4)

trig_fraction_mc = ROOT.TH1F("Trig_fraction_mc", "Trig_fraction_mc", 100, -0.005, 1.005)
trig_fraction_mc.SetLineColor(kBlue)
trig_fraction_mc.SetLineWidth(4)

MEC_mc = ROOT.TH1F("MEC_mc", "MEC_mc", 100, -0.5, 99.5)
MEC_mc.SetLineColor(kBlue)
MEC_mc.SetLineWidth(4)

MEC2_mc = ROOT.TH1F("SecMEC_mc", "SecMEC_mc", 100, -0.5, 99.5)
MEC2_mc.SetLineColor(kBlue)
MEC2_mc.SetLineWidth(4)

clustetasize_mc = ROOT.TH1F("clusteta237", "clusteta237", 100, 0, .1)
clustetasize_mc.SetLineColor(kBlue)
clustetasize_mc.SetLineWidth(4)

clustphisize_mc = ROOT.TH1F("clustphi237_mc", "clustphi237_mc", 100, 0, .1)
clustphisize_mc.SetLineColor(kBlue)
clustphisize_mc.SetLineWidth(4)

frac_vs_eta_mc = ROOT.TH2F("frac_vs_eta_mc", "frac_vs_eta_mc", 101, -0.005, 1.005, 100, -3, 3)
frac_vs_eta_mc.SetMarkerStyle(20)
frac_vs_eta_mc.SetMarkerColor(kBlue)

frac_vs_phi_mc = ROOT.TH2F("frac_vs_phi_mc", "frac_vs_phi_mc", 101, -0.005, 1.005, 100, -3, 3)
frac_vs_phi_mc.SetMarkerStyle(20)
frac_vs_phi_mc.SetMarkerColor(kBlue)

MEC2_and_MEC1_mc = ROOT.TH1F("MEC2+MEC1_mc", "MEC2+MEC1_mc", 100, -0.5, 99.5)
MEC2_and_MEC1_mc.SetLineColor(kBlue)
MEC2_and_MEC1_mc.SetLineWidth(4)

Sum_MEC_123_inEM1_mc= ROOT.TH1F("Sum_MEC123_mc", "Sum_MEC123_mc", 100, 0.05, 99.5)
Sum_MEC_123_inEM1_mc.SetLineWidth(4)
Sum_MEC_123_inEM1_mc.SetLineColor(kBlue)

Frac_MEC_123_inEM1_mc = ROOT.TH1F("Frac_MEC123_mc", "Frac_MEC123_mc", 101, -0.005, 1.005)
Frac_MEC_123_inEM1_mc.SetLineColor(kBlue)
Frac_MEC_123_inEM1_mc.SetLineWidth(4)

w_variable = []
for number in range(0, 5):
    hist = ROOT.TH1F("w-variable for "+str(number)+" cells ", "w-variable for "+str(number)+ " cells", 101, -0.005, 1.005)
    w_variable.append(hist)
    

for mc_ev in range(0, chain_mc.GetEntries()):
    tmpmc = chain_mc.GetEntry(mc_ev)
    if chain_mc.nHTTRT2 > 0 and chain_mc.EF_g_nocut_hiptrtL2[0] == 1:
        nClusters_mc.Fill(chain_mc.nClust[0])
        for hits_mc in range(0, chain_mc.nClust[0]):
            if abs(chain_mc.clustEta[hits_mc]) < 3:
                nHTTRT_mc.Fill(chain_mc.nHTTRT2[hits_mc])
                nTRT_mc.Fill(chain_mc.nTRT2[hits_mc])
                nTRTCl_mc.Fill(chain_mc.nTRT[hits_mc])
                nHTTRTCl_mc.Fill(chain_mc.nHTTRT[hits_mc])
                if chain_mc.nTRT2[hits_mc] > 0:
                    frac_mc.Fill(chain_mc.nHTTRT2[hits_mc]/float(chain_mc.nTRT2[hits_mc]))
                    frac_vs_eta_mc.Fill(chain_mc.nHTTRT2[hits_mc]/float(chain_mc.nTRT2[hits_mc]), chain_mc.clustEta[hits_mc])
                    frac_vs_phi_mc.Fill(chain_mc.nHTTRT2[hits_mc]/float(chain_mc.nTRT2[hits_mc]), chain_mc.clustPhi[hits_mc])
                if chain_mc.nTRT[hits_mc] > 0:
                    frac_mcCl.Fill(chain_mc.nHTTRT[hits_mc]/float(chain_mc.nTRT[hits_mc]))
                if chain_mc.clustEta237[hits_mc] > 0:
                    size_eta_mc = chain_mc.clustEta237[hits_mc]
                    size_phi_mc = chain_mc.clustPhi237[hits_mc]
                    clustetasize_mc.Fill(size_eta_mc)
                    clustphisize_mc.Fill(size_phi_mc)
                    delta_r_mc = math.sqrt(size_eta_mc*size_eta_mc + size_phi_mc*size_phi_mc)
                    delta_R_mc.Fill(delta_r_mc)
                    if chain_mc.nTRT2[hits_mc] > 0:
                        deltaR_vs_frac_mc.Fill(delta_r_mc, chain_mc.nHTTRT2[hits_mc]/float(chain_mc.nTRT2[hits_mc]))
                clusterEt_mc.Fill(chain_mc.clustEt[hits_mc]/1000.)
                clusterE_mc.Fill(chain_mc.clustE[hits_mc]/1000.)
                clusterEta_mc.Fill(chain_mc.clustEta[hits_mc])
                clusterPhi_mc.Fill(chain_mc.clustPhi[hits_mc])
                clustTime_mc.Fill(chain_mc.clustTime[hits_mc])
                PreEnergy_mc.Fill(chain_mc.TotalEnergyinPre[hits_mc]/1000.)
                EM1Energy_mc.Fill(chain_mc.TotalEnergyinEM1[hits_mc]/1000.)
                EM2Energy_mc.Fill(chain_mc.TotalEnergyinEM2[hits_mc]/1000.)
                BEMEnergy_mc.Fill(chain_mc.TotalEnergyBeyondEM[hits_mc]/1000.)
                MECinEM1_mc.Fill(chain_mc.hotCellinEM1Energy[hits_mc][0]/1000.)
                MEC2inEM1_mc.Fill((chain_mc.hotCellinEM1Energy[hits_mc][1]+chain_mc.hotCellinEM1Energy[hits_mc][0])/1000.)
                MEC_mc.Fill(chain_mc.hotCellEnergy[hits_mc][0]/1000.)
                MEC2_mc.Fill(chain_mc.hotCellEnergy[hits_mc][1]/1000.)
                MEC2_and_MEC1_mc.Fill((chain_mc.hotCellEnergy[hits_mc][0]+chain_mc.hotCellEnergy[hits_mc][1])/1000.)
                Sum_MEC_123_inEM1_mc.Fill((chain_mc.hotCellinEM1Energy[hits_mc][0]+chain_mc.hotCellinEM1Energy[hits_mc][1]+chain_mc.hotCellinEM1Energy[hits_mc][2])/1000.)

                if chain_mc.TotalEnergyinEM1[hits_mc] > 0:
                    FracEinMECinEM1_mc.Fill(chain_mc.hotCellinEM1Energy[hits_mc][0]/float(chain_mc.TotalEnergyinEM1[hits_mc]))
                    FracEinMEC2inEM1_mc.Fill((chain_mc.hotCellinEM1Energy[hits_mc][1]+chain_mc.hotCellinEM1Energy[hits_mc][0])/float(chain_mc.TotalEnergyinEM1[hits_mc]))
                    Frac_MEC_123_inEM1_mc.Fill((chain_mc.hotCellinEM1Energy[hits_mc][0]+chain_mc.hotCellinEM1Energy[hits_mc][1]+chain_mc.hotCellinEM1Energy[hits_mc][2])/chain_mc.TotalEnergyinEM1[hits_mc])

                for trig_ev_mc in range(0, len(chain_mc.HT_hits_in_phi0)):
                    trig_nTRT_mc.Fill(chain_mc.Total_number_of_hits_in_phi0[trig_ev_mc])
                    trig_nHTTRT_mc.Fill(chain_mc.HT_hits_in_phi0[trig_ev_mc])
                    trig_fraction_mc.Fill(chain_mc.HT_hits_in_phi0[trig_ev_mc]/float(chain_mc.Total_number_of_hits_in_phi0[trig_ev_mc]))


                # max_E_mc= 0
                # index_mc = 0
                # if chain_mc.nClust[0] > 1 and chain_mc.EF_g_nocut_hiptrtL2[0] == 1:
                #     for clustmc in range(0, chain_mc.nClust[0]):
                #         if max_E_mc < chain_mc.clustE[clustmc]:
                #             max_E_mc = chain_mc.clustE[clustmc]
                #             index_mc = clustmc
                #     MEClustEnergy_mc.Fill(chain_mc.clustE[index_mc]/1000.)
                #     MEClustEta_mc.Fill(chain_mc.clustEta[index_mc])
                #     MEClustPhi_mc.Fill(chain_mc.clustPhi[index_mc])
                #     for clustnewmc in range(0, chain_mc.nClust[0]):
                #         if clustnewmc ==0:
                #             for clustnewmc2 in range(1, chain_mc.nClust[0]):
                #                 Eprofile_mc.Fill(chain_mc.clustE[clustnewmc]/1000., chain_mc.clustE[clustnewmc2]/1000.)
                #                 if clustnewmc == index_mc :
                #                     continue
                #                 eta_diff_max_mc.Fill(chain_mc.clustEta[index_mc]-chain_mc.clustEta[clustnewmc])
                #                 OtherEnergyclust_mc.Fill(chain_mc.clustE[clustnewmc]/1000.)
                #                 Eprofile_max_mc.Fill(chain_mc.clustE[index_mc]/1000. , chain_mc.clustE[clustnewmc]/1000.)


## ready plots for comparisons

Legend = ROOT.TLegend(0.63, 0.8,0.9,0.9)
Legend.AddEntry(clusterE, args.DATALABEL, "L")
Legend.AddEntry(clusterE_mc, args.MCLABEL, "L")

canvas_Text = ROOT.TPaveText(0.5,0.6,0.9,0.7, "NDC")

CanvasMaker(clusterE, clusterE_mc, Legend, norm=args.NORM, save=args.SAVE, Canvasname="energy_canvas", histogram_Title="Cluster_Energy", Text='ClusterEnergy')
CanvasMaker(clusterEt, clusterEt_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="Et_canvas", histogram_Title="Transverse_Energy", Text='ClusterEt')
CanvasMaker(nHTTRT_data, nHTTRT_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="nHTTRT_canvas", histogram_Title="nHTTRT2", Text="Number of HT TRT hits in a road within cone")
CanvasMaker(nTRT_data, nTRT_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="nTRT_canvas", histogram_Title="nTRT2", Text="Total number of TRT hits in road within cone")
CanvasMaker(frac_data, frac_mc,  Legend, norm=args.NORM, save=args.SAVE,Canvasname="frac_canavs", histogram_Title="Fraction_of_HT_TRT_hits_in_Road_within_cone")
CanvasMaker(deltaR, delta_R_mc,   Legend, norm=args.NORM, save=False, Canvasname="deltaR_canvas", histogram_Title="DeltaR_Variable")
CanvasMaker(PreEnergy, PreEnergy_mc,  Legend,  norm=args.NORM, save=args.SAVE, Canvasname="PreEnergy_Canvas", histogram_Title="TotalEnergyinPre", Text="Energy deposited in Presampler")
CanvasMaker(EM1Energy, EM1Energy_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="EnergyinEM1_Canvas", histogram_Title="TotalEnergyinEM1", Text="Energy deposited in first layer of Calorimeter")
CanvasMaker(EM2Energy, EM2Energy_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="EnergyinEM2_canvas", histogram_Title="TotalEnergyinEM2", Text="Energy deposited in second layer of Calorimeter")
CanvasMaker(BEMEnergy, BEMEnergy_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="BEMEnergy_canvas", histogram_Title="TotalEnergyBeyondEM", Text="Energy desposition beyond the EM Calorimeter")
CanvasMaker(MECinEM1, MECinEM1_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="MECinEM1_canvas", histogram_Title="hotCellinEM1Energy[clustnum][0]", Text="Energy of the most energetic cell over in EM1")
CanvasMaker(MEC2inEM1, MEC2inEM1_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="MEC2inEM1_canvas", histogram_Title="hotCellinEM1Energy[clustnum][0]+hotCellinEM1Energy[clustnum][1]", Text="Sum of Energy of two most energetic cells in EM1")
CanvasMaker(FracEinMECinEM1, FracEinMECinEM1_mc,  Legend,norm=args.NORM, save = False, Canvasname="FracEinMECinEM1_canvas", histogram_Title="", Text="Fraction of Energy in most energetic cell over total energy in EM1")
CanvasMaker(FracEinMEC2inEM1, FracEinMEC2inEM1_mc, Legend, norm=args.NORM, save = False, Canvasname="FracEinMEC2inEM1_canvas", histogram_Title="", Text="Fraction of Energy in two most energetic cell over total energy in EM1")
CanvasMaker(nClusters, nClusters_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="nClusters_canvas", histogram_Title="nClust", Text="Number of clusters in event")
CanvasMaker(clustTime_data, clustTime_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="clusttime_canvas", histogram_Title="clustTime", Text="cluster timing variable")
CanvasMaker(nHTTRTCl_data, nHTTRTCl_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="nHTTRTCl_canvas", histogram_Title="nHTTRT", Text="number of HT TRT hits from cone around cluster")
CanvasMaker(nTRTCl_data, nTRTCl_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="nTRTCl_canvas", histogram_Title="nTRT", Text="Total number of TRT hits from cone")
CanvasMaker(frac_dataCl, frac_mcCl,  Legend, norm=args.NORM, save=False, Canvasname="FracCl_canvas", histogram_Title="", Text="Fraction of TRT hits from cone around cluster (nHTTRT/nTRT)")
CanvasMaker(trig_nHTTRT, trig_nHTTRT_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="trig_nHTTRT_canvas", histogram_Title="HT_hits_in_phi0", Text="Number of HT TRT hits from HIP trigger")
CanvasMaker(trig_nTRT, trig_nTRT_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="trig_nTRT_canvas", histogram_Title="Total_number_of_hits_in_phi0", Text="Total number of TRT hits from HIP Trigger")
CanvasMaker(trig_fraction, trig_fraction_mc,  Legend,  norm=args.NORM, save=args.SAVE, Canvasname="trig_Fraction_Canvas", histogram_Title="trig_Fraction", Text="Fraction of TRT hits from HIP trigger")
CanvasMaker(MEC, MEC_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="MEC_canvas", histogram_Title="hotCellEnergy[clustnum][0]", Text="Energy of the most energetic cell over of the cluster")
CanvasMaker(MEC2_and_MEC1, MEC2_and_MEC1_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="MEC2_and_MEC1_canvas", histogram_Title="hotCellEnergy[clustnum][0]+hotCellEnergy[clustnum][1]", Text="Sum of Energy of two most energetic cells of the cluster")
CanvasMaker(MEC2, MEC2_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="MEC2_canvas", histogram_Title="hotCellEnergy[clustnum][1]", Text="Energy of the second most energetic cell over of the cluster")

CanvasMaker(clustetasize, clustetasize_mc,  Legend,norm=args.NORM,save=args.SAVE, Canvasname="clustetasize_canvas", histogram_Title="clustEta237", Text="Resolution of size in eta of cluster with a phi cut")
CanvasMaker(clustphisize, clustphisize_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="clustphisize_canvas", histogram_Title="clustPhi237", Text="Resolution of size in phi of cluster with a phi cut")
CanvasMaker(clusterEta, clusterEta_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="clustEta_canvas", histogram_Title="clustEta", Text="Eta of all clusters passing the trigger and preselection")
CanvasMaker(clusterPhi, clusterPhi_mc,  Legend, norm=args.NORM, save=args.SAVE, Canvasname="clustphi_canvas", histogram_Title="clustPhi", Text="Phi of all clusters passing the trigger and preselection")

CanvasMaker(Sum_MEC_123_inEM1, Sum_MEC_123_inEM1_mc, Legend, norm=args.NORM, save=args.SAVE, Canvasname="Sum_MEC_123_canvas", histogram_Title="hotCellinEM1Energy[clustnum][0] + hotCellinEM1Energy[clustnum][1] + hotCellinEM1Energy[clustnum][0]", Text="sum of three most energetic cells in EM1")
CanvasMaker(Frac_MEC_123_inEM1, Frac_MEC_123_inEM1_mc, Legend, norm=args.NORM, save=args.SAVE, Canvasname="Frac_MEC_123_canvas", histogram_Title="", Text="Fraction of energy deposited in three most energetic cells in EM1")

code.interact()
# ##Energy profile data
# Eprofile_data_canvas = ROOT.TCanvas()
# Eprofile_data_canvas.cd()
# Eprofile.SetMarkerStyle(20)
# Eprofile.Draw()
# #Legend.Draw()

# ##Energy MAX profile data
# Eprofile_max_data_canvas = ROOT.TCanvas()
# Eprofile_max_data_canvas.cd()
# Eprofile_max.SetMarkerStyle(20)
# Eprofile_max.Draw()
# #Legend.Draw()

# ##Energy profile mc
# Eprofile_mc_canvas = ROOT.TCanvas()
# Eprofile_mc_canvas.cd()
# Eprofile_mc.SetMarkerStyle(20)
# Eprofile_mc.Draw()
# #Legend.Draw()

# ##Energy MAX profile MC
# Eprofile_max_mc_canvas = ROOT.TCanvas()
# Eprofile_max_mc_canvas.cd()
# Eprofile_max_mc.SetMarkerStyle(20)
# Eprofile_max_mc.Draw()
# #Legend.Draw()
