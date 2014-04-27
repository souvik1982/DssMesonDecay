#include <iostream>
#include "TCanvas.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "THStack.h"
#include "TF1.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TFile.h"
#include "TPaveStats.h"
#include <set>
#include <map>

double pi=3.14159265358979;

bool DsMass_optimize=false;
bool mBC_optimize=false;
bool DeltaM_optimize=false;
bool diffD0_optimize=false;
bool dPhi_optimize=false;
bool pi0_optimize=false;

std::string decay="pieta";
// options: KKpi, KsK, pieta, pietaprime, KKpipi0, pipipi, KsKmpipi, pipi0eta, pietaprimerho
int decayNumber;

double dsPlusMCut_center, dsPlusMCut_range;
double mbcCut_center, mbcCut_range;
double deltaMCut_center, deltaMCut_range;
double diffD0Cut;
double dPhiCutLess;
double pi0MassCut_center=0.135, pi0MassCut_range=0;
double electronEnergyThreshold=0.15;
double branchingFr_mode;

double nSignalSample=9988*2;
double genericScale=1/19.2;
double continuScale=1/5.;

double luminosity=586; // /pb
double prodCrossSection_DsDss=948;
double branchingFr_signal=0.942*0.0065;

std::string signalFileName_Dsp="/nfs/cor/an2/souvik/MC_vtosll_Dsp_";
std::string signalFileName_Dsm="/nfs/cor/an2/souvik/MC_vtosll_Dsm_";

void setValues()
{
  if (decay=="KKpi")
  {
    decayNumber=401;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.011;
    mbcCut_center=2.112; mbcCut_range=0.005;
    deltaMCut_center=0.155; deltaMCut_range=0.009;
    diffD0Cut=-0.002;
    dPhiCutLess=0.06;
    branchingFr_mode=0.055;
  }
  else if (decay=="KsK")
  {
    decayNumber=400;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.012;
    mbcCut_center=2.112; mbcCut_range=0.006;
    deltaMCut_center=0.158; deltaMCut_range=0.01;
    diffD0Cut=-0.002;
    dPhiCutLess=0.09;
    branchingFr_mode=0.0149;
  }
  else if (decay=="pieta")
  {
    decayNumber=440;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.015;
    mbcCut_center=2.112; mbcCut_range=0.007;
    deltaMCut_center=0.158; deltaMCut_range=0.013;
    diffD0Cut=-0.007;
    dPhiCutLess=0.07;
    branchingFr_mode=0.0158*0.3931;
  }
  else if (decay=="pietaprime")
  {
    decayNumber=460;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.011;
    mbcCut_center=2.112; mbcCut_range=0.011;
    deltaMCut_center=0.158; deltaMCut_range=0.013;
    diffD0Cut=-0.003;
    dPhiCutLess=0.07;
    branchingFr_mode=0.038*0.446*0.3931;
  }
  else if (decay=="KKpipi0")
  {
    decayNumber=404;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.009;
    mbcCut_center=2.112; mbcCut_range=0.007;
    deltaMCut_center=0.158; deltaMCut_range=0.011;
    diffD0Cut=-0.002;
    dPhiCutLess=0.07;
    branchingFr_mode=0.056; // check
  }
  else if (decay=="pipipi")
  {
    decayNumber=421;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.013;
    mbcCut_center=2.112; mbcCut_range=0.005;
    deltaMCut_center=0.158; deltaMCut_range=0.009;
    diffD0Cut=-0.001;
    dPhiCutLess=0.06;
    branchingFr_mode=0.0111;
  }
  else if (decay=="KsKmpipi")
  {
    decayNumber=406;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.007;
    mbcCut_center=2.112; mbcCut_range=0.007;
    deltaMCut_center=0.158; deltaMCut_range=0.009;
    diffD0Cut=-0.004;
    dPhiCutLess=0.07;
    branchingFr_mode=0.0164;
  }
  else if (decay=="pipi0eta")
  {
    decayNumber=441;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.014;
    mbcCut_center=2.112; mbcCut_range=0.006;
    deltaMCut_center=0.158; deltaMCut_range=0.009;
    diffD0Cut=-0.003;
    dPhiCutLess=0.07;
    branchingFr_mode=0.131*0.995*0.3931;
  }
  else if (decay=="pietaprimerho")
  {
    decayNumber=480;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.018;
    mbcCut_center=2.112; mbcCut_range=0.004;
    deltaMCut_center=0.155; deltaMCut_range=0.008;
    diffD0Cut=-0.004;
    dPhiCutLess=0.09;
    branchingFr_mode=0.038*0.294;
  }
}

void setFileNames()
{
  if (decay=="KKpi")
  {
    signalFileName_Dsp=signalFileName_Dsp+"KKpi/DsTaggedDecaysProc_MC_vtosll_Dsp_KKpi.root";
    signalFileName_Dsm=signalFileName_Dsm+"KKpi/DsTaggedDecaysProc_MC_vtosll_Dsm_KKpi.root";
  }
  else if (decay=="KsK")
  {
    signalFileName_Dsp=signalFileName_Dsp+"KsK/DsTaggedDecaysProc_MC_vtosll_Dsp_KsK.root";
    signalFileName_Dsm=signalFileName_Dsm+"KsK/DsTaggedDecaysProc_MC_vtosll_Dsm_KsK.root";
  }
  else if (decay=="pieta")
  {
    signalFileName_Dsp=signalFileName_Dsp+"pieta/DsTaggedDecaysProc_MC_vtosll_Dsp_pieta.root";
    signalFileName_Dsm=signalFileName_Dsm+"pieta/DsTaggedDecaysProc_MC_vtosll_Dsm_pieta.root";
  }
  else if (decay=="pietaprime")
  {
    signalFileName_Dsp=signalFileName_Dsp+"pietaprime/DsTaggedDecaysProc_MC_vtosll_Dsp_pietaprime.root";
    signalFileName_Dsm=signalFileName_Dsm+"pietaprime/DsTaggedDecaysProc_MC_vtosll_Dsm_pietaprime.root";
  }
  else if (decay=="KKpipi0")
  {
    signalFileName_Dsp=signalFileName_Dsp+"KKpipi0/DsTaggedDecaysProc_MC_vtosll_Dsp_KKpipi0.root";
    signalFileName_Dsm=signalFileName_Dsm+"KKpipi0/DsTaggedDecaysProc_MC_vtosll_Dsm_KKpipi0.root";
  }
  else if (decay=="pipipi")
  {
    signalFileName_Dsp=signalFileName_Dsp+"pipipi/DsTaggedDecaysProc_MC_vtosll_Dsp_pipipi.root";
    signalFileName_Dsm=signalFileName_Dsm+"pipipi/DsTaggedDecaysProc_MC_vtosll_Dsm_pipipi.root";
  }
  else if (decay=="KsKmpipi")
  {
    signalFileName_Dsp=signalFileName_Dsp+"KsKmpipi/DsTaggedDecaysProc_MC_vtosll_Dsp_KsKmpipi.root";
    signalFileName_Dsm=signalFileName_Dsm+"KsKmpipi/DsTaggedDecaysProc_MC_vtosll_Dsm_KsKmpipi.root";
  }
  else if (decay=="pipi0eta")
  {
    signalFileName_Dsp=signalFileName_Dsp+"pipi0eta/DsTaggedDecaysProc_MC_vtosll_Dsp_pipi0eta.root";
    signalFileName_Dsm=signalFileName_Dsm+"pipi0eta/DsTaggedDecaysProc_MC_vtosll_Dsm_pipi0eta.root";
  }
  else if (decay=="pietaprimerho")
  {
    signalFileName_Dsp=signalFileName_Dsp+"pietaprimerho/DsTaggedDecaysProc_MC_vtosll_Dsp_pietaprimerho.root";
    signalFileName_Dsm=signalFileName_Dsm+"pietaprimerho/DsTaggedDecaysProc_MC_vtosll_Dsm_pietaprimerho.root";
  }
}

double deltaPhi(double phi1, double phi2)
{
  double dphi=phi1-phi2;
  if (dphi<-pi) dphi=2*pi+dphi;
  if (dphi>pi) dphi=2*pi-dphi;
  return dphi;
}

bool dsPlusMCut(float dsPlusM)
{
  return (fabs(dsPlusM-dsPlusMCut_center)<dsPlusMCut_range);
}

bool MBCCut(float MBC)
{
  return (fabs(MBC-mbcCut_center)<mbcCut_range);
}

bool DeltaMCut(float DeltaM)
{
  return (fabs(DeltaM-deltaMCut_center)<deltaMCut_range);
}

bool dD0(float diff)
{
  return ((diff)>diffD0Cut);
}

bool dPhiCut(float dPhi)
{
  return (dPhi<dPhiCutLess);
}

bool electronEnergyCut(double electronEnergy, double positronEnergy)
{
  return (electronEnergy<electronEnergyThreshold && positronEnergy<electronEnergyThreshold);
}

bool pi0MassCut(float pi0Mass)
{
  return (fabs(pi0Mass-pi0MassCut_center)>pi0MassCut_range);
}

struct NEvents {
  float noCut;
  float tagCut;
  float dsPlusMCut;
  float deltaECut;
  float mbcCut;
  float deltaMCut;
  float diffD0Cut;
  float dPhiCut;
  float pi0MassCut;
};

struct EventNumber {
  float noCut_run;      float noCut_event;
  float tagCut_run;     float tagCut_event;
  float dsPlusMCut_run; float dsPlusMCut_event;
  float deltaECut_run;  float deltaECut_event;
  float mbcCut_run;     float mbcCut_event;
  float deltaMCut_run;  float deltaMCut_event;
  float diffD0Cut_run;  float diffD0Cut_event;
  float dPhiCut_run;    float dPhiCut_event;
  float pi0MassCut_run; float pi0MassCut_event;
};

int DsTaggedAnalysis_Optimize_pionFit() // (bool f1, bool f2, bool f3, bool f4, bool f5)
{
  /*
  DsMass_optimize=f1;
  mBC_optimize=f2;
  DeltaM_optimize=f3;
  diffD0_optimize=f4;
  dPhi_optimize=f5;
  */
  
  setFileNames();
  setValues();
  
  TChain *signalTree=new TChain("DsTaggedDecaysProc/nt1");
  signalTree->Add(signalFileName_Dsp.c_str());
  signalTree->Add(signalFileName_Dsm.c_str());  // Make sure they're in the same ratio, or the ratio is at least known.
  
  TChain *genericTree = new TChain("DsTaggedDecaysProc/nt1");
  genericTree->Add("/nfs/cor/an3/souvik/Dataset39/DsTaggedProc_GenericMC_213586_214863.root");
  genericTree->Add("/nfs/cor/an3/souvik/Dataset40/DsTaggedProc_GenericMC_215307_217385.root");
  genericTree->Add("/nfs/cor/an3/souvik/Dataset41/DsTaggedProc_GenericMC_217687_219721.root");
  genericTree->Add("/nfs/cor/an3/souvik/Dataset47/DsTaggedProc_GenericMC_230474_232255.root");
  genericTree->Add("/nfs/cor/an3/souvik/Dataset48/DsTaggedProc_GenericMC_232264_234607.root");
  
  TChain *continuTree = new TChain("DsTaggedDecaysProc/nt1");
  continuTree->Add("/nfs/cor/an3/souvik/Dataset39/DsTaggedProc_ContinuumMC_213586_214863.root");
  continuTree->Add("/nfs/cor/an3/souvik/Dataset40/DsTaggedProc_ContinuumMC_215307_217385.root");
  continuTree->Add("/nfs/cor/an3/souvik/Dataset41/DsTaggedProc_ContinuumMC_217687_219721.root");
  continuTree->Add("/nfs/cor/an3/souvik/Dataset47/DsTaggedProc_ContinuumMC_230474_232255.root");
  continuTree->Add("/nfs/cor/an3/souvik/Dataset48/DsTaggedProc_ContinuumMC_232264_234607.root");
  
  NEvents signalEvents={0,0,0,0,0,0,0,0,0};
  NEvents genericEvents={0,0,0,0,0,0,0,0,0};
  NEvents continuEvents={0,0,0,0,0,0,0,0,0};
  
  EventNumber signalNumber={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  EventNumber genericNumber={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  EventNumber continuNumber={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  float runNumber_signal, eventNumber_signal;
  float runNumber_generic, eventNumber_generic;
  float runNumber_continu, eventNumber_continu;
  float dsPlusM_signal, dsPlusCharge_signal, DeltaE_signal, MBC_signal, DeltaM_signal, decayMode_signal;
  float dsPlusM_generic, dsPlusCharge_generic, DeltaE_generic, MBC_generic, DeltaM_generic, decayMode_generic;
  float dsPlusM_continu, dsPlusCharge_continu, DeltaE_continu, MBC_continu, DeltaM_continu, decayMode_continu;
  float d0_e_signal, d0_p_signal, z0_e_signal, z0_p_signal, px_e_signal, py_e_signal, pz_e_signal, px_p_signal, py_p_signal, pz_p_signal, E_e_signal, E_p_signal, curv_e_signal, curv_p_signal, E_e_signal_MC, E_p_signal_MC;
  float d0_e_generic, d0_p_generic, z0_e_generic, z0_p_generic, px_e_generic, py_e_generic, pz_e_generic, px_p_generic, py_p_generic, pz_p_generic, E_e_generic, E_p_generic, curv_e_generic, curv_p_generic;
  float d0_e_continu, d0_p_continu, z0_e_continu, z0_p_continu, px_e_continu, py_e_continu, pz_e_continu, px_p_continu, py_p_continu, pz_p_continu, E_e_continu, E_p_continu, curv_e_continu, curv_p_continu;
  float pi0Mass_signal;
  float pi0Mass_generic;
  float pi0Mass_continu;
  
  TH1D *h_dsPlusM_signal = new TH1D("h_dsPlusM_signal", "m_{D_{S}^{#pm}} Signal Sample; m_{D_{S}^{#pm}} (GeV); # Events", 100, 1.94, 2.0); h_dsPlusM_signal->SetLineColor(kCyan);
  TH1D *h_dsPlusM_generic = new TH1D("h_dsPlusM_generic", "m_{D_{S}^{#pm}} Generic MC Background Sample; m_{D_{S}^{#pm}} (GeV); # Events", 100, 1.9, 2.05); h_dsPlusM_generic->SetLineColor(kGreen);
  TH1D *h_dsPlusM_continu = new TH1D("h_dsPlusM_continu", "m_{D_{S}^{#pm}} Continuum MC Background Sample; m_{D_{S}^{#pm}} (GeV); # Events", 100, 1.9, 2.05); h_dsPlusM_continu->SetLineColor(kBlue);
  TH1D *h_dsPlusM_signal_range = new TH1D("h_dsPlusM_signal_range", "m_{D_{S}^{#pm}} Signal Sample vs Cut Width; Cut Width (GeV); # Events", 20, 0., 0.02); h_dsPlusM_signal_range->SetLineColor(kCyan);
  TH1D *h_dsPlusM_generic_range = new TH1D("h_dsPlusM_generic_range", "m_{D_{S}^{#pm}} Generic MC Background Sample vs Cut Width; m_{D_{S}^{#pm}} GeV; # Events", 20, 0., 0.02); h_dsPlusM_generic_range->SetLineColor(kGreen);
  TH1D *h_dsPlusM_continu_range = new TH1D("h_dsPlusM_continu_range", "m_{D_{S}^{#pm}} Continuum MC Background Sample vs Cut Width; Cut Width GeV; # Events", 20, 0., 0.02); h_dsPlusM_continu_range->SetLineColor(kBlue);
  TH1D *h_dsPlusM_significance = new TH1D("h_dsPlusM_significance", "m_{D_{S}^{#pm}} Signal Significance vs Cut Width; Cut Width (GeV); #sigma", 20, 0., 0.02); h_dsPlusM_significance->SetLineColor(kBlack);
  TH1D *h_dsPlusM_precision = new TH1D("h_dsPlusM_precision", "m_{D_{S}^{#pm}} Signal Precision vs Cut Width; Cut Width (GeV)", 20, 0., 0.02); h_dsPlusM_precision->SetLineColor(kBlack);
  /*
  TH1D *h_dsPlusM_generic = new TH1D("h_dsPlusM_generic", "m_{D_{S}^{+}} Generic MC Background Sample; GeV", 100, 1.9, 2.05); h_dsPlusM_generic->SetLineColor(kBlue);
  TH1D *h_dsPlusM_continu = new TH1D("h_dsPlusM_continu", "m_{D_{S}^{+}} Continuum MC Background Sample; GeV", 100, 1.9, 2.05); h_dsPlusM_continu->SetLineColor(kGreen);
  TH1D *h_dsPlusM_signal_range = new TH1D("h_dsPlusM_signal_range", "m_{D_{S}^{+}} Signal Sample vs Range; GeV", 20, 0., 0.02); h_dsPlusM_signal_range->SetLineColor(kRed);
  TH1D *h_dsPlusM_generic_range = new TH1D("h_dsPlusM_generic_range", "m_{D_{S}^{+}} Generic MC Background Sample vs Range; GeV", 20, 0., 0.02); h_dsPlusM_generic_range->SetLineColor(kBlue);
  TH1D *h_dsPlusM_continu_range = new TH1D("h_dsPlusM_continu_range", "m_{D_{S}^{+}} Continuum MC Background Sample; GeV", 20, 0., 0.02); h_dsPlusM_continu_range->SetLineColor(kGreen);
  TH1D *h_dsPlusM_significance = new TH1D("h_dsPlusM_significance", "m_{D_{S}^{+}} Signal Significance vs Range", 20, 0., 0.02); h_dsPlusM_significance->SetLineColor(kBlack);
  TH1D *h_dsPlusM_precision = new TH1D("h_dsPlusM_precision", "m_{D_{S}^{+}} Signal Precision", 20, 0., 0.02); h_dsPlusM_precision->SetLineColor(kBlack);
  */
  
  TH1D *h_MBC_signal = new TH1D("h_MBC_signal", "m_{BC} Signal Sample; m_{BC} (GeV); # Events", 100, 2.1, 2.14); h_MBC_signal->SetLineColor(kCyan);
  TH1D *h_MBC_generic = new TH1D("h_MBC_generic", "m_{BC} Generic MC Background; m_{BC} (GeV); # Events", 100, 2.04, 2.16); h_MBC_generic->SetLineColor(kGreen);
  TH1D *h_MBC_continu = new TH1D("h_MBC_continu", "m_{BC} Continuum MC Background; m_{BC} (GeV); # Events", 100, 2.04, 2.16); h_MBC_continu->SetLineColor(kBlue);
  TH1D *h_MBC_signal_range = new TH1D("h_MBC_signal_range", "m_{BC} Signal Sample vs Cut Width; Cut Width (GeV); # Events", 20, 0., 0.02); h_MBC_signal_range->SetLineColor(kCyan);
  TH1D *h_MBC_generic_range = new TH1D("h_MBC_generic_range", "m_{BC} Generic MC Background vs Cut Width; Cut Width (GeV); # Events", 20, 0., 0.02); h_MBC_generic_range->SetLineColor(kGreen);
  TH1D *h_MBC_continu_range = new TH1D("h_MBC_continu_range", "m_{BC} Continuum MC Background vs Cut Width; Cut Width (GeV); # Events", 20, 0., 0.02); h_MBC_continu_range->SetLineColor(kBlue);
  TH1D *h_MBC_significance = new TH1D("h_MBC_significance", "m_{BC} Signal Significance vs Cut Width; Cut Width (GeV); # Events", 20, 0., 0.02); h_MBC_significance->SetLineColor(kBlack);
  TH1D *h_MBC_precision = new TH1D("h_MBC_precision", "m_{BC} Signal Precision vs Cut Width; Cut Width (GeV); # Events", 20, 0., 0.02); h_MBC_precision->SetLineColor(kBlack);
  /*
  TH1D *h_MBC_signal = new TH1D("h_MBC_signal", "m_{BC} Signal Sample; GeV", 100, 2.1, 2.14); h_MBC_signal->SetLineColor(kRed);
  TH1D *h_MBC_generic = new TH1D("h_MBC_generic", "m_{BC} Generic MC Background; GeV", 100, 2.04, 2.16); h_MBC_generic->SetLineColor(kBlue);
  TH1D *h_MBC_continu = new TH1D("h_MBC_continu", "m_{BC} Continuum MC Background; GeV", 100, 2.04, 2.16); h_MBC_continu->SetLineColor(kGreen);
  TH1D *h_MBC_signal_range = new TH1D("h_MBC_signal_range", "m_{BC} Signal Sample vs Range; GeV", 20, 0., 0.02); h_MBC_signal_range->SetLineColor(kRed);
  TH1D *h_MBC_generic_range = new TH1D("h_MBC_generic_range", "m_{BC} Generic MC Background vs Range; GeV", 20, 0., 0.02); h_MBC_generic_range->SetLineColor(kBlue);
  TH1D *h_MBC_continu_range = new TH1D("h_MBC_continu_range", "m_{BC} Continuum MC Background vs Range; GeV", 20, 0., 0.02); h_MBC_continu_range->SetLineColor(kGreen);
  TH1D *h_MBC_significance = new TH1D("h_MBC_significance", "m_{BC} Signal Significance vs Range", 20, 0., 0.02); h_MBC_significance->SetLineColor(kBlack);
  TH1D *h_MBC_precision = new TH1D("h_MBC_precision", "m_{BC} Signal Precision vs Range", 20, 0., 0.02); h_MBC_precision->SetLineColor(kBlack);
  */
  
  TH1D *h_DeltaM_signal = new TH1D("h_DeltaM_signal", "#deltam Signal Sample; #deltam (GeV); # Events", 100, 0.12, 0.16); h_DeltaM_signal->SetLineColor(kCyan);
  TH1D *h_DeltaM_generic = new TH1D("h_DeltaM_generic", "#deltam Generic MC Background Sample; #deltam (GeV); # Events", 100, 0.08, 0.2); h_DeltaM_generic->SetLineColor(kGreen);
  TH1D *h_DeltaM_continu = new TH1D("h_DeltaM_continu", "#deltam Contiuum MC Background Sample; #deltam (GeV); # Events", 100, 0.08, 0.2); h_DeltaM_continu->SetLineColor(kBlue);
  TH1D *h_DeltaM_signal_range = new TH1D("h_DeltaM_signal_range", "#deltam Signal Sample vs Cut Width; Cut Width (GeV); # Events", 20, 0.0, 0.02); h_DeltaM_signal_range->SetLineColor(kCyan);
  TH1D *h_DeltaM_generic_range = new TH1D("h_DeltaM_generic_range", "#deltam Generic MC Background vs Cut Width; Cut Width (GeV); # Events", 20, 0.0, 0.02); h_DeltaM_generic_range->SetLineColor(kGreen);
  TH1D *h_DeltaM_continu_range = new TH1D("h_DeltaM_continu_range", "#deltam Continuum MC Background vs Cut Width; Cut Width (GeV); # Events", 20, 0.0, 0.02); h_DeltaM_continu_range->SetLineColor(kBlue);
  TH1D *h_DeltaM_significance = new TH1D("h_DeltaM_significance", "#deltam Signal Significance vs Cut Width; Cut Width (GeV); # Events", 20, 0.0, 0.02); h_DeltaM_significance->SetLineColor(kBlack);
  TH1D *h_DeltaM_precision = new TH1D("h_DeltaM_precision", "#deltam Signal Precision vs Cut Width; Cut Width (GeV); # Events", 20, 0.0, 0.02); h_DeltaM_precision->SetLineColor(kBlack);
  /*
  TH1D *h_DeltaM_signal = new TH1D("h_DeltaM_signal", "#deltaM Signal Sample; GeV", 100, 0.12, 0.2); h_DeltaM_signal->SetLineColor(kRed);
  TH1D *h_DeltaM_generic = new TH1D("h_DeltaM_generic", "#deltaM Generic MC Background Sample; GeV", 100, 0.08, 0.2); h_DeltaM_generic->SetLineColor(kBlue);
  TH1D *h_DeltaM_continu = new TH1D("h_DeltaM_continu", "#deltaM Contiuum MC Background Sample; GeV", 100, 0.08, 0.2); h_DeltaM_continu->SetLineColor(kGreen);
  TH1D *h_DeltaM_signal_range = new TH1D("h_DeltaM_signal_range", "#deltaM Signal Sample vs Range; GeV", 20, 0.0, 0.02); h_DeltaM_signal_range->SetLineColor(kRed);
  TH1D *h_DeltaM_generic_range = new TH1D("h_DeltaM_generic_range", "#deltaM Generic MC Background vs Range; GeV", 20, 0.0, 0.02); h_DeltaM_generic_range->SetLineColor(kBlue);
  TH1D *h_DeltaM_continu_range = new TH1D("h_DeltaM_continu_range", "#deltaM Continuum MC Background vs Range; GeV", 20, 0.0, 0.02); h_DeltaM_continu_range->SetLineColor(kGreen);
  TH1D *h_DeltaM_significance = new TH1D("h_DeltaM_significance", "#deltaM Signal Significance vs Range; GeV", 20, 0.0, 0.02); h_DeltaM_significance->SetLineColor(kBlack);
  TH1D *h_DeltaM_precision = new TH1D("h_DeltaM_precision", "#deltaM Signal Precision vs Range; GeV", 20, 0.0, 0.02); h_DeltaM_precision->SetLineColor(kBlack);
  */
  
  TH1D *h_diffD0_signal = new TH1D("h_diffD0_signal", "#Deltad_{0} Signal Sample; #Deltad_{0} (m); # Events", 50, -0.01, 0.01); h_diffD0_signal->SetLineColor(kCyan);
  TH1D *h_diffD0_generic = new TH1D("h_diffD0_generic", "#Deltad_{0} Generic MC Background Sample; #Deltad_{0} (m); # Events", 50, -0.01, 0.01); h_diffD0_generic->SetLineColor(kGreen);
  TH1D *h_diffD0_continu = new TH1D("h_diffD0_continu", "#Deltad_{0} Continuum MC Background Sample; #Deltad_{0} (m); # Events", 50, -0.01, 0.01); h_diffD0_continu->SetLineColor(kBlue);
  TH1D *h_diffD0_signal_range = new TH1D("h_diffD0_signal_range", "#Deltad_{0} Signal Sample vs Cut Width; Cut Width (m); # Events", 20, -0.01, 0.01); h_diffD0_signal_range->SetLineColor(kCyan);
  TH1D *h_diffD0_generic_range = new TH1D("h_diffD0_generic_range", "#Deltad_{0} Generic MC Background Sample vs Cut Width; Cut Width (m); # Events", 20, -0.01, 0.01); h_diffD0_generic_range->SetLineColor(kGreen);
  TH1D *h_diffD0_continu_range = new TH1D("h_diffD0_continu_range", "#Deltad_{0} Continuum MC Background Sample vs Cut Width; Cut Width (m); # Events", 20, -0.01, 0.01); h_diffD0_continu_range->SetLineColor(kBlue);
  TH1D *h_diffD0_significance = new TH1D("h_diffD0_significance", "#Deltad_{0} Significance vs Cut Width; Cut Width (m); # Events", 20, -0.01, 0.01); h_diffD0_significance->SetLineColor(kBlack);
  TH1D *h_diffD0_precision = new TH1D("h_diffD0_precision", "#Deltad_{0} Precision vs Cut Width; Cut Width (m); # Events", 20, -0.01, 0.01); h_diffD0_precision->SetLineColor(kBlack);
  /*
  TH1D *h_diffD0_signal = new TH1D("h_diffD0_signal", "#Deltad_{0} Signal Sample; m", 50, -0.01, 0.01); h_diffD0_signal->SetLineColor(kRed); h_diffD0_signal->GetXaxis()->SetLabelSize(0.02);
  TH1D *h_diffD0_generic = new TH1D("h_diffD0_generic", "#Deltad_{0} Generic MC Background Sample; m", 50, -0.01, 0.01); h_diffD0_generic->SetLineColor(kBlue); h_diffD0_generic->SetLabelSize(0.02);
  TH1D *h_diffD0_continu = new TH1D("h_diffD0_continu", "#Deltad_{0} Continuum MC Background Sample; m", 50, -0.01, 0.01); h_diffD0_continu->SetLineColor(kGreen); h_diffD0_continu->SetLabelSize(0.02);
  TH1D *h_diffD0_signal_range = new TH1D("h_diffD0_signal_range", "#Deltad_{0} Signal Sample vs Range; m", 20, -0.01, 0.01); h_diffD0_signal_range->SetLineColor(kRed);
  TH1D *h_diffD0_generic_range = new TH1D("h_diffD0_generic_range", "#Deltad_{0} Generic MC Background Sample vs Range; m", 20, -0.01, 0.01); h_diffD0_generic_range->SetLineColor(kBlue);
  TH1D *h_diffD0_continu_range = new TH1D("h_diffD0_continu_range", "#Deltad_{0} Continuum MC Background Sample vs Range; m", 20, -0.01, 0.01); h_diffD0_continu_range->SetLineColor(kGreen);
  TH1D *h_diffD0_significance = new TH1D("h_diffD0_significance", "#Deltad_{0} Significance vs Range; m", 20, -0.01, 0.01); h_diffD0_significance->SetLineColor(kBlack);
  TH1D *h_diffD0_precision = new TH1D("h_diffD0_precision", "#Deltad_{0} Precision vs Range; m", 20, -0.01, 0.01); h_diffD0_precision->SetLineColor(kBlack);
  */
  
  TH1D *h_dPhi_signal = new TH1D("h_dPhi_signal", "#Delta#Phi Signal Sample; #Delta#Phi (GeV); # Events", 50, -pi/2, pi/2); h_dPhi_signal->SetLineColor(kCyan);
  TH1D *h_dPhi_generic = new TH1D("h_dPhi_generic", "#Delta#Phi Generic MC Background; #Delta#Phi (GeV); # Events", 50, -pi/2, pi/2); h_dPhi_generic->SetLineColor(kGreen);
  TH1D *h_dPhi_continu = new TH1D("h_dPhi_continu", "#Delta#Phi Continuum MC Background; #Delta#Phi (GeV); # Events", 50, -pi/2, pi/2); h_dPhi_continu->SetLineColor(kBlue);
  TH1D *h_dPhi_signal_range = new TH1D("h_dPhi_signal_range", "#Delta#Phi Signal Sample vs Cut Width; Cut Width (GeV); # Events", 30, -0.1, 0.2); h_dPhi_signal_range->SetLineColor(kRed);
  TH1D *h_dPhi_generic_range = new TH1D("h_dPhi_generic_range", "#Delta#Phi Generic MC Background vs Cut Width; Cut Width (GeV); # Events", 30, -0.1, 0.2); h_dPhi_generic_range->SetLineColor(kBlue);
  TH1D *h_dPhi_continu_range = new TH1D("h_dPhi_continu_range", "#Delta#Phi Continuum MC Background vs Cut Width; Cut Width (GeV); # Events", 30, -0.1, 0.2); h_dPhi_continu_range->SetLineColor(kGreen);
  TH1D *h_dPhi_significance = new TH1D("h_dPhi_significance", "#Delta#Phi Significance vs Cut Width; Cut Width (GeV); # Events", 30, -0.1, 0.2); h_dPhi_significance->SetLineColor(kBlack);
  TH1D *h_dPhi_precision = new TH1D("h_dPhi_precision", "#Delta#Phi Precision vs Cut Width; Cut Width (GeV); # Events", 30, -0.1, 0.2); h_dPhi_precision->SetLineColor(kBlack);
  /*
  TH1D *h_dPhi_signal = new TH1D("h_dPhi_signal", "#Delta#Phi Signal Sample", 50, -1., 1.); h_dPhi_signal->SetLineColor(kRed);
  TH1D *h_dPhi_generic = new TH1D("h_dPhi_generic", "#Delta#Phi Generic MC Background", 50, -1., 1.); h_dPhi_generic->SetLineColor(kBlue);
  TH1D *h_dPhi_continu = new TH1D("h_dPhi_continu", "#Delta#Phi Continuum MC Background", 50, -1., 1.); h_dPhi_continu->SetLineColor(kGreen);
  TH1D *h_dPhi_signal_range = new TH1D("h_dPhi_signal_range", "#Delta#Phi Signal Sample vs Range", 30, -0.1, 0.2); h_dPhi_signal_range->SetLineColor(kRed);
  TH1D *h_dPhi_generic_range = new TH1D("h_dPhi_generic_range", "#Delta#Phi Generic MC Background vs Range", 30, -0.1, 0.2); h_dPhi_generic_range->SetLineColor(kBlue);
  TH1D *h_dPhi_continu_range = new TH1D("h_dPhi_continu_range", "#Delta#Phi Continuum MC Background vs Range", 30, -0.1, 0.2); h_dPhi_continu_range->SetLineColor(kGreen);
  TH1D *h_dPhi_significance = new TH1D("h_dPhi_significance", "#Delta#Phi Significance vs Range", 30, -0.1, 0.2); h_dPhi_significance->SetLineColor(kBlack);
  TH1D *h_dPhi_precision = new TH1D("h_dPhi_precision", "#Delta#Phi Precision vs Range", 30, -0.1, 0.2); h_dPhi_precision->SetLineColor(kBlack);
  */
  
  TH1D *h_pi0_signal = new TH1D("h_pi0_signal", "Pion Mass; GeV", 50, 0.035, .235); h_pi0_signal->SetLineColor(kRed);
  TH1D *h_pi0_generic = new TH1D("h_pi0_generic", "Pion Mass; GeV", 50, 0.035, .235); h_pi0_generic->SetLineColor(kBlue);
  TH1D *h_pi0_continu = new TH1D("h_pi0_continu", "Pion Mass; Gev", 50, 0.035, .235); h_pi0_continu->SetLineColor(kGreen);
  TH1D *h_pi0_signal_range = new TH1D("h_pi0_signal_range", "Pion Mass vs Range; Gev", 20, 0., .02); h_pi0_signal_range->SetLineColor(kRed);
  TH1D *h_pi0_generic_range = new TH1D("h_pi0_generic_range", "Pion Mass vs Range; Gev", 20, 0., .02); h_pi0_generic_range->SetLineColor(kBlue);
  TH1D *h_pi0_continu_range = new TH1D("h_pi0_continu_range", "Pion Mass vs Ranve; Gev", 20, 0., .02); h_pi0_continu_range->SetLineColor(kGreen);
  TH1D *h_pi0_significance = new TH1D("h_pi0_significance", "Pion Mass Cut Significance vs Range; GeV", 20, 0., .02); h_pi0_significance->SetLineColor(kGreen);
  TH1D *h_pi0_precision = new TH1D("h_pi0_precision", "Pion Mass Cut Precision vs Range; GeV", 20, 0., .02); h_pi0_precision->SetLineColor(kGreen);
  
  TH2D *h_electronE = new TH2D("h_electronE", "Electron Energy Resolution; MC Electron Energy (GeV); Reco-MC Electron Energy (GeV)", 100, 0.0, 0.16, 100, -0.03, 0.03);
  
  signalTree->SetBranchAddress("Run", &(runNumber_signal));
  signalTree->SetBranchAddress("Event", &(eventNumber_signal));
  signalTree->SetBranchAddress("dsPlusM", &(dsPlusM_signal));
  signalTree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_signal));
  signalTree->SetBranchAddress("DecayMode", &(decayMode_signal));
  signalTree->SetBranchAddress("DeltaE", &(DeltaE_signal));
  signalTree->SetBranchAddress("MBC", &(MBC_signal));
  signalTree->SetBranchAddress("DeltaM", &(DeltaM_signal));
  signalTree->SetBranchAddress("kElectron1D0_reco", &(d0_e_signal));
  signalTree->SetBranchAddress("kElectron2D0_reco", &(d0_p_signal));
  signalTree->SetBranchAddress("kElectron1Z0_reco", &(z0_e_signal));
  signalTree->SetBranchAddress("kElectron2Z0_reco", &(z0_p_signal));
  signalTree->SetBranchAddress("kElectron1Px_reco", &(px_e_signal));
  signalTree->SetBranchAddress("kElectron1Py_reco", &(py_e_signal));
  signalTree->SetBranchAddress("kElectron1Pz_reco", &(pz_e_signal));
  signalTree->SetBranchAddress("kElectron2Px_reco", &(px_p_signal));
  signalTree->SetBranchAddress("kElectron2Py_reco", &(py_p_signal));
  signalTree->SetBranchAddress("kElectron2Pz_reco", &(pz_p_signal));
  signalTree->SetBranchAddress("kElectron1E_reco", &(E_e_signal));
  signalTree->SetBranchAddress("kElectron2E_reco", &(E_p_signal));
  signalTree->SetBranchAddress("kElectron1Curv_reco", &(curv_e_signal));
  signalTree->SetBranchAddress("kElectron2Curv_reco", &(curv_p_signal));
  signalTree->SetBranchAddress("kPi0Mass_reco", &(pi0Mass_signal));
  signalTree->SetBranchAddress("kElectron1E_MC", &(E_e_signal_MC));
  signalTree->SetBranchAddress("kElectron2E_MC", &(E_p_signal_MC));
  
  genericTree->SetBranchAddress("Run", &(runNumber_generic));
  genericTree->SetBranchAddress("Event", &(eventNumber_generic));
  genericTree->SetBranchAddress("dsPlusM", &(dsPlusM_generic));
  genericTree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_generic));
  genericTree->SetBranchAddress("DecayMode", &(decayMode_generic));
  genericTree->SetBranchAddress("DeltaE", &(DeltaE_generic));
  genericTree->SetBranchAddress("MBC", &(MBC_generic));
  genericTree->SetBranchAddress("DeltaM", &(DeltaM_generic));
  genericTree->SetBranchAddress("kElectron1D0_reco", &(d0_e_generic));
  genericTree->SetBranchAddress("kElectron2D0_reco", &(d0_p_generic));
  genericTree->SetBranchAddress("kElectron1Z0_reco", &(z0_e_generic));
  genericTree->SetBranchAddress("kElectron2Z0_reco", &(z0_p_generic));
  genericTree->SetBranchAddress("kElectron1Px_reco", &(px_e_generic));
  genericTree->SetBranchAddress("kElectron1Py_reco", &(py_e_generic));
  genericTree->SetBranchAddress("kElectron1Pz_reco", &(pz_e_generic));
  genericTree->SetBranchAddress("kElectron2Px_reco", &(px_p_generic));
  genericTree->SetBranchAddress("kElectron2Py_reco", &(py_p_generic));
  genericTree->SetBranchAddress("kElectron2Pz_reco", &(pz_p_generic));
  genericTree->SetBranchAddress("kElectron1E_reco", &(E_e_generic));
  genericTree->SetBranchAddress("kElectron2E_reco", &(E_p_generic));
  genericTree->SetBranchAddress("kElectron1Curv_reco", &(curv_e_generic));
  genericTree->SetBranchAddress("kElectron2Curv_reco", &(curv_p_generic));
  genericTree->SetBranchAddress("kPi0Mass_reco", &(pi0Mass_generic));
  
  continuTree->SetBranchAddress("Run", &(runNumber_continu));
  continuTree->SetBranchAddress("Event", &(eventNumber_continu));
  continuTree->SetBranchAddress("dsPlusM", &(dsPlusM_continu));
  continuTree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_continu));
  continuTree->SetBranchAddress("DecayMode", &(decayMode_continu));
  continuTree->SetBranchAddress("DeltaE", &(DeltaE_continu));
  continuTree->SetBranchAddress("MBC", &(MBC_continu));
  continuTree->SetBranchAddress("DeltaM", &(DeltaM_continu));
  continuTree->SetBranchAddress("kElectron1D0_reco", &(d0_e_continu));
  continuTree->SetBranchAddress("kElectron2D0_reco", &(d0_p_continu));
  continuTree->SetBranchAddress("kElectron1Z0_reco", &(z0_e_continu));
  continuTree->SetBranchAddress("kElectron2Z0_reco", &(z0_p_continu));
  continuTree->SetBranchAddress("kElectron1Px_reco", &(px_e_continu));
  continuTree->SetBranchAddress("kElectron1Py_reco", &(py_e_continu));
  continuTree->SetBranchAddress("kElectron2Px_reco", &(px_p_continu));
  continuTree->SetBranchAddress("kElectron2Py_reco", &(py_p_continu));
  continuTree->SetBranchAddress("kElectron1E_reco", &(E_e_continu));
  continuTree->SetBranchAddress("kElectron2E_reco", &(E_p_continu));
  continuTree->SetBranchAddress("kElectron1Curv_reco", &(curv_e_continu));
  continuTree->SetBranchAddress("kElectron2Curv_reco", &(curv_p_continu));
  continuTree->SetBranchAddress("kPi0Mass_reco", &(pi0Mass_continu));
  

  int nSignalEvents=signalTree->GetEntries();
  for (int i=0; i<nSignalEvents; ++i)
  {
    signalTree->GetEvent(i);
    
    double phi_e=atan2(py_e_signal, px_e_signal);
    double phi_p=atan2(py_p_signal, px_p_signal);
    double dPhi=deltaPhi(phi_e, phi_p);
    
    if (runNumber_signal!=signalNumber.noCut_run || eventNumber_signal!=signalNumber.noCut_event)
    {
      ++signalEvents.noCut;
      signalNumber.noCut_run=runNumber_signal;
      signalNumber.noCut_event=eventNumber_signal;
    }
    
    if (fabs(d0_e_signal)<0.005 && fabs(d0_p_signal)<0.005 &&
        fabs(z0_e_signal)<0.05 && fabs(z0_p_signal)<0.05 &&
        electronEnergyCut(E_e_signal, E_p_signal))
    {
      if (decayMode_signal==decayNumber)
      {
        if (!DsMass_optimize) h_dsPlusM_signal->Fill(dsPlusM_signal);
        if (DsMass_optimize || dsPlusMCut(dsPlusM_signal))
        {
          if (!mBC_optimize) h_MBC_signal->Fill(MBC_signal);
          if (mBC_optimize || MBCCut(MBC_signal))
          {
            if (!DeltaM_optimize) h_DeltaM_signal->Fill(DeltaM_signal);
            if (DeltaM_optimize || DeltaMCut(DeltaM_signal))
            {
              if (!diffD0_optimize) h_diffD0_signal->Fill(d0_e_signal-d0_p_signal);
              if (diffD0_optimize || dD0(d0_e_signal-d0_p_signal))
              {
                if (!dPhi_optimize) h_dPhi_signal->Fill(dPhi);
                if (dPhi_optimize || dPhiCut(dPhi))
                {
                  if (!pi0_optimize) h_pi0_signal->Fill(pi0Mass_signal);
                  if (pi0_optimize || pi0MassCut(pi0Mass_signal))
                  {
                    if (runNumber_signal!=signalNumber.pi0MassCut_run || eventNumber_signal!=signalNumber.pi0MassCut_event)
                    {
                      ++signalEvents.pi0MassCut;
                      signalNumber.pi0MassCut_run=runNumber_signal;
                      signalNumber.pi0MassCut_event=eventNumber_signal;
                    }
                    if (DsMass_optimize)
                    {
                      h_dsPlusM_signal->Fill(dsPlusM_signal);
                      for (dsPlusMCut_range=0.0005; dsPlusMCut_range<0.0195; dsPlusMCut_range+=0.001)
                      {
                        if (dsPlusMCut(dsPlusM_signal)) h_dsPlusM_signal_range->Fill(dsPlusMCut_range);
                      }
                    }
                    if (mBC_optimize)
                    {
                      h_MBC_signal->Fill(MBC_signal);
                      for (mbcCut_range=0.0005; mbcCut_range<0.0195; mbcCut_range+=0.001)
                      {
                        if (MBCCut(MBC_signal)) h_MBC_signal_range->Fill(mbcCut_range);
                      }
                    }
                    if (DeltaM_optimize)
                    {
                      h_DeltaM_signal->Fill(DeltaM_signal);
                      for (deltaMCut_range=0.0005; deltaMCut_range<0.0195; deltaMCut_range+=0.001)
                      {
                        if (DeltaMCut(DeltaM_signal)) h_DeltaM_signal_range->Fill(deltaMCut_range);
                      }
                    }
                    if (diffD0_optimize)
                    {
                      h_diffD0_signal->Fill(d0_e_signal-d0_p_signal);
                      for (diffD0Cut=-0.0095; diffD0Cut<=0.0095; diffD0Cut+=0.001)
                      {
                        if (dD0(d0_e_signal-d0_p_signal)) h_diffD0_signal_range->Fill(diffD0Cut);
                      }
                    }
                    if (dPhi_optimize)
                    {
                      h_dPhi_signal->Fill(dPhi);
                      for (dPhiCutLess=-0.095; dPhiCutLess<=0.195; dPhiCutLess+=0.01)
                      {
                        if (dPhiCut(dPhi)) h_dPhi_signal_range->Fill(dPhiCutLess);
                      }
                    }
                    if (pi0_optimize)
                    {
                      h_pi0_signal->Fill(pi0Mass_signal);
                      for (pi0MassCut_range=0.0005; pi0MassCut_range<0.1995; pi0MassCut_range+=0.001)
                      {
                        if (pi0MassCut(pi0Mass_signal)) h_pi0_signal_range->Fill(pi0MassCut_range);
                      }
                    }
                    h_electronE->Fill(E_e_signal_MC, E_e_signal-E_e_signal_MC);
                    h_electronE->Fill(E_p_signal_MC, E_p_signal-E_p_signal_MC);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  int nGenericEvents=genericTree->GetEntries();
  for (int i=0; i<nGenericEvents; ++i)
  {
    genericTree->GetEvent(i);
    
    double phi_e=atan2(py_e_generic, px_e_generic);
    double phi_p=atan2(py_p_generic, px_p_generic);
    double dPhi=deltaPhi(phi_e, phi_p);
    
    if (runNumber_generic!=genericNumber.noCut_run || eventNumber_generic!=genericNumber.noCut_event)
    {
      ++genericEvents.noCut;
      genericNumber.noCut_run=runNumber_generic;
      genericNumber.noCut_event=eventNumber_generic;
    }
    
    if (fabs(d0_e_generic)<0.005 && fabs(d0_p_generic)<0.005 &&
        fabs(z0_e_generic)<0.05 && fabs(z0_p_generic)<0.05 &&
        electronEnergyCut(E_e_generic, E_p_generic))
    {
      if (decayMode_generic==decayNumber)
      {
        if (!DsMass_optimize) h_dsPlusM_generic->Fill(dsPlusM_generic);
        if (DsMass_optimize || dsPlusMCut(dsPlusM_generic))
        {
          if (!mBC_optimize) h_MBC_generic->Fill(MBC_generic);
          if (mBC_optimize || MBCCut(MBC_generic))
          {
            if (!DeltaM_optimize) h_DeltaM_generic->Fill(DeltaM_generic);
            if (DeltaM_optimize || DeltaMCut(DeltaM_generic))
            {
              if (!diffD0_optimize) h_diffD0_generic->Fill(d0_e_generic-d0_p_generic);
              if (diffD0_optimize || dD0(d0_e_generic-d0_p_generic))
              {
                if (!dPhi_optimize) h_dPhi_generic->Fill(dPhi);
                if (dPhi_optimize || dPhiCut(dPhi))
                {
                  if (!pi0_optimize) h_pi0_generic->Fill(pi0Mass_generic);
                  if (pi0_optimize || pi0MassCut(pi0Mass_generic))
                  {
                    if (runNumber_generic!=genericNumber.pi0MassCut_run || eventNumber_generic!=genericNumber.pi0MassCut_event)
                    {
                      ++genericEvents.pi0MassCut;
                      genericNumber.pi0MassCut_run=runNumber_generic;
                      genericNumber.pi0MassCut_event=eventNumber_generic;
                    }
                    if (DsMass_optimize)
                    {
                      h_dsPlusM_generic->Fill(dsPlusM_generic);
                      for (dsPlusMCut_range=0.0005; dsPlusMCut_range<0.0195; dsPlusMCut_range+=0.001)
                      {
                        if (dsPlusMCut(dsPlusM_generic)) h_dsPlusM_generic_range->Fill(dsPlusMCut_range);
                      }
                    }
                    if (mBC_optimize)
                    {
                      h_MBC_generic->Fill(MBC_generic);
                      for (mbcCut_range=0.0005; mbcCut_range<0.0195; mbcCut_range+=0.001)
                      {
                        if (MBCCut(MBC_generic)) h_MBC_generic_range->Fill(mbcCut_range);
                      }
                    }
                    if (DeltaM_optimize)
                    {
                      h_DeltaM_generic->Fill(DeltaM_generic);
                      for (deltaMCut_range=0.0005; deltaMCut_range<0.0195; deltaMCut_range+=0.001)
                      {
                        if (DeltaMCut(DeltaM_generic)) h_DeltaM_generic_range->Fill(deltaMCut_range);
                      }
                    }
                    if (diffD0_optimize)
                    {
                      h_diffD0_generic->Fill(d0_e_generic-d0_p_generic);
                      for (diffD0Cut=-0.0095; diffD0Cut<=0.0095; diffD0Cut+=0.001)
                      {
                        if (dD0(d0_e_generic-d0_p_generic)) h_diffD0_generic_range->Fill(diffD0Cut);
                      }
                    }
                    if (dPhi_optimize)
                    {
                      h_dPhi_generic->Fill(dPhi);
                      for (dPhiCutLess=-0.095; dPhiCutLess<=0.195; dPhiCutLess+=0.01)
                      {
                        if (dPhiCut(dPhi)) h_dPhi_generic_range->Fill(dPhiCutLess);
                      }
                    }
                    if (pi0_optimize)
                    {
                      h_pi0_generic->Fill(pi0Mass_generic);
                      for (pi0MassCut_range=0.0005; pi0MassCut_range<0.1995; pi0MassCut_range+=0.001)
                      {
                        if (pi0MassCut(pi0Mass_generic)) h_pi0_generic_range->Fill(pi0MassCut_range);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  int nContinuEvents=continuTree->GetEntries();
  for (int i=0; i<nContinuEvents; ++i)
  {
    continuTree->GetEvent(i);
    
    double phi_e=atan2(py_e_continu, px_e_continu);
    double phi_p=atan2(py_p_continu, px_p_continu);
    double dPhi=deltaPhi(phi_e, phi_p);
    
    if (runNumber_continu!=continuNumber.noCut_run || eventNumber_continu!=continuNumber.noCut_event)
    {
      ++continuEvents.noCut;
      continuNumber.noCut_run=runNumber_continu;
      continuNumber.noCut_event=eventNumber_continu;
    }
    
    if (fabs(d0_e_continu)<0.005 && fabs(d0_p_continu)<0.005 &&
        fabs(z0_e_continu)<0.05 && fabs(z0_p_continu)<0.05 &&
        electronEnergyCut(E_e_continu, E_p_continu))
    {
      if (decayMode_continu==decayNumber)
      {
        if (!DsMass_optimize) h_dsPlusM_continu->Fill(dsPlusM_continu);
        if (DsMass_optimize || dsPlusMCut(dsPlusM_continu))
        {
          if (!mBC_optimize) h_MBC_continu->Fill(MBC_continu);
          if (mBC_optimize || MBCCut(MBC_continu))
          {
            if (!DeltaM_optimize) h_DeltaM_continu->Fill(DeltaM_continu);
            if (DeltaM_optimize || DeltaMCut(DeltaM_continu))
            {
              if (!diffD0_optimize) h_diffD0_continu->Fill(d0_e_continu-d0_p_continu);
              if (diffD0_optimize || dD0(d0_e_continu-d0_p_continu))
              {
                if (!dPhi_optimize) h_dPhi_continu->Fill(dPhi);
                if (dPhi_optimize || dPhiCut(dPhi))
                {
                  if (!pi0_optimize) h_pi0_continu->Fill(pi0Mass_continu);
                  if (pi0_optimize || pi0MassCut(pi0Mass_continu))
                  {
                    if (runNumber_continu!=continuNumber.pi0MassCut_run || eventNumber_continu!=continuNumber.pi0MassCut_event)
                    {
                      ++continuEvents.pi0MassCut;
                      continuNumber.pi0MassCut_run=runNumber_continu;
                      continuNumber.pi0MassCut_event=eventNumber_continu;
                    }
                    if (DsMass_optimize)
                    {
                      h_dsPlusM_continu->Fill(dsPlusM_continu);
                      for (dsPlusMCut_range=0.0005; dsPlusMCut_range<0.0195; dsPlusMCut_range+=0.001)
                      {
                        if (dsPlusMCut(dsPlusM_continu)) h_dsPlusM_continu_range->Fill(dsPlusMCut_range);
                      }
                    }
                    if (mBC_optimize)
                    {
                      h_MBC_continu->Fill(MBC_continu);
                      for (mbcCut_range=0.0005; mbcCut_range<0.0195; mbcCut_range+=0.001)
                      {
                        if (MBCCut(MBC_continu)) h_MBC_continu_range->Fill(mbcCut_range);
                      }
                    }
                    if (DeltaM_optimize)
                    {
                      h_DeltaM_continu->Fill(DeltaM_continu);
                      for (deltaMCut_range=0.0005; deltaMCut_range<0.0195; deltaMCut_range+=0.001)
                      {
                        if (DeltaMCut(DeltaM_continu)) h_DeltaM_continu_range->Fill(deltaMCut_range);
                      }
                    }
                    if (diffD0_optimize)
                    {
                      h_diffD0_continu->Fill(d0_e_continu-d0_p_continu);
                      for (diffD0Cut=-0.0095; diffD0Cut<=0.0095; diffD0Cut+=0.001)
                      {
                        if (dD0(d0_e_continu-d0_p_continu)) h_diffD0_continu_range->Fill(diffD0Cut);
                      }
                    }
                    if (dPhi_optimize)
                    {
                      h_dPhi_continu->Fill(dPhi);
                      for (dPhiCutLess=-0.095; dPhiCutLess<=0.195; dPhiCutLess+=0.01)
                      {
                        if (dPhiCut(dPhi)) h_dPhi_continu_range->Fill(dPhiCutLess);
                      }
                    }
                    if (pi0_optimize)
                    {
                      h_pi0_continu->Fill(pi0Mass_continu);
                      for (pi0MassCut_range=0.0005; pi0MassCut_range<0.1995; pi0MassCut_range+=0.001)
                      {
                        if (pi0MassCut(pi0Mass_continu)) h_pi0_continu_range->Fill(pi0MassCut_range);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  // Create Significance and Precision plots
  
  setValues();
  
  double signalScale=(luminosity*prodCrossSection_DsDss*branchingFr_signal*branchingFr_mode)/nSignalSample;
  std::cout<<"Generated signal sample events = "<<nSignalSample<<std::endl;
  std::cout<<"Expected number of signal events in this mode is "<<luminosity<<"(luminosity)*";
  std::cout<<prodCrossSection_DsDss<<"(prodCrossSection_DsDss)*";
  std::cout<<branchingFr_signal<<"(branchingFr_signal)*";
  std::cout<<branchingFr_mode<<"(branchingFr_mode)=";  
  std::cout<<(luminosity*prodCrossSection_DsDss*branchingFr_signal*branchingFr_mode)<<std::endl;
  std::cout<<"Signal Scale = "<<signalScale<<std::endl;
  
  // Scale plots appropriately
  h_dsPlusM_signal->Scale(signalScale);
  h_dsPlusM_generic->Scale(genericScale);
  h_dsPlusM_continu->Scale(continuScale);
  h_dsPlusM_signal_range->Scale(signalScale);
  h_dsPlusM_generic_range->Scale(genericScale);
  h_dsPlusM_continu_range->Scale(continuScale);
  h_MBC_signal->Scale(signalScale);
  h_MBC_generic->Scale(genericScale);
  h_MBC_continu->Scale(continuScale);
  h_MBC_signal_range->Scale(signalScale);
  h_MBC_generic_range->Scale(genericScale);
  h_MBC_continu_range->Scale(continuScale);
  h_DeltaM_signal->Scale(signalScale);
  h_DeltaM_generic->Scale(genericScale);
  h_DeltaM_continu->Scale(continuScale);
  h_DeltaM_signal_range->Scale(signalScale);
  h_DeltaM_generic_range->Scale(genericScale);
  h_DeltaM_continu_range->Scale(continuScale);
  h_diffD0_signal->Scale(signalScale);
  h_diffD0_generic->Scale(genericScale);
  h_diffD0_continu->Scale(continuScale);
  h_diffD0_signal_range->Scale(signalScale);
  h_diffD0_generic_range->Scale(genericScale);
  h_diffD0_continu_range->Scale(continuScale);
  h_dPhi_signal->Scale(signalScale);
  h_dPhi_generic->Scale(genericScale);
  h_dPhi_continu->Scale(continuScale);
  h_dPhi_signal_range->Scale(signalScale);
  h_dPhi_generic_range->Scale(genericScale);
  h_dPhi_continu_range->Scale(continuScale);
  h_pi0_signal->Scale(signalScale);
  h_pi0_generic->Scale(genericScale);
  h_pi0_continu->Scale(continuScale);
  h_pi0_signal_range->Scale(signalScale);
  h_pi0_generic_range->Scale(genericScale);
  h_pi0_continu_range->Scale(continuScale);
  
  TLine *line;
  float xmin, xmax, ymin, ymax;
  gROOT->SetStyle("Plain");
  
  if (DsMass_optimize)
  {
    for (int i=1; i<=h_dsPlusM_significance->GetNbinsX(); ++i)
    {
      double denom=pow(h_dsPlusM_generic_range->GetBinContent(i)+h_dsPlusM_continu_range->GetBinContent(i), 0.5);
      if (denom>0) h_dsPlusM_significance->SetBinContent(i, h_dsPlusM_signal_range->GetBinContent(i)/denom);
      else h_dsPlusM_significance->SetBinContent(i, 0);
      
      double denom_precision=pow(h_dsPlusM_signal_range->GetBinContent(i)+h_dsPlusM_generic_range->GetBinContent(i)+h_dsPlusM_continu_range->GetBinContent(i), 0.5);
      if (denom_precision>0) h_dsPlusM_precision->SetBinContent(i, h_dsPlusM_signal_range->GetBinContent(i)/denom_precision);
      else h_dsPlusM_precision->SetBinContent(i,0);
    }
    TCanvas *c_DsPlusM = new TCanvas("c_DsPlusM", "", 600, 1200);
    xmin=dsPlusMCut_center-dsPlusMCut_range;
    xmax=dsPlusMCut_center+dsPlusMCut_range;
    c_DsPlusM->Divide(2,4);
    c_DsPlusM->cd(1);
    h_dsPlusM_signal->Draw();
    ymax=(h_dsPlusM_signal->GetMaximum())*0.95;
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    c_DsPlusM->cd(3);
    h_dsPlusM_generic->Draw();
    ymax=(h_dsPlusM_generic->GetMaximum())*0.95;
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    c_DsPlusM->cd(5);
    h_dsPlusM_continu->Draw();
    ymax=(h_dsPlusM_continu->GetMaximum())*0.95;
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    c_DsPlusM->cd(2);
    h_dsPlusM_signal_range->Draw();
    ymax=(h_dsPlusM_signal_range->GetBinContent(h_dsPlusM_signal_range->FindBin(dsPlusMCut_range)));
    line = new TLine(dsPlusMCut_range,ymax,dsPlusMCut_range,h_dsPlusM_signal_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(dsPlusMCut_range,ymax,0,ymax); line->Draw();
    c_DsPlusM->cd(4);
    h_dsPlusM_generic_range->Draw();
    ymax=(h_dsPlusM_generic_range->GetBinContent(h_dsPlusM_generic_range->FindBin(dsPlusMCut_range)));
    line = new TLine(dsPlusMCut_range,ymax,dsPlusMCut_range,h_dsPlusM_generic_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(dsPlusMCut_range,ymax,0,ymax); line->Draw();
    c_DsPlusM->cd(6);
    h_dsPlusM_continu_range->Draw();
    ymax=(h_dsPlusM_continu_range->GetBinContent(h_dsPlusM_continu_range->FindBin(dsPlusMCut_range)));
    line = new TLine(dsPlusMCut_range,ymax,dsPlusMCut_range,h_dsPlusM_continu_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(dsPlusMCut_range,ymax,0,ymax); line->Draw();
    c_DsPlusM->cd(7);
    h_dsPlusM_significance->Draw();
    ymax=(h_dsPlusM_significance->GetBinContent(h_dsPlusM_significance->FindBin(dsPlusMCut_range)));
    line = new TLine(dsPlusMCut_range,ymax,dsPlusMCut_range,h_dsPlusM_significance->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(dsPlusMCut_range,ymax,0,ymax); line->Draw();
    c_DsPlusM->cd(8);
    h_dsPlusM_precision->Draw();
    ymax=(h_dsPlusM_precision->GetBinContent(h_dsPlusM_precision->FindBin(dsPlusMCut_range)));
    line = new TLine(dsPlusMCut_range,ymax,dsPlusMCut_range,h_dsPlusM_precision->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(dsPlusMCut_range,ymax,0,ymax); line->Draw();
    std::string filename=decay+"_Optimize_DsPlusM_pionFit.eps";
    c_DsPlusM->SaveAs(filename.c_str());
  }
  
  if (mBC_optimize)
  {
    for (int i=1; i<h_MBC_significance->GetNbinsX(); ++i)
    {
      double denom=pow(h_MBC_generic_range->GetBinContent(i)+h_MBC_continu_range->GetBinContent(i), 0.5);
      if (denom>0) h_MBC_significance->SetBinContent(i, h_MBC_signal_range->GetBinContent(i)/denom);
      else h_MBC_significance->SetBinContent(i, 0);
      
      double denom_precision=pow(h_MBC_signal_range->GetBinContent(i)+h_MBC_generic_range->GetBinContent(i)+h_MBC_continu_range->GetBinContent(i), 0.5);
      if (denom_precision>0) h_MBC_precision->SetBinContent(i, h_MBC_signal_range->GetBinContent(i)/denom_precision);
      else h_MBC_precision->SetBinContent(i,0);
    }
    TCanvas *c_MBC = new TCanvas("c_MBC", "", 600, 1200);
    xmin=mbcCut_center-mbcCut_range; xmax=mbcCut_center+mbcCut_range;
    c_MBC->Divide(2,4);
    c_MBC->cd(1);
    h_MBC_signal->Draw();
    ymax=(h_MBC_signal->GetMaximum())*0.95;
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    c_MBC->cd(3);
    h_MBC_generic->Draw();
    ymax=(h_MBC_generic->GetMaximum())*0.95;
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    c_MBC->cd(5);
    h_MBC_continu->Draw();
    ymax=(h_MBC_continu->GetMaximum())*0.95;
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    c_MBC->cd(2);
    h_MBC_signal_range->Draw();
    ymax=(h_MBC_signal_range->GetBinContent(h_MBC_signal_range->FindBin(mbcCut_range)));
    line = new TLine(mbcCut_range,ymax,mbcCut_range,h_MBC_signal_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(mbcCut_range,ymax,0,ymax); line->Draw();
    c_MBC->cd(4);
    h_MBC_generic_range->Draw();
    ymax=(h_MBC_generic_range->GetBinContent(h_MBC_generic_range->FindBin(mbcCut_range)));
    line = new TLine(mbcCut_range,ymax,mbcCut_range,h_MBC_generic_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(mbcCut_range,ymax,0,ymax); line->Draw();
    c_MBC->cd(6);
    h_MBC_continu_range->Draw();
    ymax=(h_MBC_continu_range->GetBinContent(h_MBC_continu_range->FindBin(mbcCut_range)));
    line = new TLine(mbcCut_range,ymax,mbcCut_range,h_MBC_continu_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(mbcCut_range,ymax,0,ymax); line->Draw();
    c_MBC->cd(7);
    h_MBC_significance->Draw();
    ymax=(h_MBC_significance->GetBinContent(h_MBC_significance->FindBin(mbcCut_range)));
    line = new TLine(mbcCut_range,ymax,mbcCut_range,h_MBC_significance->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(mbcCut_range,ymax,0,ymax); line->Draw();
    c_MBC->cd(8);
    h_MBC_precision->Draw();
    ymax=(h_MBC_precision->GetBinContent(h_MBC_precision->FindBin(mbcCut_range)));
    line = new TLine(mbcCut_range,ymax,mbcCut_range,h_MBC_precision->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(mbcCut_range,ymax,0,ymax); line->Draw();
    c_MBC->SaveAs((decay+"_Optimize_MBC_pionFit.eps").c_str());
  }
  
  if (DeltaM_optimize)
  {
    for (int i=1; i<h_DeltaM_significance->GetNbinsX(); ++i)
    {
      double denom=pow(h_DeltaM_generic_range->GetBinContent(i)+h_DeltaM_continu_range->GetBinContent(i), 0.5);
      if (denom>0) h_DeltaM_significance->SetBinContent(i, h_DeltaM_signal_range->GetBinContent(i)/denom);
      else h_DeltaM_significance->SetBinContent(i, 0);
      
      double denom_precision=pow(h_DeltaM_signal_range->GetBinContent(i)+h_DeltaM_generic_range->GetBinContent(i)+h_DeltaM_continu_range->GetBinContent(i), 0.5);
      if (denom_precision>0) h_DeltaM_precision->SetBinContent(i, h_DeltaM_signal_range->GetBinContent(i)/denom_precision);
      else h_DeltaM_precision->SetBinContent(i,0);
    }
    TCanvas *c_DeltaM = new TCanvas("c_DeltaM", "", 600, 1200);
    xmin=deltaMCut_center-deltaMCut_range; xmax=deltaMCut_center+deltaMCut_range;
    c_DeltaM->Divide(2,4);
    c_DeltaM->cd(1);
    h_DeltaM_signal->Draw();
    ymax=(h_DeltaM_signal->GetMaximum())*0.95;
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    c_DeltaM->cd(3);
    h_DeltaM_generic->Draw();
    ymax=(h_DeltaM_generic->GetMaximum())*0.95;
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    c_DeltaM->cd(5);
    h_DeltaM_continu->Draw();
    ymax=(h_DeltaM_continu->GetMaximum())*0.95;
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    c_DeltaM->cd(2);
    h_DeltaM_signal_range->Draw();
    ymax=(h_DeltaM_signal_range->GetBinContent(h_DeltaM_signal_range->FindBin(deltaMCut_range)));
    line = new TLine(deltaMCut_range,ymax,deltaMCut_range,h_DeltaM_signal_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(deltaMCut_range,ymax,0,ymax); line->Draw();
    c_DeltaM->cd(4);
    h_DeltaM_generic_range->Draw();
    ymax=(h_DeltaM_generic_range->GetBinContent(h_DeltaM_generic_range->FindBin(deltaMCut_range)));
    line = new TLine(deltaMCut_range,ymax,deltaMCut_range,h_DeltaM_generic_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(deltaMCut_range,ymax,0,ymax); line->Draw();
    c_DeltaM->cd(6);
    h_DeltaM_continu_range->Draw();
    ymax=(h_DeltaM_continu_range->GetBinContent(h_DeltaM_continu_range->FindBin(deltaMCut_range)));
    line = new TLine(deltaMCut_range,ymax,deltaMCut_range,h_DeltaM_continu_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(deltaMCut_range,ymax,0,ymax); line->Draw();
    c_DeltaM->cd(7);
    h_DeltaM_significance->Draw();
    ymax=(h_DeltaM_significance->GetBinContent(h_DeltaM_significance->FindBin(deltaMCut_range)));
    line = new TLine(deltaMCut_range,ymax,deltaMCut_range,h_DeltaM_significance->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(deltaMCut_range,ymax,0,ymax); line->Draw();
    c_DeltaM->cd(8);
    h_DeltaM_precision->Draw();
    ymax=(h_DeltaM_precision->GetBinContent(h_DeltaM_precision->FindBin(deltaMCut_range)));
    line = new TLine(deltaMCut_range,ymax,deltaMCut_range,h_DeltaM_precision->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(deltaMCut_range,ymax,0,ymax); line->Draw();
    c_DeltaM->SaveAs((decay+"_Optimize_DeltaM_pionFit.eps").c_str());
  }
  
  if (diffD0_optimize)
  {
    for (int i=1; i<h_diffD0_significance->GetNbinsX(); ++i)
    {
      double denom=pow(h_diffD0_generic_range->GetBinContent(i)+h_diffD0_continu_range->GetBinContent(i), 0.5);
      if (denom>0) h_diffD0_significance->SetBinContent(i, h_diffD0_signal_range->GetBinContent(i)/denom);
      else h_diffD0_significance->SetBinContent(i, 0);
      
      double denom_precision=pow(h_diffD0_signal_range->GetBinContent(i)+h_diffD0_generic_range->GetBinContent(i)+h_diffD0_continu_range->GetBinContent(i), 0.5);
      if (denom_precision>0) h_diffD0_precision->SetBinContent(i, h_diffD0_signal_range->GetBinContent(i)/denom_precision);
      else h_diffD0_precision->SetBinContent(i,0);
    }
    TCanvas *c_diffD0 = new TCanvas("c_diffD0", "", 600, 1200);
    c_diffD0->Divide(2,4);
    c_diffD0->cd(1);
    h_diffD0_signal->Draw();
    ymax=(h_diffD0_signal->GetMaximum())*0.95;
    line = new TLine(diffD0Cut,ymax,diffD0Cut,0); line->Draw();
    c_diffD0->cd(3);
    h_diffD0_generic->Draw();
    ymax=(h_diffD0_generic->GetMaximum())*0.95;
    line = new TLine(diffD0Cut,ymax,diffD0Cut,0); line->Draw();
    c_diffD0->cd(5);
    h_diffD0_continu->Draw();
    ymax=(h_diffD0_continu->GetMaximum())*0.95;
    line = new TLine(diffD0Cut,ymax,diffD0Cut,0); line->Draw();
    c_diffD0->cd(2);
    h_diffD0_signal_range->Draw();
    ymax=(h_diffD0_signal_range->GetBinContent(h_diffD0_signal_range->FindBin(diffD0Cut)));
    line = new TLine(diffD0Cut,ymax,diffD0Cut,h_diffD0_signal_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(diffD0Cut,ymax,h_diffD0_signal_range->GetXaxis()->GetXmin(),ymax); line->Draw();
    c_diffD0->cd(4);
    h_diffD0_generic_range->Draw();
    ymax=(h_diffD0_generic_range->GetBinContent(h_diffD0_generic_range->FindBin(diffD0Cut)));
    line = new TLine(diffD0Cut,ymax,diffD0Cut,h_diffD0_generic_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(diffD0Cut,ymax,h_diffD0_generic_range->GetXaxis()->GetXmin(),ymax); line->Draw();
    c_diffD0->cd(6);
    h_diffD0_continu_range->Draw();
    ymax=(h_diffD0_continu_range->GetBinContent(h_diffD0_continu_range->FindBin(diffD0Cut)));
    line = new TLine(diffD0Cut,ymax,diffD0Cut,h_diffD0_continu_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(diffD0Cut,ymax,h_diffD0_continu_range->GetXaxis()->GetXmin(),ymax); line->Draw();
    c_diffD0->cd(7);
    h_diffD0_significance->Draw();
    ymax=(h_diffD0_significance->GetBinContent(h_diffD0_significance->FindBin(diffD0Cut)));
    line = new TLine(diffD0Cut,ymax,diffD0Cut,h_diffD0_significance->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(diffD0Cut,ymax,h_diffD0_significance->GetXaxis()->GetXmin(),ymax); line->Draw();
    c_diffD0->cd(8);
    h_diffD0_precision->Draw();
    ymax=(h_diffD0_precision->GetBinContent(h_diffD0_precision->FindBin(diffD0Cut)));
    line = new TLine(diffD0Cut,ymax,diffD0Cut,h_diffD0_precision->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(diffD0Cut,ymax,h_diffD0_precision->GetXaxis()->GetXmin(),ymax); line->Draw();
    c_diffD0->SaveAs((decay+"_Optimize_diffD0_pionFit.eps").c_str());
  }
  
  if (dPhi_optimize)
  {
    for (int i=1; i<=h_dPhi_significance->GetNbinsX(); ++i)
    {
      double denom=pow(h_dPhi_generic_range->GetBinContent(i)+h_dPhi_continu_range->GetBinContent(i), 0.5);
      if (denom>0) h_dPhi_significance->SetBinContent(i, h_dPhi_signal_range->GetBinContent(i)/denom);
      else h_dPhi_significance->SetBinContent(i, 0);
      
      double denom_precision=pow(h_dPhi_signal_range->GetBinContent(i)+h_dPhi_generic_range->GetBinContent(i)+h_dPhi_continu_range->GetBinContent(i), 0.5);
      if (denom_precision>0) h_dPhi_precision->SetBinContent(i, h_dPhi_signal_range->GetBinContent(i)/denom_precision);
      else h_dPhi_precision->SetBinContent(i,0);
    }
    TCanvas *c_dPhi = new TCanvas("c_dPhi", "", 600, 1200);
    c_dPhi->Divide(2,4);
    c_dPhi->cd(1);
    h_dPhi_signal->Draw();
    ymax=(h_dPhi_signal->GetMaximum())*0.95;
    line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
    c_dPhi->cd(3);
    h_dPhi_generic->Draw();
    ymax=(h_dPhi_generic->GetMaximum())*0.95;
    line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
    c_dPhi->cd(5);
    h_dPhi_continu->Draw();
    ymax=(h_dPhi_continu->GetMaximum())*0.95;
    line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
    c_dPhi->cd(2);
    h_dPhi_signal_range->Draw();
    ymax=(h_dPhi_signal_range->GetBinContent(h_dPhi_signal_range->FindBin(dPhiCutLess)));
    line = new TLine(dPhiCutLess,ymax,dPhiCutLess,h_dPhi_signal_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(dPhiCutLess,ymax,h_dPhi_signal_range->GetXaxis()->GetXmin(),ymax); line->Draw();
    c_dPhi->cd(4);
    h_dPhi_generic_range->Draw();
    ymax=(h_dPhi_generic_range->GetBinContent(h_dPhi_generic_range->FindBin(dPhiCutLess)));
    line = new TLine(dPhiCutLess,ymax,dPhiCutLess,h_dPhi_generic_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(dPhiCutLess,ymax,h_dPhi_generic_range->GetXaxis()->GetXmin(),ymax); line->Draw();
    c_dPhi->cd(6);
    h_dPhi_continu_range->Draw();
    ymax=(h_dPhi_continu_range->GetBinContent(h_dPhi_continu_range->FindBin(dPhiCutLess)));
    line = new TLine(dPhiCutLess,ymax,dPhiCutLess,h_dPhi_continu_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(dPhiCutLess,ymax,h_dPhi_continu_range->GetXaxis()->GetXmin(),ymax); line->Draw();
    c_dPhi->cd(7);
    h_dPhi_significance->Draw();
    ymax=(h_dPhi_significance->GetBinContent(h_dPhi_significance->FindBin(dPhiCutLess)));
    line = new TLine(dPhiCutLess,ymax,dPhiCutLess,h_dPhi_significance->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(dPhiCutLess,ymax,h_dPhi_continu_range->GetXaxis()->GetXmin(),ymax); line->Draw();
    c_dPhi->cd(8);
    h_dPhi_precision->Draw();
    ymax=(h_dPhi_precision->GetBinContent(h_dPhi_precision->FindBin(dPhiCutLess)));
    line = new TLine(dPhiCutLess,ymax,dPhiCutLess,h_dPhi_precision->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(dPhiCutLess,ymax,h_dPhi_precision->GetXaxis()->GetXmin(),ymax); line->Draw();
    c_dPhi->SaveAs((decay+"_Optimize_dPhi_pionFit.eps").c_str());
  }
  
  if (pi0_optimize)
  {
    for (int i=1; i<=h_pi0_significance->GetNbinsX(); ++i)
    {
      double denom=pow(h_pi0_generic_range->GetBinContent(i)+h_pi0_continu_range->GetBinContent(i), 0.5);
      if (denom>0) h_pi0_significance->SetBinContent(i, h_pi0_signal_range->GetBinContent(i)/denom);
      else h_pi0_significance->SetBinContent(i, 0);
      
      double denom_precision=pow(h_pi0_signal_range->GetBinContent(i)+h_pi0_generic_range->GetBinContent(i)+h_pi0_continu_range->GetBinContent(i), 0.5);
      if (denom_precision>0) h_pi0_precision->SetBinContent(i, h_pi0_signal_range->GetBinContent(i)/denom_precision);
      else h_pi0_precision->SetBinContent(i,0);
    }
    TCanvas *c_pi0 = new TCanvas("c_pi0", "", 600, 1200);
    xmin=pi0MassCut_center-pi0MassCut_range; xmax=pi0MassCut_center+pi0MassCut_range;
    c_pi0->Divide(2,4);
    c_pi0->cd(1);
    h_pi0_signal->Draw();
    ymax=(h_pi0_signal->GetMaximum())*0.95;
    line = new TLine(xmin, ymax, xmin, 0); line->Draw();
    line = new TLine(xmax, ymax, xmax, 0); line->Draw();
    c_pi0->cd(3);
    h_pi0_generic->Draw();
    ymax=(h_pi0_generic->GetMaximum())*0.95;
    line = new TLine(xmin, ymax, xmin, 0); line->Draw();
    line = new TLine(xmax, ymax, xmax, 0); line->Draw();
    c_pi0->cd(5);
    h_pi0_continu->Draw();
    ymax=(h_pi0_continu->GetMaximum())*0.95;
    line = new TLine(xmin, ymax, xmin, 0); line->Draw();
    line = new TLine(xmax, ymax, xmax, 0); line->Draw();
    c_pi0->cd(2);
    h_pi0_signal_range->Draw();
    ymax=(h_pi0_signal_range->GetBinContent(h_pi0_signal_range->FindBin(pi0MassCut_range)));
    line = new TLine(pi0MassCut_range,ymax,pi0MassCut_range,h_pi0_signal_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(pi0MassCut_range,ymax,0,ymax); line->Draw();
    c_pi0->cd(4);
    h_pi0_generic_range->Draw();
    ymax=(h_pi0_generic_range->GetBinContent(h_pi0_generic_range->FindBin(pi0MassCut_range)));
    line = new TLine(pi0MassCut_range,ymax,pi0MassCut_range,h_pi0_generic_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(pi0MassCut_range,ymax,0,ymax); line->Draw();
    c_pi0->cd(6);
    h_pi0_continu_range->Draw();
    ymax=(h_pi0_continu_range->GetBinContent(h_pi0_continu_range->FindBin(pi0MassCut_range)));
    line = new TLine(pi0MassCut_range,ymax,pi0MassCut_range,h_pi0_continu_range->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(pi0MassCut_range,ymax,0,ymax); line->Draw();
    c_pi0->cd(7);
    h_pi0_significance->Draw();
    ymax=(h_pi0_significance->GetBinContent(h_pi0_significance->FindBin(pi0MassCut_range)));
    line = new TLine(pi0MassCut_range,ymax,pi0MassCut_range,h_pi0_significance->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(pi0MassCut_range,ymax,0,ymax); line->Draw();
    c_pi0->cd(8);
    h_pi0_precision->Draw();
    ymax=(h_pi0_precision->GetBinContent(h_pi0_precision->FindBin(pi0MassCut_range)));
    line = new TLine(pi0MassCut_range,ymax,pi0MassCut_range,h_pi0_precision->GetYaxis()->GetXmin()); line->Draw();
    line = new TLine(pi0MassCut_range,ymax,0,ymax); line->Draw();
    c_pi0->SaveAs((decay+"_Optimize_pi0_pionFit.eps").c_str());
  }
  
  if ( !DsMass_optimize &&
       !mBC_optimize    &&
       !DeltaM_optimize &&
       !diffD0_optimize &&
       !dPhi_optimize   &&
       !pi0_optimize)
  {
    std::cout<<"Trying to stack"<<std::endl;
    h_dsPlusM_signal->SetFillColor(kRed);   
    h_dsPlusM_continu->SetFillColor(kGreen);
    h_dsPlusM_generic->SetFillColor(kBlue); 
    h_MBC_signal->SetFillColor(kRed);       
    h_MBC_continu->SetFillColor(kGreen);    
    h_MBC_generic->SetFillColor(kBlue);     
    h_DeltaM_signal->SetFillColor(kRed);    
    h_DeltaM_continu->SetFillColor(kGreen); 
    h_DeltaM_generic->SetFillColor(kBlue);  
    h_diffD0_signal->SetFillColor(kRed);    
    h_diffD0_continu->SetFillColor(kGreen); 
    h_diffD0_generic->SetFillColor(kBlue);  
    h_dPhi_signal->SetFillColor(kRed);      
    h_dPhi_continu->SetFillColor(kGreen);   
    h_dPhi_generic->SetFillColor(kBlue);    
    h_pi0_signal->SetFillColor(kRed);       
    h_pi0_continu->SetFillColor(kGreen);    
    h_pi0_generic->SetFillColor(kBlue);     
  
    TCanvas *c_Stacks = new TCanvas("c_Stacks", decay.c_str(), 800, 1200); //, "Stacked Signal Region"
    c_Stacks->Divide(2,3);
    c_Stacks->cd(1);
    THStack *s_DsPlusM_stacked = new THStack("s_DsPlusM_stacked", "DsPlusM");
    s_DsPlusM_stacked->Add(h_dsPlusM_continu);
    s_DsPlusM_stacked->Add(h_dsPlusM_generic);
    s_DsPlusM_stacked->Add(h_dsPlusM_signal);
    s_DsPlusM_stacked->Draw();
    xmin=dsPlusMCut_center-dsPlusMCut_range;
    xmax=dsPlusMCut_center+dsPlusMCut_range;
    ymax=(s_DsPlusM_stacked->GetMaximum())*0.95;
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    c_Stacks->cd(2);
    THStack *s_MBC_stacked = new THStack("s_MBC_stacked", "MBC");
    s_MBC_stacked->Add(h_MBC_continu);
    s_MBC_stacked->Add(h_MBC_generic);
    s_MBC_stacked->Add(h_MBC_signal);
    s_MBC_stacked->Draw();
    xmin=mbcCut_center-mbcCut_range;
    xmax=mbcCut_center+mbcCut_range;
    ymax=(s_MBC_stacked->GetMaximum())*0.95;
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    c_Stacks->cd(3);
    THStack *s_DeltaM_stacked = new THStack("s_DeltaM_stacked", "DeltaM");
    s_DeltaM_stacked->Add(h_DeltaM_continu);
    s_DeltaM_stacked->Add(h_DeltaM_generic);
    s_DeltaM_stacked->Add(h_DeltaM_signal);
    s_DeltaM_stacked->Draw();
    xmin=deltaMCut_center-deltaMCut_range;
    xmax=deltaMCut_center+deltaMCut_range;
    ymax=(s_DeltaM_stacked->GetMaximum())*0.95;
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    c_Stacks->cd(4);
    THStack *s_diffD0_stacked = new THStack("s_diffD0_stacked", "diffD0");
    s_diffD0_stacked->Add(h_diffD0_continu);
    s_diffD0_stacked->Add(h_diffD0_generic);
    s_diffD0_stacked->Add(h_diffD0_signal);
    //s_diffD0_stacked->GetXaxis()->SetLabelSize(0.02);
    s_diffD0_stacked->Draw();
    ymax=(s_diffD0_stacked->GetMaximum())*0.95;
    line = new TLine(diffD0Cut,ymax,diffD0Cut,0); line->Draw();
    c_Stacks->cd(5);
    THStack *s_dPhi_stacked = new THStack("s_dPhi_stacked", "dPhi");
    s_dPhi_stacked->Add(h_dPhi_continu);
    s_dPhi_stacked->Add(h_dPhi_generic);
    s_dPhi_stacked->Add(h_dPhi_signal);
    s_dPhi_stacked->Draw();
    ymax=(s_dPhi_stacked->GetMaximum())*0.95;
    line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
    std::cout<<"dPhiCutLess = "<<dPhiCutLess<<std::endl;
    c_Stacks->cd(6);
    THStack *s_pi0_stacked = new THStack("s_pi0_stacked", "pi0");
    s_pi0_stacked->Add(h_pi0_continu);
    s_pi0_stacked->Add(h_pi0_generic);
    s_pi0_stacked->Add(h_pi0_signal);
    s_pi0_stacked->Draw();
    xmin=pi0MassCut_center-pi0MassCut_range;
    xmax=pi0MassCut_center+pi0MassCut_range;
    ymax=(s_pi0_stacked->GetMaximum())*0.95;
    line = new TLine(xmin, ymax, xmin, 0); line->Draw();
    line = new TLine(xmax, ymax, xmax, 0); line->Draw();
    c_Stacks->SaveAs((decay+"_Optimize_signalRegion_pionFit.eps").c_str());
    
    std::cout<<" === Spill the Numbers in "<<decay<<" Channel = pion-fit =="<<std::endl;
    std::cout<<" Expected # of Signal Events = "<<signalEvents.pi0MassCut*signalScale<<std::endl;
    std::cout<<" Expected # of Generic MC Events = "<<genericEvents.pi0MassCut*genericScale<<std::endl;
    std::cout<<" Expected # of Continuum MC Events = "<<continuEvents.pi0MassCut*continuScale<<std::endl;
    double totalBackground=genericEvents.pi0MassCut*genericScale+continuEvents.pi0MassCut*continuScale;
    std::cout<<" Expected Total Background Events = "<<totalBackground<<std::endl;
    std::cout<<" Significance = "<<signalEvents.pi0MassCut*signalScale/pow(totalBackground, 0.5)<<std::endl;
    std::cout<<" === "<<std::endl;
    
    gStyle->SetOptStat(0);
    
    TCanvas *electronEResolution = new TCanvas("electronEResolution", "electronEResolution", 500, 500);
    gPad->SetLeftMargin(0.18);
    h_electronE->GetYaxis()->SetTitleOffset(2.3);
    h_electronE->Draw("box");
    
    TCanvas *electronEResolution_p = new TCanvas("electronEResolution_p", "electronEResolution_p", 500, 500);
    TProfile *p_electronE=h_electronE->ProfileX();
    // TF1 *f1=new TF1("f1", "[0]+[1]/(x*x)", 0.05, 0.13);
    // p_electronE->Fit("f1", "R");
    gPad->SetLeftMargin(0.18);
    p_electronE->GetYaxis()->SetTitleOffset(2.3);
    p_electronE->GetYaxis()->SetTitle("Reco - MC Electron Energy (GeV)");
    p_electronE->Draw();
    
  }
  
  return 0;
}

/*
int DsTaggedAnalysis_Optimize_pionFit()
{
  DsTaggedAnalysis_Optimize_pionFit_Criterion(true, false, false, false, false);
  DsTaggedAnalysis_Optimize_pionFit_Criterion(false, true, false, false, false);
  DsTaggedAnalysis_Optimize_pionFit_Criterion(false, false, true, false, false);
  DsTaggedAnalysis_Optimize_pionFit_Criterion(false, false, false, true, false);
  DsTaggedAnalysis_Optimize_pionFit_Criterion(false, false, false, false, true);
  
  return 0;
}
*/
