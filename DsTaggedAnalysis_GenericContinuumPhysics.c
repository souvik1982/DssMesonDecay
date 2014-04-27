#include <iostream>
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TFile.h"
#include "TChain.h"
#include "TPaveStats.h"
#include <set>
#include <map>
#include <iomanip>

bool DeltaM_sideband=false;
bool mBC_sideband=false;
bool DsMass_sideband=true;

std::string decay="pietaprimerho";
// options: KKpi, KsK, pieta, pietaprime, KKpipi0, pipipi, KsKmpipi, pipi0eta, pietaprimerho

double pi=3.14159265358979;
float decayNumber;
double dsPlusMCut_center, dsPlusMCut_range;
double mbcCut_center, mbcCut_range;
double deltaMCut_center=0.150, deltaMCut_range=0.013;
double diffD0Cut;
double dPhiCutLess;
double pi0MassCut_center, pi0MassCut_range;
double deltaECut_center=0.012, deltaECut_range=0.019;
double electronEnergyThreshold=0.15;


std::string dsPlusM_generic_fileName="dsPlusM_generic", DeltaE_generic_fileName="DeltaE_generic", MBC_generic_fileName="mBC_generic", DeltaM_generic_fileName="DeltaM_generic", diffD0_generic_fileName="diffD0_generic", dPhi_generic_fileName="dPhi_generic", dPhi_diffD0_generic_fileName="dPhi_diffD0_generic", pi0_generic_fileName="pi0_generic_";
std::string dsPlusM_continu_fileName="dsPlusM_continu", DeltaE_continu_fileName="DeltaE_continu", MBC_continu_fileName="mBC_continu", DeltaM_continu_fileName="DeltaM_continu", diffD0_continu_fileName="diffD0_continu", dPhi_continu_fileName="dPhi_continu", dPhi_diffD0_continu_fileName="dPhi_diffD0_continu", pi0_continu_fileName="pi0_continu_";
std::string dsPlusM_physics1_fileName="dsPlusM_physics1", DeltaE_physics1_fileName="DeltaE_physics1", MBC_physics1_fileName="mBC_physics1", DeltaM_physics1_fileName="DeltaM_physics1", diffD0_physics1_fileName="diffD0_physics1", dPhi_physics1_fileName="dPhi_physics1", dPhi_diffD0_physics1_fileName="dPhi_diffD0_physics1", pi0_physics1_fileName="pi0_physics1_";
std::string dsPlusM_physics3_fileName="dsPlusM_physics3", DeltaE_physics3_fileName="DeltaE_physics3", MBC_physics3_fileName="mBC_physics3", DeltaM_physics3_fileName="DeltaM_physics3", diffD0_physics3_fileName="diffD0_physics3", dPhi_physics3_fileName="dPhi_physics3", dPhi_diffD0_physics3_fileName="dPhi_diffD0_physics3", pi0_physics3_fileName="pi0_physics3_";

void setValues()
{
  if (decay=="KKpi")
  {
    decayNumber=401;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.011;
    mbcCut_center=2.112; mbcCut_range=0.005;
    // deltaMCut_center=0.155; deltaMCut_range=0.008;
    diffD0Cut=-0.002;
    dPhiCutLess=0.01;
    pi0MassCut_center=0.135, pi0MassCut_range=0.0;
  }
  else if (decay=="KsK")
  {
    decayNumber=400;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.012;
    mbcCut_center=2.112; mbcCut_range=0.006;
    // deltaMCut_center=0.158; deltaMCut_range=0.01;
    diffD0Cut=-0.002;
    dPhiCutLess=0.11;
    pi0MassCut_center=0.135, pi0MassCut_range=0.008;
  }
  else if (decay=="pieta")
  {
    decayNumber=440;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.015;
    mbcCut_center=2.112; mbcCut_range=0.009;
    // deltaMCut_center=0.158; deltaMCut_range=0.008;
    diffD0Cut=-0.007;
    dPhiCutLess=0.07;
    pi0MassCut_center=0.135, pi0MassCut_range=0.005;
  }
  else if (decay=="pietaprime")
  {
    decayNumber=460;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.011;
    mbcCut_center=2.112; mbcCut_range=0.016;
    // deltaMCut_center=0.158; deltaMCut_range=0.013;
    diffD0Cut=-0.003;
    dPhiCutLess=0.07;
    pi0MassCut_center=0.135, pi0MassCut_range=0.012;
  }
  else if (decay=="KKpipi0")
  {
    decayNumber=404;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.009;
    mbcCut_center=2.112; mbcCut_range=0.007;
    // deltaMCut_center=0.158; deltaMCut_range=0.011;
    diffD0Cut=-0.002;
    dPhiCutLess=0.07;
    pi0MassCut_center=0.135, pi0MassCut_range=0.009;
  }
  else if (decay=="pipipi")
  {
    decayNumber=421;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.013;
    mbcCut_center=2.112; mbcCut_range=0.005;
    // deltaMCut_center=0.158; deltaMCut_range=0.009;
    diffD0Cut=-0.001;
    dPhiCutLess=0.06;
    pi0MassCut_center=0.135, pi0MassCut_range=0.0;
  }
  else if (decay=="KsKmpipi")
  {
    decayNumber=406;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.007;
    mbcCut_center=2.112; mbcCut_range=0.007;
    // deltaMCut_center=0.158; deltaMCut_range=0.008;
    diffD0Cut=-0.005;
    dPhiCutLess=0.07;
    pi0MassCut_center=0.135, pi0MassCut_range=0.011;
  }
  else if (decay=="pipi0eta")
  {
    decayNumber=441;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.014;
    mbcCut_center=2.112; mbcCut_range=0.01;
    // deltaMCut_center=0.158; deltaMCut_range=0.009;
    diffD0Cut=-0.003;
    dPhiCutLess=0.06;
    pi0MassCut_center=0.135, pi0MassCut_range=0.007;
  }
  else if (decay=="pietaprimerho")
  {
    decayNumber=480;
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.017;
    mbcCut_center=2.112; mbcCut_range=0.003;
    // deltaMCut_center=0.155; deltaMCut_range=0.008;
    diffD0Cut=-0.004;
    dPhiCutLess=0.09;
    pi0MassCut_center=0.135, pi0MassCut_range=0.015;
  }
  
  dsPlusM_generic_fileName+=decay; DeltaE_generic_fileName+=decay; MBC_generic_fileName+=decay; DeltaM_generic_fileName+=decay; diffD0_generic_fileName+=decay; dPhi_generic_fileName+=decay; dPhi_diffD0_generic_fileName+=decay; pi0_generic_fileName+=decay;
  dsPlusM_continu_fileName+=decay; DeltaE_continu_fileName+=decay; MBC_continu_fileName+=decay; DeltaM_continu_fileName+=decay; diffD0_continu_fileName+=decay; dPhi_continu_fileName+=decay; dPhi_diffD0_continu_fileName+=decay; pi0_continu_fileName+=decay;
  dsPlusM_physics1_fileName+=decay; DeltaE_physics1_fileName+=decay; MBC_physics1_fileName+=decay; DeltaM_physics1_fileName+=decay; diffD0_physics1_fileName+=decay; dPhi_physics1_fileName+=decay; dPhi_diffD0_physics1_fileName+=decay; pi0_physics1_fileName+=decay;
  dsPlusM_physics3_fileName+=decay; DeltaE_physics3_fileName+=decay; MBC_physics3_fileName+=decay; DeltaM_physics3_fileName+=decay; diffD0_physics3_fileName+=decay; dPhi_physics3_fileName+=decay; dPhi_diffD0_physics3_fileName+=decay; pi0_physics3_fileName+=decay;
  dsPlusM_generic_fileName+=".png"; DeltaE_generic_fileName+=".png"; MBC_generic_fileName+=".png"; DeltaM_generic_fileName+=".png"; diffD0_generic_fileName+=".png"; dPhi_generic_fileName+=".png"; dPhi_diffD0_generic_fileName+=".png"; pi0_generic_fileName+=".png";
  dsPlusM_continu_fileName+=".png"; DeltaE_continu_fileName+=".png"; MBC_continu_fileName+=".png"; DeltaM_continu_fileName+=".png"; diffD0_continu_fileName+=".png"; dPhi_continu_fileName+=".png"; dPhi_diffD0_continu_fileName+=".png"; pi0_continu_fileName+=".png";
  dsPlusM_physics1_fileName+=".png"; DeltaE_physics1_fileName+=".png"; MBC_physics1_fileName+=".png"; DeltaM_physics1_fileName+=".png"; diffD0_physics1_fileName+=".png"; dPhi_physics1_fileName+=".png"; dPhi_diffD0_physics1_fileName+=".png"; pi0_physics1_fileName+=".png";
  dsPlusM_physics3_fileName+=".png"; DeltaE_physics3_fileName+=".png"; MBC_physics3_fileName+=".png"; DeltaM_physics3_fileName+=".png"; diffD0_physics3_fileName+=".png"; dPhi_physics3_fileName+=".png"; dPhi_diffD0_physics3_fileName+=".png"; pi0_physics3_fileName+=".png";
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

bool DeltaECut(float DeltaE)
{
  return (fabs(DeltaE-deltaECut_center)<deltaECut_range);
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

bool pi0MassCut(float pi0Mass)
{
  return (fabs(pi0Mass-pi0MassCut_center)>pi0MassCut_range);
}

bool electronEnergyCut(double electronEnergy, double positronEnergy)
{
  return (electronEnergy<electronEnergyThreshold && positronEnergy<electronEnergyThreshold);
}

double convert(double r1, double r2, double d1, double d2, double phi1, double phi2)
{
  double a1=r1+d1;
  double a2=r2+d2;
  double theta=phi2-phi1;
  double R=a1*a2*sin(theta)/pow(a1*a1+a2*a2+2*a1*a2*cos(theta), 0.5);
  return R;
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

int DsTaggedAnalysis_GenericContinuumPhysics()
{
  setValues();

  TChain *genericTree = new TChain("DsTaggedDecaysProc/nt1");
  // genericTree->Add("/nfs/cor/an2/souvik/Dataset39/DsTaggedProc_GenericMC_213586_214863.root");
  // genericTree->Add("/nfs/cor/an2/souvik/Dataset40/DsTaggedProc_GenericMC_215307_217385.root");
  // genericTree->Add("/nfs/cor/an2/souvik/Dataset41/DsTaggedProc_GenericMC_217687_219721.root");
  genericTree->Add("/nfs/cor/an2/souvik/Dataset47/DsTaggedProc_GenericMC_230474_232255.root");
  genericTree->Add("/nfs/cor/an2/souvik/Dataset48/DsTaggedProc_GenericMC_232264_234607.root");
  
  TChain *continuTree = new TChain("DsTaggedDecaysProc/nt1");
  // continuTree->Add("/nfs/cor/an2/souvik/Dataset39/DsTaggedProc_ContinuumMC_213586_214863.root");
  // continuTree->Add("/nfs/cor/an2/souvik/Dataset40/DsTaggedProc_ContinuumMC_215307_217385.root");
  // continuTree->Add("/nfs/cor/an2/souvik/Dataset41/DsTaggedProc_ContinuumMC_217687_219721.root");
  continuTree->Add("/nfs/cor/an2/souvik/Dataset47/DsTaggedProc_ContinuumMC_230474_232255.root");
  continuTree->Add("/nfs/cor/an2/souvik/Dataset48/DsTaggedProc_ContinuumMC_232264_234607.root");
  
  TChain *physics1Tree = new TChain("DsTaggedDecaysProc/nt1");
  // physics1Tree->Add("/nfs/cor/an2/souvik/Dataset39/DsTaggedProc_Data_213586_214863.root");
  // physics1Tree->Add("/nfs/cor/an2/souvik/Dataset40/DsTaggedProc_Data_215307_217385.root");
  // physics1Tree->Add("/nfs/cor/an2/souvik/Dataset41/DsTaggedProc_Data_217687_219721.root");
  physics1Tree->Add("/nfs/cor/an2/souvik/Dataset47/DsTaggedProc_Data_230474_232255.root");
  physics1Tree->Add("/nfs/cor/an2/souvik/Dataset48/DsTaggedDecaysProc_OriginalData_data48.root");
  
  TChain *physics3Tree = new TChain("DsTaggedDecaysProc/nt3");
  // physics3Tree->Add("/nfs/cor/an2/souvik/Dataset39/DsTaggedProc_Data_213586_214863.root");
  // physics3Tree->Add("/nfs/cor/an2/souvik/Dataset40/DsTaggedProc_Data_215307_217385.root");
  // physics3Tree->Add("/nfs/cor/an2/souvik/Dataset41/DsTaggedProc_Data_217687_219721.root");
  physics3Tree->Add("/nfs/cor/an2/souvik/Dataset47/DsTaggedDecaysProc_ReTaggedData_230474_232255.root");
  physics3Tree->Add("/nfs/cor/an2/souvik/Dataset48/DsTaggedDecaysProc_ReTaggedData_232264_234607.root");
  
  NEvents genericEvents={0,0,0,0,0,0,0,0};
  NEvents continuEvents={0,0,0,0,0,0,0,0};
  NEvents physics1Events={0,0,0,0,0,0,0,0};
  NEvents physics3Events={0,0,0,0,0,0,0,0};
  
  EventNumber genericNumber={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  EventNumber continuNumber={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  EventNumber physics1Number={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  EventNumber physics3Number={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  float runNumber_generic, eventNumber_generic;
  float runNumber_continu, eventNumber_continu;
  float runNumber_physics1, eventNumber_physics1;
  float runNumber_physics3, eventNumber_physics3;
  float dsPlusM_generic, dsPlusCharge_generic, DeltaE_generic, MBC_generic, DeltaM_generic, decayMode_generic;
  float dsPlusM_continu, dsPlusCharge_continu, DeltaE_continu, MBC_continu, DeltaM_continu, decayMode_continu;
  float dsPlusM_physics1, dsPlusCharge_physics1, DeltaE_physics1, MBC_physics1, DeltaM_physics1, decayMode_physics1;
  float dsPlusM_physics3, dsPlusCharge_physics3, DeltaE_physics3, MBC_physics3, DeltaM_physics3, decayMode_physics3;
  float d0_e_generic, d0_p_generic, z0_e_generic, z0_p_generic, px_e_generic, py_e_generic, pz_e_generic, px_p_generic, py_p_generic, pz_p_generic, E_e_generic, E_p_generic, curv_e_generic, curv_p_generic;
  float d0_e_continu, d0_p_continu, z0_e_continu, z0_p_continu, px_e_continu, py_e_continu, pz_e_continu, px_p_continu, py_p_continu, pz_p_continu, E_e_continu, E_p_continu, curv_e_continu, curv_p_continu;
  float d0_e_physics1, d0_p_physics1, z0_e_physics1, z0_p_physics1, px_e_physics1, py_e_physics1, pz_e_physics1, px_p_physics1, py_p_physics1, pz_p_physics1, E_e_physics1, E_p_physics1, curv_e_physics1, curv_p_physics1;
  float d0_e_physics3, d0_p_physics3, z0_e_physics3, z0_p_physics3, px_e_physics3, py_e_physics3, pz_e_physics3, px_p_physics3, py_p_physics3, pz_p_physics3, E_e_physics3, E_p_physics3, curv_e_physics3, curv_p_physics3;
  float px_e_generic_MC, py_e_generic_MC, pz_e_generic_MC, E_e_generic_MC, px_p_generic_MC, py_p_generic_MC, pz_p_generic_MC, E_p_generic_MC;
  float px_e_continu_MC, py_e_continu_MC, pz_e_continu_MC, E_e_continu_MC, px_p_continu_MC, py_p_continu_MC, pz_p_continu_MC, E_p_continu_MC;
  float pi0Mass_generic;
  float pi0Mass_continu;
  float pi0Mass_physics1;
  float pi0Mass_physics3;
  
  typedef std::map<int, float> DecayMap;
  DecayMap decayFrequency_generic, decayFrequency_continu, decayFrequency_physics1, decayFrequency_physics3;
  
  TH1D *h_dsPlusM_generic = new TH1D("h_dsPlusM_generic", "m_{D_{S}^{+}} generic Sample; GeV", 100, 1.9, 2.1); h_dsPlusM_generic->SetLineColor(kRed);
  TH1D *h_dsPlusM_continu = new TH1D("h_dsPlusM_continu", "m_{D_{S}^{+}} continuum Background Sample; GeV", 100, 1.9, 2.1); h_dsPlusM_continu->SetLineColor(kBlue);
  TH1D *h_dsPlusM_physics1 = new TH1D("h_dsPlusM_physics1", "m_{D_{S}^{+}} Data; GeV", 100, 1.9, 2.1); h_dsPlusM_physics1->SetLineColor(kGreen);
  TH1D *h_dsPlusM_physics3 = new TH1D("h_dsPlusM_physics3", "m_{D_{S}^{+}} Data; GeV", 100, 1.9, 2.1); h_dsPlusM_physics3->SetLineColor(kGreen);
  // TH1D *h_DeltaE_generic = new TH1D("h_DeltaE_generic", "#DeltaE generic Sample; GeV", 100, -0.1, 0.2); h_DeltaE_generic->SetLineColor(kRed);
  // TH1D *h_DeltaE_continu = new TH1D("h_DeltaE_continu", "#DeltaE continuum Background Sample; GeV", 100, -0.1, 0.2); h_DeltaE_continu->SetLineColor(kBlue);
  // TH1D *h_DeltaE_physics1 = new TH1D("h_DeltaE_physics1", "#DeltaE Data; GeV", 100, -0.1, 0.2); h_DeltaE_physics1->SetLineColor(kGreen);
  // TH1D *h_DeltaE_physics3 = new TH1D("h_DeltaE_physics3", "#DeltaE Data; GeV", 100, -0.1, 0.2); h_DeltaE_physics3->SetLineColor(kGreen);
  TH1D *h_MBC_generic = new TH1D("h_MBC_generic", "m_{BC} generic Sample; GeV", 100, 2., 2.2); h_MBC_generic->SetLineColor(kRed);
  TH1D *h_MBC_continu = new TH1D("h_MBC_continu", "m_{BC} continuum Background; GeV", 100, 2., 2.2); h_MBC_continu->SetLineColor(kBlue);
  TH1D *h_MBC_physics1 = new TH1D("h_MBC_physics1", "m_{BC} Data; GeV", 100, 2., 2.2); h_MBC_physics1->SetLineColor(kGreen);
  TH1D *h_MBC_physics3 = new TH1D("h_MBC_physics3", "m_{BC} Data; GeV", 100, 2., 2.2); h_MBC_physics3->SetLineColor(kGreen);
  TH1D *h_DeltaM_generic = new TH1D("h_DeltaM_generic", "#deltaM generic Sample", 100, 0.0, 0.5); h_DeltaM_generic->SetLineColor(kRed);
  TH1D *h_DeltaM_continu = new TH1D("h_DeltaM_continu", "#deltaM continuum Background Sample; GeV", 100, 0.0, 0.5); h_DeltaM_continu->SetLineColor(kBlue);
  TH1D *h_DeltaM_physics1 = new TH1D("h_DeltaM_physics1", "#deltaM Data; GeV", 100, 0.0, 0.5); h_DeltaM_physics1->SetLineColor(kGreen);
  TH1D *h_DeltaM_physics3 = new TH1D("h_DeltaM_physics3", "#deltaM Data; GeV", 100, 0.0, 0.5); h_DeltaM_physics3->SetLineColor(kGreen);
  
  // TH2D *h_DeltaE_MBC_generic = new TH2D("h_DeltaE_MBC_generic", "h_DeltaE_MBC_generic", 100, -0.1, 0.2, 100, 2., 2.2);
  // TH2D *h_DeltaE_MBC_continu = new TH2D("h_DeltaE_MBC_continu", "h_DeltaE_MBC_continu", 100, -0.1, 0.2, 100, 2., 2.2);
  // TH2D *h_DeltaE_MBC_physics1 = new TH2D("h_DeltaE_MBC_physics1", "h_DeltaE_MBC_physics1", 100, -0.1, 0.2, 100, 2., 2.2);
  // TH2D *h_DeltaE_DeltaM_generic = new TH2D("h_DeltaE_DeltaM_generic", "h_DeltaE_DeltaM_generic", 100, -0.1, 0.2, 100, 0.0, 0.2);
  // TH2D *h_DeltaE_DeltaM_continu = new TH2D("h_DeltaE_DeltaM_continu", "h_DeltaE_DeltaM_continu", 100, -0.1, 0.2, 100, 0.0, 0.2);
  // TH2D *h_DeltaE_DeltaM_physics1 = new TH2D("h_DeltaE_DeltaM_physics1", "h_DeltaE_DeltaM_physics1", 100, -0.1, 0.2, 100, 0.0, 0.2);
  // TH2D *h_MBC_DeltaM_generic = new TH2D("h_MBC_DeltaM_generic", "h_MBC_DeltaM_generic", 100, 2., 2.2, 100, 0.0, 0.2);
  // TH2D *h_MBC_DeltaM_continu = new TH2D("h_MBC_DeltaM_continu", "h_MBC_DeltaM_continu", 100, 2., 2.2, 100, 0.0, 0.2);
  // TH2D *h_MBC_DeltaM_physics1 = new TH2D("h_MBC_DeltaM_physics1", "h_MBC_DeltaM_physics1", 100, 2., 2.2, 100, 0.0, 0.2);
  
  // TH1D *h_d0_e_generic = new TH1D("h_d0_e_generic", "d0_e", 50, -0.01, 0.01); h_d0_e_generic->SetLineColor(kRed);
  // TH1D *h_d0_e_continu = new TH1D("h_d0_e_continu", "d0_e", 50, -0.01, 0.01); h_d0_e_continu->SetLineColor(kBlue);
  // TH1D *h_d0_e_physics1 = new TH1D("h_d0_e_physics1", "d0_e", 50, -0.01, 0.01); h_d0_e_physics1->SetLineColor(kGreen);
  // TH1D *h_d0_p_generic = new TH1D("h_d0_p_generic", "d0_p", 50, -0.01, 0.01); h_d0_p_generic->SetLineColor(kRed);
  // TH1D *h_d0_p_continu = new TH1D("h_d0_p_continu", "d0_p", 50, -0.01, 0.01); h_d0_p_continu->SetLineColor(kBlue);
  // TH1D *h_d0_p_physics1 = new TH1D("h_d0_p_physics1", "d0_p", 50, -0.01, 0.01); h_d0_p_physics1->SetLineColor(kGreen);
  TH1D *h_diffD0_generic = new TH1D("h_diffD0_generic", "#Deltad_{0} generic Sample; m", 50, -0.01, 0.01); h_diffD0_generic->SetLineColor(kRed);
  TH1D *h_diffD0_continu = new TH1D("h_diffD0_continu", "#Deltad_{0} continuum Background Sample; m", 50, -0.01, 0.01); h_diffD0_continu->SetLineColor(kBlue);
  TH1D *h_diffD0_physics1 = new TH1D("h_diffD0_physics1", "#Deltad_{0} Data; m", 50, -0.01, 0.01); h_diffD0_physics1->SetLineColor(kGreen);
  TH1D *h_diffD0_physics3 = new TH1D("h_diffD0_physics3", "#Deltad_{0} Data; m", 50, -0.01, 0.01); h_diffD0_physics3->SetLineColor(kGreen);
  TH1D *h_dPhi_generic = new TH1D("h_dPhi_generic", "#Delta#Phi generic Sample", 50, -2., 2.); h_dPhi_generic->SetLineColor(kRed);
  TH1D *h_dPhi_continu = new TH1D("h_dPhi_continu", "#Delta#Phi continuum Background Sample", 50, -2., 2.); h_dPhi_continu->SetLineColor(kBlue);
  TH1D *h_dPhi_physics1 = new TH1D("h_dPhi_physics1", "#Delta#Phi Data", 50, -2., 2.); h_dPhi_physics1->SetLineColor(kGreen);
  TH1D *h_dPhi_physics3 = new TH1D("h_dPhi_physics3", "#Delta#Phi Data", 50, -2., 2.); h_dPhi_physics3->SetLineColor(kGreen);
  // TH2D *h_d0_phi_e_generic = new TH2D("h_d0_phi_e_generic", "h_d0_phi_e_generic", 50, -pi, pi, 50, -0.01, 0.01);
  // TH2D *h_d0_phi_p_generic = new TH2D("h_d0_phi_p_generic", "h_d0_phi_p_generic", 50, -pi, pi, 50, -0.01, 0.01);
  // TH2D *h_d0_phi_e_continu = new TH2D("h_d0_phi_e_continu", "h_d0_phi_e_continu", 50, -pi, pi, 50, -0.01, 0.01);
  // TH2D *h_d0_phi_p_continu = new TH2D("h_d0_phi_p_continu", "h_d0_phi_p_continu", 50, -pi, pi, 50, -0.01, 0.01);
  // TH2D *h_d0_phi_e_physics1 = new TH2D("h_d0_phi_e_physics1", "h_d0_phi_e_physics1", 50, -pi, pi, 50, -0.01, 0.01);
  // TH2D *h_d0_phi_p_physics1 = new TH2D("h_d0_phi_p_physics1", "h_d0_phi_p_physics1", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_dPhi_diffD0_generic = new TH2D("h_dPhi_diffD0_generic", "#Delta#Phi vs #Deltad_{0} generic Sample; ; m", 50, -2., 2., 50, -0.01, 0.01); h_dPhi_diffD0_generic->SetLineColor(kRed);
  TH2D *h_dPhi_diffD0_continu = new TH2D("h_dPhi_diffD0_continu", "#Delta#Phi vs #Deltad_{0} continuum Background Sample; ; m", 50, -2., 2., 50, -0.01, 0.01); h_dPhi_diffD0_continu->SetLineColor(kBlue);
  TH2D *h_dPhi_diffD0_physics1 = new TH2D("h_dPhi_diffD0_physics1", "#Delta#Phi vs #Deltad_{0} Data; ; m", 50, -2., 2., 50, -0.01, 0.01); h_dPhi_diffD0_physics1->SetLineColor(kGreen);
  TH2D *h_dPhi_diffD0_physics3 = new TH2D("h_dPhi_diffD0_physics3", "#Delta#Phi vs #Deltad_{0} Data; ; m", 50, -2., 2., 50, -0.01, 0.01); h_dPhi_diffD0_physics3->SetLineColor(kGreen);
  //TH2D *h_electronE = new TH2D("h_electronE", "Electron Energy Resolution; Reco Electron Energy (GeV); reco-MC", 100, 0., 0.2, 100, -0.1, 0.1);
  //TH2D *h_electronPx = new TH2D("h_electronPx", "Electron Px Resolution; Reco Electron Px (GeV)", 100, 0.0, 0.14, 100, -0.03, 0.03); 
  //TH2D *h_electronPy = new TH2D("h_electronPy", "Electron Py Resolution; Reco Electron Py (GeV)", 100, 0.0, 0.14, 100, -0.03, 0.03); 
  //TH2D *h_electronPz = new TH2D("h_electronPz", "Electron Pz Resolution; Reco Electron Pz (GeV)", 100, 0.0, 0.14, 100, -0.03, 0.03);
  TH1D *h_electronE;
  if (DeltaM_sideband) h_electronE = new TH1D("h_electronE", "Electron Energy; Reco Electron Energy (GeV)", 100, 0.0, 0.5); 
  else h_electronE = new TH1D("h_electronE", "Electron Energy; Reco Electron Energy (GeV)", 100, 0.0, 0.2); 
  h_electronE->SetLineColor(kRed);
  
  TH1D *h_pi0_generic = new TH1D("h_pi0_generic", "Pion Mass", 50, 0.035, .235); h_pi0_generic->SetLineColor(kRed);
  TH1D *h_pi0_continu = new TH1D("h_pi0_continu", "Pion Mass", 50, 0.035, .235); h_pi0_continu->SetLineColor(kBlue);
  TH1D *h_pi0_physics1 = new TH1D("h_pi0_physics1", "Pion Mass", 50, 0.035, .235); h_pi0_physics1->SetLineColor(kGreen);
  TH1D *h_pi0_physics3 = new TH1D("h_pi0_physics3", "Pion Mass", 50, 0.035, .235); h_pi0_physics3->SetLineColor(kGreen);
  
  genericTree->SetBranchAddress("Run", &(runNumber_generic));
  genericTree->SetBranchAddress("Event", &(eventNumber_generic));
  genericTree->SetBranchAddress("dsPlusM", &(dsPlusM_generic));
  genericTree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_generic));
  genericTree->SetBranchAddress("DecayMode", &(decayMode_generic));
  genericTree->SetBranchAddress("MBC", &(MBC_generic));
  genericTree->SetBranchAddress("DeltaM", &(DeltaM_generic));
  genericTree->SetBranchAddress("kElectron1E_MC", &(E_e_generic_MC));  
  genericTree->SetBranchAddress("kElectron1Px_MC", &(px_e_generic_MC));
  genericTree->SetBranchAddress("kElectron1Py_MC", &(py_e_generic_MC));
  genericTree->SetBranchAddress("kElectron1Pz_MC", &(pz_e_generic_MC));
  genericTree->SetBranchAddress("kElectron2E_MC", &(E_p_generic_MC));
  genericTree->SetBranchAddress("kElectron2Px_MC", &(px_p_generic_MC));
  genericTree->SetBranchAddress("kElectron2Py_MC", &(py_p_generic_MC));
  genericTree->SetBranchAddress("kElectron2Pz_MC", &(pz_p_generic_MC));
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
  continuTree->SetBranchAddress("MBC", &(MBC_continu));
  continuTree->SetBranchAddress("DeltaM", &(DeltaM_continu));
  continuTree->SetBranchAddress("kElectron1E_MC", &(E_e_continu_MC));  
  continuTree->SetBranchAddress("kElectron1Px_MC", &(px_e_continu_MC));
  continuTree->SetBranchAddress("kElectron1Py_MC", &(py_e_continu_MC));
  continuTree->SetBranchAddress("kElectron1Pz_MC", &(pz_e_continu_MC));
  continuTree->SetBranchAddress("kElectron2E_MC", &(E_p_continu_MC));
  continuTree->SetBranchAddress("kElectron2Px_MC", &(px_p_continu_MC));
  continuTree->SetBranchAddress("kElectron2Py_MC", &(py_p_continu_MC));
  continuTree->SetBranchAddress("kElectron2Pz_MC", &(pz_p_continu_MC));
  continuTree->SetBranchAddress("kElectron1D0_reco", &(d0_e_continu));
  continuTree->SetBranchAddress("kElectron2D0_reco", &(d0_p_continu));
  continuTree->SetBranchAddress("kElectron1Z0_reco", &(z0_e_continu));
  continuTree->SetBranchAddress("kElectron2Z0_reco", &(z0_p_continu));
  continuTree->SetBranchAddress("kElectron1Px_reco", &(px_e_continu));
  continuTree->SetBranchAddress("kElectron1Py_reco", &(py_e_continu));
  continuTree->SetBranchAddress("kElectron1Pz_reco", &(pz_e_continu));
  continuTree->SetBranchAddress("kElectron2Px_reco", &(px_p_continu));
  continuTree->SetBranchAddress("kElectron2Py_reco", &(py_p_continu));
  continuTree->SetBranchAddress("kElectron2Pz_reco", &(pz_p_continu));
  continuTree->SetBranchAddress("kElectron1E_reco", &(E_e_continu));
  continuTree->SetBranchAddress("kElectron2E_reco", &(E_p_continu));
  continuTree->SetBranchAddress("kElectron1Curv_reco", &(curv_e_continu));
  continuTree->SetBranchAddress("kElectron2Curv_reco", &(curv_p_continu));
  continuTree->SetBranchAddress("kPi0Mass_reco", &(pi0Mass_continu));
  
  physics1Tree->SetBranchAddress("Run", &(runNumber_physics1));
  physics1Tree->SetBranchAddress("Event", &(eventNumber_physics1));
  physics1Tree->SetBranchAddress("dsPlusM", &(dsPlusM_physics1));
  physics1Tree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_physics1));
  physics1Tree->SetBranchAddress("DecayMode", &(decayMode_physics1));
  physics1Tree->SetBranchAddress("MBC", &(MBC_physics1));
  physics1Tree->SetBranchAddress("DeltaM", &(DeltaM_physics1));
  physics1Tree->SetBranchAddress("kElectron1D0_reco", &(d0_e_physics1));
  physics1Tree->SetBranchAddress("kElectron2D0_reco", &(d0_p_physics1));
  physics1Tree->SetBranchAddress("kElectron1Z0_reco", &(z0_e_physics1));
  physics1Tree->SetBranchAddress("kElectron2Z0_reco", &(z0_p_physics1));
  physics1Tree->SetBranchAddress("kElectron1Px_reco", &(px_e_physics1));
  physics1Tree->SetBranchAddress("kElectron1Py_reco", &(py_e_physics1));
  physics1Tree->SetBranchAddress("kElectron1Pz_reco", &(pz_e_physics1));
  physics1Tree->SetBranchAddress("kElectron2Px_reco", &(px_p_physics1));
  physics1Tree->SetBranchAddress("kElectron2Py_reco", &(py_p_physics1));
  physics1Tree->SetBranchAddress("kElectron2Pz_reco", &(pz_p_physics1));
  physics1Tree->SetBranchAddress("kElectron1E_reco", &(E_e_physics1));
  physics1Tree->SetBranchAddress("kElectron2E_reco", &(E_p_physics1));
  physics1Tree->SetBranchAddress("kElectron1Curv_reco", &(curv_e_physics1));
  physics1Tree->SetBranchAddress("kElectron2Curv_reco", &(curv_p_physics1));
  physics1Tree->SetBranchAddress("kPi0Mass_reco", &(pi0Mass_physics1));
  
  physics3Tree->SetBranchAddress("Run", &(runNumber_physics3));
  physics3Tree->SetBranchAddress("Event", &(eventNumber_physics3));
  physics3Tree->SetBranchAddress("dsPlusM", &(dsPlusM_physics3));
  physics3Tree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_physics3));
  physics3Tree->SetBranchAddress("DecayMode", &(decayMode_physics3));
  physics3Tree->SetBranchAddress("MBC", &(MBC_physics3));
  physics3Tree->SetBranchAddress("DeltaM", &(DeltaM_physics3));
  physics3Tree->SetBranchAddress("kElectron1D0_reco", &(d0_e_physics3));
  physics3Tree->SetBranchAddress("kElectron2D0_reco", &(d0_p_physics3));
  physics3Tree->SetBranchAddress("kElectron1Z0_reco", &(z0_e_physics3));
  physics3Tree->SetBranchAddress("kElectron2Z0_reco", &(z0_p_physics3));
  physics3Tree->SetBranchAddress("kElectron1Px_reco", &(px_e_physics3));
  physics3Tree->SetBranchAddress("kElectron1Py_reco", &(py_e_physics3));
  physics3Tree->SetBranchAddress("kElectron1Pz_reco", &(pz_e_physics3));
  physics3Tree->SetBranchAddress("kElectron2Px_reco", &(px_p_physics3));
  physics3Tree->SetBranchAddress("kElectron2Py_reco", &(py_p_physics3));
  physics3Tree->SetBranchAddress("kElectron2Pz_reco", &(pz_p_physics3));
  physics3Tree->SetBranchAddress("kElectron1E_reco", &(E_e_physics3));
  physics3Tree->SetBranchAddress("kElectron2E_reco", &(E_p_physics3));
  physics3Tree->SetBranchAddress("kElectron1Curv_reco", &(curv_e_physics3));
  physics3Tree->SetBranchAddress("kElectron2Curv_reco", &(curv_p_physics3));
  physics3Tree->SetBranchAddress("kPi0Mass_reco", &(pi0Mass_physics3));
  
  int ngenericEvents=genericTree->GetEntries();
  int oppRecon_generic=0;
  for (int i=0; i<ngenericEvents; ++i)
  {
    genericTree->GetEvent(i);
    
    double phi_e=atan2(py_e_generic, px_e_generic);
    double phi_p=atan2(py_p_generic, px_p_generic);
    double dPhi=deltaPhi(phi_e, phi_p);
    
    if (runNumber_generic!=genericNumber.noCut_run || eventNumber_generic!=genericNumber.noCut_event)
    {      
      if (decayMode_generic>-1) decayFrequency_generic[int(decayMode_generic)]+=1;
      ++genericEvents.noCut;
      genericNumber.noCut_run=runNumber_generic;
      genericNumber.noCut_event=eventNumber_generic;
    }
    
    if (fabs(d0_e_generic)<0.005 && fabs(d0_p_generic)<0.005 &&
        fabs(z0_e_generic)<0.05 && fabs(z0_p_generic)<0.05)
    {
    
    if (decayMode_generic==decayNumber)
    {
      if (runNumber_generic!=genericNumber.tagCut_run || eventNumber_generic!=genericNumber.tagCut_event)
      {      
        ++genericEvents.tagCut;
        genericNumber.tagCut_run=runNumber_generic;
        genericNumber.tagCut_event=eventNumber_generic;
      }
      if (!DsMass_sideband) h_dsPlusM_generic->Fill(dsPlusM_generic);
      if (DsMass_sideband || dsPlusMCut(dsPlusM_generic))
      {
        if (runNumber_generic!=genericNumber.dsPlusMCut_run || eventNumber_generic!=genericNumber.dsPlusMCut_event)
        {
          ++genericEvents.dsPlusMCut;
          genericNumber.dsPlusMCut_run=runNumber_generic;
          genericNumber.dsPlusMCut_event=eventNumber_generic;
        }
        if (!mBC_sideband) h_MBC_generic->Fill(MBC_generic);
        if (mBC_sideband || MBCCut(MBC_generic))
        {
          if (runNumber_generic!=genericNumber.mbcCut_run || eventNumber_generic!=genericNumber.mbcCut_event)
          {
            ++genericEvents.mbcCut;
            genericNumber.mbcCut_run=runNumber_generic;
            genericNumber.mbcCut_event=eventNumber_generic;
          }
          if (!DeltaM_sideband) h_DeltaM_generic->Fill(DeltaM_generic);
          if (DeltaM_sideband || DeltaMCut(DeltaM_generic))
          {
            if (runNumber_generic!=genericNumber.deltaMCut_run || eventNumber_generic!=genericNumber.deltaMCut_event)
            {
              ++genericEvents.deltaMCut;
              genericNumber.deltaMCut_run=runNumber_generic;
              genericNumber.deltaMCut_event=eventNumber_generic;
            }
            h_diffD0_generic->Fill(d0_e_generic-d0_p_generic);
            h_dPhi_diffD0_generic->Fill(dPhi, d0_e_generic-d0_p_generic);
            if (dD0(d0_e_generic-d0_p_generic))
            {
              if (runNumber_generic!=genericNumber.diffD0Cut_run || eventNumber_generic!=genericNumber.diffD0Cut_event)
              {
                ++genericEvents.diffD0Cut;
                genericNumber.diffD0Cut_run=runNumber_generic;
                genericNumber.diffD0Cut_event=eventNumber_generic;
              }
              h_dPhi_generic->Fill(dPhi);
              if (dPhiCut(dPhi))
              {
                if (runNumber_generic!=genericNumber.dPhiCut_run || eventNumber_generic!=genericNumber.dPhiCut_event)
                {
                  ++genericEvents.dPhiCut;
                  genericNumber.dPhiCut_run=runNumber_generic;
                  genericNumber.dPhiCut_event=eventNumber_generic;
                  //std::cout<<runNumber_generic<<" "<<eventNumber_generic<<std::endl;
                }
                h_pi0_generic->Fill(pi0Mass_generic);
                if (pi0MassCut(pi0Mass_generic))
                {
                  if (runNumber_generic!=genericNumber.pi0MassCut_run || eventNumber_generic!=genericNumber.pi0MassCut_event)
                  {
                    ++genericEvents.pi0MassCut;
                    genericNumber.pi0MassCut_run=runNumber_generic;
                    genericNumber.pi0MassCut_event=eventNumber_generic;
                  }
                  if (dsPlusCharge_generic<0) oppRecon_generic+=1;
                  if (DeltaM_sideband)
                  {
                    if (electronEnergyCut(E_e_generic, E_p_generic))
                    {
                      h_DeltaM_generic->Fill(DeltaM_generic);
                      if (!DeltaMCut(DeltaM_generic))
                      {
                        h_electronE->Fill(E_e_generic);
                        h_electronE->Fill(E_p_generic);
                      }
                    }
                  }
                  else if (mBC_sideband)
                  {
                    h_MBC_generic->Fill(MBC_generic);
                    if (!MBCCut(MBC_generic))
                    {
                      h_electronE->Fill(E_e_generic);
                      h_electronE->Fill(E_p_generic);
                    }
                  }
                  else if (DsMass_sideband)
                  {
                    h_dsPlusM_generic->Fill(dsPlusM_generic);
                    if (!dsPlusMCut(dsPlusM_generic))
                    {
                      h_electronE->Fill(E_e_generic);
                      h_electronE->Fill(E_p_generic);
                    }
                  }
                  else
                  {
                    h_electronE->Fill(E_e_generic);
                    h_electronE->Fill(E_p_generic);
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
  
  
  int ncontinuEvents=continuTree->GetEntries();
  int oppRecon_continu=0;
  for (int i=0; i<ncontinuEvents; ++i)
  {
    continuTree->GetEvent(i);
    
    double phi_e=atan2(py_e_continu, px_e_continu);
    double phi_p=atan2(py_p_continu, px_p_continu);
    double dPhi=deltaPhi(phi_e, phi_p);
    
    if (runNumber_continu!=continuNumber.noCut_run || eventNumber_continu!=continuNumber.noCut_event)
    {      
      if (decayMode_continu>-1) decayFrequency_continu[int(decayMode_continu)]+=1;
      ++continuEvents.noCut;
      continuNumber.noCut_run=runNumber_continu;
      continuNumber.noCut_event=eventNumber_continu;
    }
    
    if (fabs(d0_e_continu)<0.005 && fabs(d0_p_continu)<0.005 &&
        fabs(z0_e_continu)<0.05 && fabs(z0_p_continu)<0.05)
    {
    
    if (decayMode_continu==decayNumber)
    {
      if (runNumber_continu!=continuNumber.tagCut_run || eventNumber_continu!=continuNumber.tagCut_event)
      {      
        ++continuEvents.tagCut;
        continuNumber.tagCut_run=runNumber_continu;
        continuNumber.tagCut_event=eventNumber_continu;
      }
      if (!DsMass_sideband) h_dsPlusM_continu->Fill(dsPlusM_continu);
      if (DsMass_sideband || dsPlusMCut(dsPlusM_continu))
      {
        if (runNumber_continu!=continuNumber.dsPlusMCut_run || eventNumber_continu!=continuNumber.dsPlusMCut_event)
        {
          ++continuEvents.dsPlusMCut;
          continuNumber.dsPlusMCut_run=runNumber_continu;
          continuNumber.dsPlusMCut_event=eventNumber_continu;
        }
        if (!mBC_sideband) h_MBC_continu->Fill(MBC_continu);
        if (mBC_sideband || MBCCut(MBC_continu))
        {
          if (runNumber_continu!=continuNumber.mbcCut_run || eventNumber_continu!=continuNumber.mbcCut_event)
          {
            ++continuEvents.mbcCut;
            continuNumber.mbcCut_run=runNumber_continu;
            continuNumber.mbcCut_event=eventNumber_continu;
          }
          if (!DeltaM_sideband) h_DeltaM_continu->Fill(DeltaM_continu);
          if (DeltaM_sideband || DeltaMCut(DeltaM_continu))
          {
            if (runNumber_continu!=continuNumber.deltaMCut_run || eventNumber_continu!=continuNumber.deltaMCut_event)
            {
              ++continuEvents.deltaMCut;
              continuNumber.deltaMCut_run=runNumber_continu;
              continuNumber.deltaMCut_event=eventNumber_continu;
            }
            h_diffD0_continu->Fill(d0_e_continu-d0_p_continu);
            h_dPhi_diffD0_continu->Fill(dPhi, d0_e_continu-d0_p_continu);
            if (dD0(d0_e_continu-d0_p_continu))
            {
              if (runNumber_continu!=continuNumber.diffD0Cut_run || eventNumber_continu!=continuNumber.diffD0Cut_event)
              {
                ++continuEvents.diffD0Cut;
                continuNumber.diffD0Cut_run=runNumber_continu;
                continuNumber.diffD0Cut_event=eventNumber_continu;
              }
              h_dPhi_continu->Fill(dPhi);
              if (dPhiCut(dPhi))
              {
                if (runNumber_continu!=continuNumber.dPhiCut_run || eventNumber_continu!=continuNumber.dPhiCut_event)
                {
                  ++continuEvents.dPhiCut;
                  continuNumber.dPhiCut_run=runNumber_continu;
                  continuNumber.dPhiCut_event=eventNumber_continu;
                }
                h_pi0_continu->Fill(pi0Mass_continu);
                if (pi0MassCut(pi0Mass_continu))
                {
                  if (runNumber_continu!=continuNumber.pi0MassCut_run || eventNumber_continu!=continuNumber.pi0MassCut_event)
                  {
                    ++continuEvents.pi0MassCut;
                    continuNumber.pi0MassCut_run=runNumber_continu;
                    continuNumber.pi0MassCut_event=eventNumber_continu;
                  }
                  if (dsPlusCharge_continu<0) oppRecon_continu+=1;
                  if (DeltaM_sideband)
                  {
                    if (electronEnergyCut(E_e_continu, E_p_continu))
                    {
                      h_DeltaM_continu->Fill(DeltaM_continu);
                    }
                  }
                  if (mBC_sideband) h_MBC_continu->Fill(MBC_continu);
                  if (DsMass_sideband) h_dsPlusM_continu->Fill(dsPlusM_continu);
                }
              }
            }
          }
        }
      }
    }
    }
  }
  
  
  int nphysics1Events=physics1Tree->GetEntries();
  int oppRecon_physics1=0;
  for (int i=0; i<nphysics1Events; ++i)
  {
    physics1Tree->GetEvent(i);
    
    double phi_e=atan2(py_e_physics1, px_e_physics1);
    double phi_p=atan2(py_p_physics1, px_p_physics1);
    double dPhi=deltaPhi(phi_e, phi_p);
    
    if (runNumber_physics1!=physics1Number.noCut_run || eventNumber_physics1!=physics1Number.noCut_event)
    {      
      if (decayMode_physics1>-1) decayFrequency_physics1[int(decayMode_physics1)]+=1;
      ++physics1Events.noCut;
      physics1Number.noCut_run=runNumber_physics1;
      physics1Number.noCut_event=eventNumber_physics1;
    }
    
    if (fabs(d0_e_physics1)<0.005 && fabs(d0_p_physics1)<0.005 && 
        fabs(z0_e_physics1)<0.05 && fabs(z0_p_physics1)<0.05)
    {
    
    if (decayMode_physics1==decayNumber)
    {
      if (runNumber_physics1!=physics1Number.tagCut_run || eventNumber_physics1!=physics1Number.tagCut_event)
      {      
        ++physics1Events.tagCut;
        physics1Number.tagCut_run=runNumber_physics1;
        physics1Number.tagCut_event=eventNumber_physics1;
      }
      if (!DsMass_sideband) h_dsPlusM_physics1->Fill(dsPlusM_physics1);
      if (DsMass_sideband || dsPlusMCut(dsPlusM_physics1))
      {
        if (runNumber_physics1!=physics1Number.dsPlusMCut_run || eventNumber_physics1!=physics1Number.dsPlusMCut_event)
        {
          ++physics1Events.dsPlusMCut;
          physics1Number.dsPlusMCut_run=runNumber_physics1;
          physics1Number.dsPlusMCut_event=eventNumber_physics1;
        }
        if (runNumber_physics1!=physics1Number.deltaECut_run || eventNumber_physics1!=physics1Number.deltaECut_event)
        {
          ++physics1Events.deltaECut;
          physics1Number.deltaECut_run=runNumber_physics1;
          physics1Number.deltaECut_event=eventNumber_physics1;
        }
        if (!mBC_sideband) h_MBC_physics1->Fill(MBC_physics1);
        if (mBC_sideband || MBCCut(MBC_physics1))
        {
          if (runNumber_physics1!=physics1Number.mbcCut_run || eventNumber_physics1!=physics1Number.mbcCut_event)
          {
            ++physics1Events.mbcCut;
            physics1Number.mbcCut_run=runNumber_physics1;
            physics1Number.mbcCut_event=eventNumber_physics1;
          }
          if (!DeltaM_sideband) h_DeltaM_physics1->Fill(DeltaM_physics1);
          if (DeltaM_sideband || DeltaMCut(DeltaM_physics1))
          {
            if (runNumber_physics1!=physics1Number.deltaMCut_run || eventNumber_physics1!=physics1Number.deltaMCut_event)
            {
              ++physics1Events.deltaMCut;
              physics1Number.deltaMCut_run=runNumber_physics1;
              physics1Number.deltaMCut_event=eventNumber_physics1;
            }
            h_diffD0_physics1->Fill(d0_e_physics1-d0_p_physics1);
            h_dPhi_diffD0_physics1->Fill(dPhi, d0_e_physics1-d0_p_physics1);
            if (dD0(d0_e_physics1-d0_p_physics1))
            {
              if (runNumber_physics1!=physics1Number.diffD0Cut_run || eventNumber_physics1!=physics1Number.diffD0Cut_event)
              {
                ++physics1Events.diffD0Cut;
                physics1Number.diffD0Cut_run=runNumber_physics1;
                physics1Number.diffD0Cut_event=eventNumber_physics1;
              }
              h_dPhi_physics1->Fill(dPhi);
              if (dPhiCut(dPhi))
              {
                if (runNumber_physics1!=physics1Number.dPhiCut_run || eventNumber_physics1!=physics1Number.dPhiCut_event)
                {
                  ++physics1Events.dPhiCut;
                  physics1Number.dPhiCut_run=runNumber_physics1;
                  physics1Number.dPhiCut_event=eventNumber_physics1;
                }
                h_pi0_physics1->Fill(pi0Mass_physics1);
                if (pi0MassCut(pi0Mass_physics1))
                {
                  if (runNumber_physics1!=physics1Number.pi0MassCut_run || eventNumber_physics1!=physics1Number.pi0MassCut_event)
                  {
                    ++physics1Events.pi0MassCut;
                    physics1Number.pi0MassCut_run=runNumber_physics1;
                    physics1Number.pi0MassCut_event=eventNumber_physics1;
                  }
                  if (dsPlusCharge_physics1<0) oppRecon_physics1+=1;
                  if (DeltaM_sideband)
                  {
                    if (electronEnergyCut(E_e_physics1, E_p_physics1))
                    {
                      h_DeltaM_physics1->Fill(DeltaM_physics1);
                    }
                  }
                  if (mBC_sideband) h_MBC_physics1->Fill(MBC_physics1);
                  if (DsMass_sideband) h_dsPlusM_physics1->Fill(dsPlusM_physics1);
                }
              }
            }
          }
        }
      }
    }
    }
  }
  
  int nphysics3Events=physics3Tree->GetEntries();
  int oppRecon_physics3=0;
  for (int i=0; i<nphysics3Events; ++i)
  {
    physics3Tree->GetEvent(i);
    
    double phi_e=atan2(py_e_physics3, px_e_physics3);
    double phi_p=atan2(py_p_physics3, px_p_physics3);
    double dPhi=deltaPhi(phi_e, phi_p);
    
    if (runNumber_physics3!=physics3Number.noCut_run || eventNumber_physics3!=physics3Number.noCut_event)
    {      
      if (decayMode_physics3>-1) decayFrequency_physics3[int(decayMode_physics3)]+=1;
      ++physics3Events.noCut;
      physics3Number.noCut_run=runNumber_physics3;
      physics3Number.noCut_event=eventNumber_physics3;
    }
    
    if (fabs(d0_e_physics3)<0.005 && fabs(d0_p_physics3)<0.005 && 
        fabs(z0_e_physics3)<0.05 && fabs(z0_p_physics3)<0.05)
    {
    
    if (decayMode_physics3==decayNumber)
    {
      if (runNumber_physics3!=physics3Number.tagCut_run || eventNumber_physics3!=physics3Number.tagCut_event)
      {
        ++physics3Events.tagCut;
        physics3Number.tagCut_run=runNumber_physics3;
        physics3Number.tagCut_event=eventNumber_physics3;
      }
      if (!DsMass_sideband) h_dsPlusM_physics3->Fill(dsPlusM_physics3);
      if (DsMass_sideband || dsPlusMCut(dsPlusM_physics3))
      {
        if (runNumber_physics3!=physics3Number.dsPlusMCut_run || eventNumber_physics3!=physics3Number.dsPlusMCut_event)
        {
          ++physics3Events.dsPlusMCut;
          physics3Number.dsPlusMCut_run=runNumber_physics3;
          physics3Number.dsPlusMCut_event=eventNumber_physics3;
        }
        if (runNumber_physics3!=physics3Number.deltaECut_run || eventNumber_physics3!=physics3Number.deltaECut_event)
        {
          ++physics3Events.deltaECut;
          physics3Number.deltaECut_run=runNumber_physics3;
          physics3Number.deltaECut_event=eventNumber_physics3;
        }
        if (!mBC_sideband) h_MBC_physics3->Fill(MBC_physics3);
        if (mBC_sideband || MBCCut(MBC_physics3))
        {
          if (runNumber_physics3!=physics3Number.mbcCut_run || eventNumber_physics3!=physics3Number.mbcCut_event)
          {
            ++physics3Events.mbcCut;
            physics3Number.mbcCut_run=runNumber_physics3;
            physics3Number.mbcCut_event=eventNumber_physics3;
          }
          if (!DeltaM_sideband) h_DeltaM_physics3->Fill(DeltaM_physics3);
          if (DeltaM_sideband || DeltaMCut(DeltaM_physics3))
          {
            if (runNumber_physics3!=physics3Number.deltaMCut_run || eventNumber_physics3!=physics3Number.deltaMCut_event)
            {
              ++physics3Events.deltaMCut;
              physics3Number.deltaMCut_run=runNumber_physics3;
              physics3Number.deltaMCut_event=eventNumber_physics3;
            }
            h_diffD0_physics3->Fill(d0_e_physics3-d0_p_physics3);
            h_dPhi_diffD0_physics3->Fill(dPhi, d0_e_physics3-d0_p_physics3);
            if (dD0(d0_e_physics3-d0_p_physics3))
            {
              if (runNumber_physics3!=physics3Number.diffD0Cut_run || eventNumber_physics3!=physics3Number.diffD0Cut_event)
              {
                ++physics3Events.diffD0Cut;
                physics3Number.diffD0Cut_run=runNumber_physics3;
                physics3Number.diffD0Cut_event=eventNumber_physics3;
              }
              h_dPhi_physics3->Fill(dPhi);
              if (dPhiCut(dPhi))
              {
                if (runNumber_physics3!=physics3Number.dPhiCut_run || eventNumber_physics3!=physics3Number.dPhiCut_event)
                {
                  ++physics3Events.dPhiCut;
                  physics3Number.dPhiCut_run=runNumber_physics3;
                  physics3Number.dPhiCut_event=eventNumber_physics3;
                }
                h_pi0_physics3->Fill(pi0Mass_physics3);
                if (pi0MassCut(pi0Mass_physics3))
                {
                  if (runNumber_physics3!=physics3Number.pi0MassCut_run || eventNumber_physics3!=physics3Number.pi0MassCut_event)
                  {
                    ++physics3Events.pi0MassCut;
                    physics3Number.pi0MassCut_run=runNumber_physics3;
                    physics3Number.pi0MassCut_event=eventNumber_physics3;
                  }
                  if (dsPlusCharge_physics3<0) oppRecon_physics3+=1;
                  if (DeltaM_sideband)
                  {
                    if (electronEnergyCut(E_e_physics3, E_p_physics3))
                    {
                      h_DeltaM_physics3->Fill(DeltaM_physics3);
                    }
                  }
                  if (mBC_sideband) h_MBC_physics3->Fill(MBC_physics3);
                  if (DsMass_sideband) h_dsPlusM_physics3->Fill(dsPlusM_physics3);
                }
              }
            }
          }
        }
      }
    }
    }
  }
  
  std::cout<<"Number of generic candidates = "<<ngenericEvents<<std::endl;
  for (DecayMap::iterator i_decay=decayFrequency_generic.begin(); i_decay!=decayFrequency_generic.end(); ++i_decay)
  {
    std::cout<<"generic decay mode "<<i_decay->first<<" had "<<i_decay->second<<" decays"<<std::endl;
  }
  std::cout<<"Number of generic events after ntuplizer = "<<genericEvents.noCut<<std::endl;
  std::cout<<"Number of generic events after tag = "<<genericEvents.tagCut<<std::endl;
  std::cout<<"Number of generic events after dsPlusMCut = "<<genericEvents.dsPlusMCut<<std::endl;
  std::cout<<"Number of generic events after deltaECut = "<<genericEvents.deltaECut<<std::endl;
  std::cout<<"Number of generic events after mbcCut = "<<genericEvents.mbcCut<<std::endl;
  std::cout<<"Number of generic events after deltaMCut = "<<genericEvents.deltaMCut<<std::endl;
  std::cout<<"Number of generic events after diffD0Cut = "<<genericEvents.diffD0Cut<<std::endl;
  std::cout<<"Number of generic events after dPhiCut = "<<genericEvents.dPhiCut<<std::endl;
  std::cout<<"Number of generic events after pi0MassCut = "<<genericEvents.pi0MassCut<<std::endl;
  std::cout<<"Number of D_s- reconstructions remaining in generic = "<<oppRecon_generic<<std::endl;
  
  std::cout<<"Number of continuum candidates = "<<ncontinuEvents<<std::endl;
  for (DecayMap::iterator i_decay=decayFrequency_continu.begin(); i_decay!=decayFrequency_continu.end(); ++i_decay)
  {
    std::cout<<"continuum decay mode "<<i_decay->first<<" had "<<i_decay->second<<" decays"<<std::endl;
  }
  std::cout<<"Number of continu events after ntuplizer = "<<continuEvents.noCut<<std::endl;
  std::cout<<"Number of continu events after tag = "<<continuEvents.tagCut<<std::endl;
  std::cout<<"Number of continu events after dsPlusMCut = "<<continuEvents.dsPlusMCut<<std::endl;
  std::cout<<"Number of continu events after deltaECut = "<<continuEvents.deltaECut<<std::endl;
  std::cout<<"Number of continu events after mbcCut = "<<continuEvents.mbcCut<<std::endl;
  std::cout<<"Number of continu events after deltaMCut = "<<continuEvents.deltaMCut<<std::endl;
  std::cout<<"Number of continu events after diffD0Cut = "<<continuEvents.diffD0Cut<<std::endl;
  std::cout<<"Number of continu events after dPhiCut = "<<continuEvents.dPhiCut<<std::endl;
  std::cout<<"Number of continu events after pi0MassCut = "<<continuEvents.pi0MassCut<<std::endl;
  std::cout<<"Number of D_s- reconstructions remaining in continu = "<<oppRecon_continu<<std::endl;
  
  std::cout<<"Number of data candidates = "<<nphysics1Events<<std::endl;
  for (DecayMap::iterator i_decay=decayFrequency_physics1.begin(); i_decay!=decayFrequency_physics1.end(); ++i_decay)
  {
    std::cout<<"physics1sion decay mode "<<i_decay->first<<" had "<<i_decay->second<<" decays"<<std::endl;
  }
  std::cout<<"Number of physics1 events after ntuplizer = "<<physics1Events.noCut<<std::endl;
  std::cout<<"Number of physics1 events after tag = "<<physics1Events.tagCut<<std::endl;
  std::cout<<"Number of physics1 events after dsPlusMCut = "<<physics1Events.dsPlusMCut<<std::endl;
  std::cout<<"Number of physics1 events after deltaECut = "<<physics1Events.deltaECut<<std::endl;
  std::cout<<"Number of physics1 events after mbcCut = "<<physics1Events.mbcCut<<std::endl;
  std::cout<<"Number of physics1 events after deltaMCut = "<<physics1Events.deltaMCut<<std::endl;
  std::cout<<"Number of physics1 events after diffD0Cut = "<<physics1Events.diffD0Cut<<std::endl;
  std::cout<<"Number of physics1 events after dPhiCut = "<<physics1Events.dPhiCut<<std::endl;
  std::cout<<"Number of physics1 events after pi0MassCut = "<<physics1Events.pi0MassCut<<std::endl;
  std::cout<<"Number of D_s- reconstructions remaining in physics1 = "<<oppRecon_physics1<<std::endl;
  
  std::cout<<"Number of data candidates = "<<nphysics3Events<<std::endl;
  for (DecayMap::iterator i_decay=decayFrequency_physics3.begin(); i_decay!=decayFrequency_physics3.end(); ++i_decay)
  {
    std::cout<<"physics3sion decay mode "<<i_decay->first<<" had "<<i_decay->second<<" decays"<<std::endl;
  }
  std::cout<<"Number of physics3 events after ntuplizer = "<<physics3Events.noCut<<std::endl;
  std::cout<<"Number of physics3 events after tag = "<<physics3Events.tagCut<<std::endl;
  std::cout<<"Number of physics3 events after dsPlusMCut = "<<physics3Events.dsPlusMCut<<std::endl;
  std::cout<<"Number of physics3 events after deltaECut = "<<physics3Events.deltaECut<<std::endl;
  std::cout<<"Number of physics3 events after mbcCut = "<<physics3Events.mbcCut<<std::endl;
  std::cout<<"Number of physics3 events after deltaMCut = "<<physics3Events.deltaMCut<<std::endl;
  std::cout<<"Number of physics3 events after diffD0Cut = "<<physics3Events.diffD0Cut<<std::endl;
  std::cout<<"Number of physics3 events after dPhiCut = "<<physics3Events.dPhiCut<<std::endl;
  std::cout<<"Number of physics3 events after pi0MassCut = "<<physics3Events.pi0MassCut<<std::endl;
  std::cout<<"Number of D_s- reconstructions remaining in physics3 = "<<oppRecon_physics3<<std::endl;
  /*
  TLine *line;
  float xmin, xmax, ymin, ymax;
  
  gROOT->SetStyle("Plain");
  
  TCanvas *dsPlusM = new TCanvas("dsPlusM", "", 300, 1200);  
  xmin=dsPlusMCut_center-dsPlusMCut_range;
  xmax=dsPlusMCut_center+dsPlusMCut_range;
  dsPlusM->Divide(1,4);
  dsPlusM->cd(1);
  h_dsPlusM_generic->Draw();
  ymax=(h_dsPlusM_generic->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print(dsPlusM_generic_fileName.c_str()); 
  dsPlusM->cd(2);
  h_dsPlusM_continu->Draw();
  ymax=(h_dsPlusM_continu->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print(dsPlusM_continu_fileName.c_str());
  dsPlusM->cd(3);
  h_dsPlusM_physics1->Draw();
  ymax=(h_dsPlusM_physics1->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print(dsPlusM_physics1_fileName.c_str());
  dsPlusM->cd(4);
  h_dsPlusM_physics3->Draw();
  ymax=(h_dsPlusM_physics3->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print(dsPlusM_physics3_fileName.c_str());
  
  TCanvas *MBC = new TCanvas("MBC", "", 300, 1200);
  xmin=mbcCut_center-mbcCut_range;
  xmax=mbcCut_center+mbcCut_range;
  MBC->Divide(1,4);
  MBC->cd(1);
  h_MBC_generic->Draw();
  ymax=(h_MBC_generic->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print(MBC_generic_fileName.c_str());
  MBC->cd(2);
  h_MBC_continu->Draw();
  ymax=(h_MBC_continu->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print(MBC_continu_fileName.c_str());
  MBC->cd(3);
  h_MBC_physics1->Draw();
  ymax=(h_MBC_physics1->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print(MBC_physics1_fileName.c_str());
  MBC->cd(4);
  h_MBC_physics3->Draw();
  ymax=(h_MBC_physics3->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print(MBC_physics3_fileName.c_str());
  
  TLine *line;
  float xmin, xmax, ymin, ymax;
  
  gROOT->SetStyle("Plain");
  
  TCanvas *DeltaM = new TCanvas("DeltaM", "", 300, 1200);
  xmin=deltaMCut_center-deltaMCut_range;
  xmax=deltaMCut_center+deltaMCut_range;
  DeltaM->Divide(1,4);
  DeltaM->cd(1);
  h_DeltaM_generic->Draw();
  ymax=(h_DeltaM_generic->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print(DeltaM_generic_fileName.c_str());
  DeltaM->cd(2);
  h_DeltaM_continu->Draw();
  ymax=(h_DeltaM_continu->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print(DeltaM_continu_fileName.c_str());
  DeltaM->cd(3);
  h_DeltaM_physics1->Draw();
  ymax=(h_DeltaM_physics1->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print(DeltaM_physics1_fileName.c_str());
  DeltaM->cd(4);
  h_DeltaM_physics3->Draw();
  ymax=(h_DeltaM_physics3->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print(DeltaM_physics3_fileName.c_str());
  
  TCanvas *diffD0 = new TCanvas("diffD0", "", 300, 1200);
  diffD0->Divide(1,4);
  diffD0->cd(1);
  h_diffD0_generic->Draw();
  ymax=(h_diffD0_generic->GetMaximum())*0.75;
  line = new TLine(diffD0Cut,ymax,diffD0Cut,0); line->Draw();
  diffD0->cd(2);
  h_diffD0_continu->Draw();
  ymax=(h_diffD0_continu->GetMaximum())*0.75;
  line = new TLine(diffD0Cut,ymax,diffD0Cut,0); line->Draw();
  diffD0->cd(3);
  h_diffD0_physics1->Draw();
  ymax=(h_diffD0_physics1->GetMaximum())*0.75;
  line = new TLine(diffD0Cut,ymax,diffD0Cut,0); line->Draw();
  diffD0->cd(4);
  h_diffD0_physics3->Draw();
  ymax=(h_diffD0_physics3->GetMaximum())*0.75;
  line = new TLine(diffD0Cut,ymax,diffD0Cut,0); line->Draw();
  
  TCanvas *dPhi_diffD0 = new TCanvas ("dPhi_diffD0", "", 300, 1200);
  xmin=-2; xmax=dPhiCutLess;
  ymin=diffD0Cut; ymax=0.01;
  dPhi_diffD0->Divide(1,4);
  dPhi_diffD0->cd(1);
  h_dPhi_diffD0_generic->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  dPhi_diffD0->cd(2);
  h_dPhi_diffD0_continu->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  dPhi_diffD0->cd(3);
  h_dPhi_diffD0_physics1->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  dPhi_diffD0->cd(4);
  h_dPhi_diffD0_physics3->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  
  TCanvas *dPhi = new TCanvas("dPhi", "", 300, 1200);
  dPhi->Divide(1,4);
  dPhi->cd(1);
  h_dPhi_generic->Draw();
  line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
  gPad->Print(dPhi_generic_fileName.c_str());
  dPhi->cd(2);  
  h_dPhi_continu->Draw();
  ymax=(h_dPhi_continu->GetMaximum())*0.75;
  line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
  gPad->Print(dPhi_continu_fileName.c_str());
  dPhi->cd(3);
  h_dPhi_physics1->Draw();
  ymax=(h_dPhi_physics1->GetMaximum())*0.75;
  line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
  gPad->Print(dPhi_physics1_fileName.c_str());
  dPhi->cd(4);
  h_dPhi_physics3->Draw();
  ymax=(h_dPhi_physics3->GetMaximum())*0.75;
  line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
  gPad->Print(dPhi_physics3_fileName.c_str());
  
  TCanvas *pi0 = new TCanvas("pi0", "", 300, 1200);
  xmin=pi0MassCut_center-pi0MassCut_range;
  xmax=pi0MassCut_center+pi0MassCut_range;
  pi0->Divide(1,4);
  pi0->cd(1);
  h_pi0_generic->Draw();
  ymax=(h_pi0_generic->GetMaximum())*0.75;
  line = new TLine(xmin, ymax, xmin, 0); line->Draw();
  line = new TLine(xmax, ymax, xmax, 0); line->Draw();
  gPad->Print(pi0_generic_fileName.c_str());
  pi0->cd(2);
  h_pi0_continu->Draw();
  ymax=(h_pi0_continu->GetMaximum())*0.75;
  line = new TLine(xmin, ymax, xmin, 0); line->Draw();
  line = new TLine(xmax, ymax, xmax, 0); line->Draw();
  gPad->Print(pi0_continu_fileName.c_str());
  pi0->cd(3);
  h_pi0_physics1->Draw();
  ymax=(h_pi0_physics1->GetMaximum())*0.75;
  line = new TLine(xmin, ymax, xmin, 0); line->Draw();
  line = new TLine(xmax, ymax, xmax, 0); line->Draw();
  gPad->Print(pi0_physics1_fileName.c_str());
  pi0->cd(4);
  h_pi0_physics3->Draw();
  ymax=(h_pi0_physics3->GetMaximum())*0.75;
  line = new TLine(xmin, ymax, xmin, 0); line->Draw();
  line = new TLine(xmax, ymax, xmax, 0); line->Draw();
  gPad->Print(pi0_physics3_fileName.c_str());
  
  TCanvas *electron = new TCanvas("electron");
  h_electronE->Draw();
  */
  
  std::cout<<"Sideband summary for "<<decay<<std::endl;
  
  if (mBC_sideband)
  {
    std::cout<<"=== mBC sideband ==="<<std::endl;
    std::cout<<"Generic MC: ";
    double generic_beforeSignal=h_MBC_generic->Integral(h_MBC_generic->GetXaxis()->GetFirst(), h_MBC_generic->GetXaxis()->FindBin(mbcCut_center-mbcCut_range));
    double generic_afterSignal=h_MBC_generic->Integral(h_MBC_generic->GetXaxis()->FindBin(mbcCut_center+mbcCut_range), h_MBC_generic->GetXaxis()->GetLast());
    std::cout<<" before signal/20 = "<<generic_beforeSignal/20;
    std::cout<<", after signal/20 = "<<generic_afterSignal/20<<std::endl;
    std::cout<<"Continuum MC: ";
    double continu_beforeSignal=h_MBC_continu->Integral(h_MBC_continu->GetXaxis()->GetFirst(), h_MBC_continu->GetXaxis()->FindBin(mbcCut_center-mbcCut_range));
    double continu_afterSignal=h_MBC_continu->Integral(h_MBC_continu->GetXaxis()->FindBin(mbcCut_center+mbcCut_range), h_MBC_continu->GetXaxis()->GetLast());
    std::cout<<" before signal/5 = "<<continu_beforeSignal/5;
    std::cout<<", after signal/5 = "<<continu_afterSignal/5<<std::endl;
    std::cout<<"Sum MC: ";
    std::cout<<" before signal = "<<generic_beforeSignal/20+continu_beforeSignal/5;
    std::cout<<", after signal = "<<generic_afterSignal/20+continu_afterSignal/5<<std::endl;
    std::cout<<"Data (pion-fitted): ";
    double physics1_beforeSignal=h_MBC_physics1->Integral(h_MBC_physics1->GetXaxis()->GetFirst(), h_MBC_physics1->GetXaxis()->FindBin(mbcCut_center-mbcCut_range));
    double physics1_afterSignal=h_MBC_physics1->Integral(h_MBC_physics1->GetXaxis()->FindBin(mbcCut_center+mbcCut_range), h_MBC_physics1->GetXaxis()->GetLast());
    std::cout<<" before signal = "<<physics1_beforeSignal<<" MC/Data = "<<(generic_beforeSignal/20+continu_beforeSignal/5)/physics1_beforeSignal;
    std::cout<<", after signal = "<<physics1_afterSignal<<" MC/Data = "<<(generic_afterSignal/20+continu_afterSignal/5)/physics1_afterSignal<<std::endl;
    std::cout<<"Data (electron-fitted): ";
    double physics3_beforeSignal=h_MBC_physics3->Integral(h_MBC_physics3->GetXaxis()->GetFirst(), h_MBC_physics3->GetXaxis()->FindBin(mbcCut_center-mbcCut_range));
    double physics3_afterSignal=h_MBC_physics3->Integral(h_MBC_physics3->GetXaxis()->FindBin(mbcCut_center+mbcCut_range), h_MBC_physics3->GetXaxis()->GetLast());
    std::cout<<" before signal = "<<physics3_beforeSignal<<" MC/Data = "<<(generic_beforeSignal/20+continu_beforeSignal/5)/physics3_beforeSignal;
    std::cout<<", after signal = "<<physics3_afterSignal<<" MC/Data = "<<(generic_afterSignal/20+continu_afterSignal/5)/physics3_afterSignal<<std::endl;
    std::cout<<"===================="<<std::endl;
    std::cout<<" TeX output "<<std::endl;
    std::cout.precision(2);
    std::cout<<" & "<<generic_beforeSignal/20<<" $\\pm$ "<<pow(generic_beforeSignal, 0.5)/20;
    std::cout<<" & "<<continu_beforeSignal/5<<" $\\pm$ "<<pow(continu_beforeSignal, 0.5)/5;
    double mc_beforeSignal=generic_beforeSignal/20+continu_beforeSignal/5;
    double mc_beforeSignal_error=pow(generic_beforeSignal/400+continu_beforeSignal/25, 0.5);
    std::cout<<" & "<<mc_beforeSignal<<" $\\pm$ "<<mc_beforeSignal_error;
    std::cout<<" & "<<physics1_beforeSignal<<" $\\pm$ "<<pow(physics1_beforeSignal, 0.5);
    std::cout<<" & "<<physics3_beforeSignal<<" $\\pm$ "<<pow(physics3_beforeSignal, 0.5);
    double data1MC_beforeSignal=physics1_beforeSignal/mc_beforeSignal;
    double data1MC_beforeSignal_error=data1MC_beforeSignal*pow(1/physics1_beforeSignal+pow(mc_beforeSignal_error/mc_beforeSignal, 2), 0.5);
    double data3MC_beforeSignal=physics3_beforeSignal/mc_beforeSignal;
    double data3MC_beforeSignal_error=data3MC_beforeSignal*pow(1/physics3_beforeSignal+pow(mc_beforeSignal_error/mc_beforeSignal, 2), 0.5);
    std::cout<<" & "<<data1MC_beforeSignal<<" $\\pm$ "<<data1MC_beforeSignal_error;
    std::cout<<" & "<<data3MC_beforeSignal<<" $\\pm$ "<<data3MC_beforeSignal_error;
    std::cout<<" & "<<generic_afterSignal/20<<" $\\pm$ "<<pow(generic_afterSignal, 0.5)/20;
    std::cout<<" & "<<continu_afterSignal/5<<" $\\pm$ "<<pow(continu_afterSignal, 0.5)/5;
    double mc_afterSignal=generic_afterSignal/20+continu_afterSignal/5;
    double mc_afterSignal_error=pow(generic_afterSignal/400+continu_afterSignal/25, 0.5);
    std::cout<<" & "<<mc_afterSignal<<" $\\pm$ "<<mc_afterSignal_error;
    std::cout<<" & "<<physics1_afterSignal<<" $\\pm$ "<<pow(physics1_afterSignal, 0.5);
    std::cout<<" & "<<physics3_afterSignal<<" $\\pm$ "<<pow(physics3_afterSignal, 0.5);
    double data1MC_afterSignal=physics1_afterSignal/mc_afterSignal;
    double data1MC_afterSignal_error=data1MC_afterSignal*pow(1/physics1_afterSignal+pow(mc_afterSignal_error/mc_afterSignal, 2), 0.5);
    double data3MC_afterSignal=physics3_afterSignal/mc_afterSignal;
    double data3MC_afterSignal_error=data3MC_afterSignal*pow(1/physics3_afterSignal+pow(mc_afterSignal_error/mc_afterSignal, 2), 0.5);
    std::cout<<" & "<<data1MC_afterSignal<<" $\\pm$ "<<data1MC_afterSignal_error;
    std::cout<<" & "<<data3MC_afterSignal<<" $\\pm$ "<<data3MC_afterSignal_error<<" \\\\  \\hline"<<std::endl;
  }
  
  if (DeltaM_sideband)
  {
    std::cout<<"=== DeltaM sideband ==="<<std::endl;
    std::cout<<"Generic MC: ";
    double generic_beforeSignal=h_DeltaM_generic->Integral(h_DeltaM_generic->GetXaxis()->GetFirst(), h_DeltaM_generic->GetXaxis()->FindBin(deltaMCut_center-deltaMCut_range));
    double generic_afterSignal=h_DeltaM_generic->Integral(h_DeltaM_generic->GetXaxis()->FindBin(deltaMCut_center+deltaMCut_range), h_DeltaM_generic->GetXaxis()->GetLast());
    std::cout<<" before signal/20 = "<<generic_beforeSignal/20;
    std::cout<<" after signal/20 = "<<generic_afterSignal/20<<std::endl;
    std::cout<<"Continuum MC: ";
    double continu_beforeSignal=h_DeltaM_continu->Integral(h_DeltaM_continu->GetXaxis()->GetFirst(), h_DeltaM_generic->GetXaxis()->FindBin(deltaMCut_center-deltaMCut_range));
    double continu_afterSignal=h_DeltaM_continu->Integral(h_DeltaM_continu->GetXaxis()->FindBin(deltaMCut_center+deltaMCut_range), h_DeltaM_generic->GetXaxis()->GetLast());
    std::cout<<" before signal/5 = "<<continu_beforeSignal/5;
    std::cout<<" after signal/5 = "<<continu_afterSignal/5<<std::endl;
    std::cout<<"Sum MC: ";
    std::cout<<" before signal = "<<generic_beforeSignal/20+continu_beforeSignal/5;
    std::cout<<" after signal = "<<generic_afterSignal/20+continu_afterSignal/5<<std::endl;
    std::cout<<"Data (pion-fitted): ";
    double physics1_beforeSignal=h_DeltaM_physics1->Integral(h_DeltaM_physics1->GetXaxis()->GetFirst(), h_DeltaM_generic->GetXaxis()->FindBin(deltaMCut_center-deltaMCut_range));
    double physics1_afterSignal=h_DeltaM_physics1->Integral(h_DeltaM_physics1->GetXaxis()->FindBin(deltaMCut_center+deltaMCut_range), h_DeltaM_generic->GetXaxis()->GetLast());
    std::cout<<" before signal = "<<physics1_beforeSignal;
    std::cout<<" after signal = "<<physics1_afterSignal<<std::endl;
    std::cout<<"Data (electron-fitted): ";
    double physics3_beforeSignal=h_DeltaM_physics3->Integral(h_DeltaM_physics3->GetXaxis()->GetFirst(), h_DeltaM_generic->GetXaxis()->FindBin(deltaMCut_center-deltaMCut_range));
    double physics3_afterSignal=h_DeltaM_physics3->Integral(h_DeltaM_physics3->GetXaxis()->FindBin(deltaMCut_center+deltaMCut_range), h_DeltaM_generic->GetXaxis()->GetLast());
    std::cout<<" before signal = "<<physics3_beforeSignal;
    std::cout<<" after signal = "<<physics3_afterSignal<<std::endl;
    std::cout<<"===================="<<std::endl;
    std::cout<<" TeX output "<<std::endl;
    std::cout<<fixed;
    std::cout<<setprecision(1);
    std::cout<<" & "<<generic_beforeSignal/20<<" $\\pm$ "<<pow(generic_beforeSignal, 0.5)/20;
    std::cout<<" & "<<continu_beforeSignal/5<<" $\\pm$ "<<pow(continu_beforeSignal, 0.5)/5;
    double mc_beforeSignal=generic_beforeSignal/20+continu_beforeSignal/5;
    double mc_beforeSignal_error=pow(generic_beforeSignal/400+continu_beforeSignal/25, 0.5);
    std::cout<<" & "<<mc_beforeSignal<<" $\\pm$ "<<mc_beforeSignal_error;
    std::cout<<" & "<<physics1_beforeSignal<<" $\\pm$ "<<pow(physics1_beforeSignal, 0.5);
    std::cout<<" & "<<physics3_beforeSignal<<" $\\pm$ "<<pow(physics3_beforeSignal, 0.5);
    double data1MC_beforeSignal=physics1_beforeSignal/mc_beforeSignal;
    double data1MC_beforeSignal_error=data1MC_beforeSignal*pow(1/physics1_beforeSignal+pow(mc_beforeSignal_error/mc_beforeSignal, 2), 0.5);
    double data3MC_beforeSignal=physics3_beforeSignal/mc_beforeSignal;
    double data3MC_beforeSignal_error=data3MC_beforeSignal*pow(1/physics3_beforeSignal+pow(mc_beforeSignal_error/mc_beforeSignal, 2), 0.5);
    std::cout<<setprecision(2);
    std::cout<<" & "<<data1MC_beforeSignal<<" $\\pm$ "<<data1MC_beforeSignal_error;
    std::cout<<" & "<<data3MC_beforeSignal<<" $\\pm$ "<<data3MC_beforeSignal_error;
    std::cout<<setprecision(1);
    std::cout<<" & "<<generic_afterSignal/20<<" $\\pm$ "<<pow(generic_afterSignal, 0.5)/20;
    std::cout<<" & "<<continu_afterSignal/5<<" $\\pm$ "<<pow(continu_afterSignal, 0.5)/5;
    double mc_afterSignal=generic_afterSignal/20+continu_afterSignal/5;
    double mc_afterSignal_error=pow(generic_afterSignal/400+continu_afterSignal/25, 0.5);
    std::cout<<" & "<<mc_afterSignal<<" $\\pm$ "<<mc_afterSignal_error;
    std::cout<<" & "<<physics1_afterSignal<<" $\\pm$ "<<pow(physics1_afterSignal, 0.5);
    std::cout<<" & "<<physics3_afterSignal<<" $\\pm$ "<<pow(physics3_afterSignal, 0.5);
    double data1MC_afterSignal=physics1_afterSignal/mc_afterSignal;
    double data1MC_afterSignal_error=data1MC_afterSignal*pow(1/physics1_afterSignal+pow(mc_afterSignal_error/mc_afterSignal, 2), 0.5);
    double data3MC_afterSignal=physics3_afterSignal/mc_afterSignal;
    double data3MC_afterSignal_error=data3MC_afterSignal*pow(1/physics3_afterSignal+pow(mc_afterSignal_error/mc_afterSignal, 2), 0.5);
    std::cout<<setprecision(2);
    std::cout<<" & "<<data1MC_afterSignal<<" $\\pm$ "<<data1MC_afterSignal_error;
    std::cout<<" & "<<data3MC_afterSignal<<" $\\pm$ "<<data3MC_afterSignal_error<<" \\\\  \\hline"<<std::endl;
  }
  
  if (DsMass_sideband)
  {
    std::cout<<"=== Ds Mass Sideband ==="<<std::endl;
    std::cout<<"Generic MC: ";
    double generic_beforeSignal=h_dsPlusM_generic->Integral(h_dsPlusM_generic->GetXaxis()->GetFirst(), h_dsPlusM_generic->GetXaxis()->FindBin(dsPlusMCut_center-dsPlusMCut_range));
    double generic_afterSignal=h_dsPlusM_generic->Integral(h_dsPlusM_generic->GetXaxis()->FindBin(dsPlusMCut_center+dsPlusMCut_range), h_dsPlusM_generic->GetXaxis()->GetLast());
    std::cout<<" before signal/20 = "<<generic_beforeSignal/20;
    std::cout<<" after signal/20 = "<<generic_afterSignal/20<<std::endl;
    std::cout<<"Continuum MC: ";
    double continu_beforeSignal=h_dsPlusM_continu->Integral(h_dsPlusM_continu->GetXaxis()->GetFirst(), h_dsPlusM_continu->GetXaxis()->FindBin(dsPlusMCut_center-dsPlusMCut_range));
    double continu_afterSignal=h_dsPlusM_continu->Integral(h_dsPlusM_continu->GetXaxis()->FindBin(dsPlusMCut_center+dsPlusMCut_range), h_dsPlusM_continu->GetXaxis()->GetLast());
    std::cout<<" before signal/5 = "<<continu_beforeSignal/5;
    std::cout<<" after signal/5 = "<<continu_afterSignal/5<<std::endl;
    std::cout<<"Data (pion-fitted): ";
    double physics1_beforeSignal=h_dsPlusM_physics1->Integral(h_dsPlusM_physics1->GetXaxis()->GetFirst(), h_dsPlusM_physics1->GetXaxis()->FindBin(dsPlusMCut_center-dsPlusMCut_range));
    double physics1_afterSignal=h_dsPlusM_physics1->Integral(h_dsPlusM_physics1->GetXaxis()->FindBin(dsPlusMCut_center+dsPlusMCut_range), h_dsPlusM_physics1->GetXaxis()->GetLast());
    std::cout<<" before signal = "<<physics1_beforeSignal;
    std::cout<<" after signal = "<<physics1_afterSignal<<std::endl;
    std::cout<<"Data (electron-fitted): ";
    double physics3_beforeSignal=h_dsPlusM_physics3->Integral(h_dsPlusM_physics3->GetXaxis()->GetFirst(), h_dsPlusM_physics3->GetXaxis()->FindBin(dsPlusMCut_center-dsPlusMCut_range));
    double physics3_afterSignal=h_dsPlusM_physics3->Integral(h_dsPlusM_physics3->GetXaxis()->FindBin(dsPlusMCut_center+dsPlusMCut_range), h_dsPlusM_physics3->GetXaxis()->GetLast());
    std::cout<<" before signal = "<<physics3_beforeSignal;
    std::cout<<" after signal = "<<physics3_afterSignal<<std::endl;
    std::cout<<"===================="<<std::endl;
    std::cout<<" TeX output "<<std::endl;
    std::cout<<fixed;
    std::cout<<setprecision(1);
    std::cout<<" & "<<generic_beforeSignal/20<<" $\\pm$ "<<pow(generic_beforeSignal, 0.5)/20;
    std::cout<<" & "<<continu_beforeSignal/5<<" $\\pm$ "<<pow(continu_beforeSignal, 0.5)/5;
    double mc_beforeSignal=generic_beforeSignal/20+continu_beforeSignal/5;
    double mc_beforeSignal_error=pow(generic_beforeSignal/400+continu_beforeSignal/25, 0.5);
    std::cout<<" & "<<mc_beforeSignal<<" $\\pm$ "<<mc_beforeSignal_error;
    std::cout<<" & "<<physics1_beforeSignal<<" $\\pm$ "<<pow(physics1_beforeSignal, 0.5);
    std::cout<<" & "<<physics3_beforeSignal<<" $\\pm$ "<<pow(physics3_beforeSignal, 0.5);
    double data1MC_beforeSignal=physics1_beforeSignal/mc_beforeSignal;
    double data1MC_beforeSignal_error=data1MC_beforeSignal*pow(1/physics1_beforeSignal+pow(mc_beforeSignal_error/mc_beforeSignal, 2), 0.5);
    double data3MC_beforeSignal=physics3_beforeSignal/mc_beforeSignal;
    double data3MC_beforeSignal_error=data3MC_beforeSignal*pow(1/physics3_beforeSignal+pow(mc_beforeSignal_error/mc_beforeSignal, 2), 0.5);
    std::cout<<setprecision(2);
    std::cout<<" & "<<data1MC_beforeSignal<<" $\\pm$ "<<data1MC_beforeSignal_error;
    std::cout<<" & "<<data3MC_beforeSignal<<" $\\pm$ "<<data3MC_beforeSignal_error;
    std::cout<<setprecision(1);
    std::cout<<" & "<<generic_afterSignal/20<<" $\\pm$ "<<pow(generic_afterSignal, 0.5)/20;
    std::cout<<" & "<<continu_afterSignal/5<<" $\\pm$ "<<pow(continu_afterSignal, 0.5)/5;
    double mc_afterSignal=generic_afterSignal/20+continu_afterSignal/5;
    double mc_afterSignal_error=pow(generic_afterSignal/400+continu_afterSignal/25, 0.5);
    std::cout<<" & "<<mc_afterSignal<<" $\\pm$ "<<mc_afterSignal_error;
    std::cout<<" & "<<physics1_afterSignal<<" $\\pm$ "<<pow(physics1_afterSignal, 0.5);
    std::cout<<" & "<<physics3_afterSignal<<" $\\pm$ "<<pow(physics3_afterSignal, 0.5);
    double data1MC_afterSignal=physics1_afterSignal/mc_afterSignal;
    double data1MC_afterSignal_error=data1MC_afterSignal*pow(1/physics1_afterSignal+pow(mc_afterSignal_error/mc_afterSignal, 2), 0.5);
    double data3MC_afterSignal=physics3_afterSignal/mc_afterSignal;
    double data3MC_afterSignal_error=data3MC_afterSignal*pow(1/physics3_afterSignal+pow(mc_afterSignal_error/mc_afterSignal, 2), 0.5);
    std::cout<<setprecision(2);
    std::cout<<" & "<<data1MC_afterSignal<<" $\\pm$ "<<data1MC_afterSignal_error;
    std::cout<<" & "<<data3MC_afterSignal<<" $\\pm$ "<<data3MC_afterSignal_error<<" \\\\  \\hline"<<std::endl;
  }
  
  if (!mBC_sideband && !DeltaM_sideband && !DsMass_sideband)
  {
    double generic_error=pow(genericEvents.pi0MassCut, 0.5)/20;
    double continu_error=pow(continuEvents.pi0MassCut, 0.5)/5;
    double physics1_error=pow(physics1Events.pi0MassCut, 0.5);
    double physics3_error=pow(physics3Events.pi0MassCut, 0.5);
    double mc=genericEvents.pi0MassCut/20+continuEvents.pi0MassCut/5;
    double mc_error=pow(genericEvents.pi0MassCut/400+continuEvents.pi0MassCut/25, 0.5);
    
    std::cout<<"=== Signal Region ==="<<std::endl;
    std::cout<<" & "<<genericEvents.pi0MassCut/20<<" $\\pm$ "<<generic_error;
    std::cout<<" & "<<continuEvents.pi0MassCut/5<<" $\\pm$ "<<continu_error;
    std::cout<<" & "<<mc<<" $\\pm$ "<<mc_error;
    std::cout<<" & "<<physics1Events.pi0MassCut<<" $\\pm$ "<<physics1_error;
    std::cout<<" & "<<physics3Events.pi0MassCut<<" $\\pm$ "<<physics3_error<<std::endl;
  }
  
  return 0;
}
