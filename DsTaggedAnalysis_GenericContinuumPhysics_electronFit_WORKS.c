#include <iostream>
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TROOT.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"
#include "TFile.h"
#include "TChain.h"
#include "TPaveStats.h"
#include <set>
#include <map>
#include <iomanip>

bool DeltaM_sideband=false;
bool mBC_sideband=true;
bool DsMass_sideband=false;
bool diffD0_converSide=false;
bool dPhi_sideband=false;
double sideband_size=4;
double sideband_dist=0.5;

double lowerSide_xmin;
double lowerSide_xmax;
double xmin;
double xmax;
double upperSide_xmin;
double upperSide_xmax;

std::string decay="KKpi";
// options: KKpi, KsK, pieta, pietaprime, KKpipi0, pipipi, KsKmpipi, pipi0eta, pietaprimerho
std::string decay_tex;

double pi=3.14159265358979;
float decayNumber;
double dsPlusMCut_center, dsPlusMCut_range;
double mbcCut_center, mbcCut_range;
double deltaMCut_center=0.150, deltaMCut_range=0.013;
double diffD0Cut;
double dPhiCutLess;
double pi0MassCut_center=0.135, pi0MassCut_range=0.0;
double electronEnergyThreshold=0.15;
double branchingFr_mode;

double nSignalSample=9988*2;
double nConversionSample_Dsp=4511222; //3753305 in backup file;
double nConversionSample_Dsm=4896941; //507839;
double genericScale=1/19.2;
double genericScale_error=0.8/(19.2*19.2); // 1/(19.2 +- 1.5) = 0.0521 +- 0.0041
double continuScale=1/5.;
double continuScale_error=0.0;

double luminosity=586; // /pb
double luminosity_error=6;
double prodCrossSection_DsDss=948;
double prodCrossSection_DsDss_error=36;
double branchingFr_Dsstgamma=0.942;
double branchingFr_Dsstgamma_error=0.007;
double branchingFr_signal=branchingFr_Dsstgamma*0.0065;

double deltam_x1, deltam_x2, deltam_y1, deltam_y2;
double mbc_x1, mbc_x2, mbc_y1, mbc_y2;

std::string signalFileName_Dsp="/nfs/cor/an2/souvik/MC_vtosll_Dsp_";
std::string signalFileName_Dsm="/nfs/cor/an2/souvik/MC_vtosll_Dsm_";
std::string converFileName_Dsp="/nfs/cor/an3/souvik/MC_gamma_Dsp_generic/DsTaggedDecaysProc_MC_gamma_Dsp_generic.root";
std::string converFileName_Dsm="/nfs/cor/an3/souvik/MC_gamma_Dsm_generic/DsTaggedDecaysProc_MC_gamma_Dsm_generic.root";

std::string itoa (int value)
{
	char buffer[50];
	sprintf(buffer,"%d",value);
	std::string str(buffer);
	return str;
}


void setValues()
{
  if (decay=="KKpi")
  {
    decayNumber=401;
    decay_tex="K^{+} K^{-} #pi^{#pm}";
    
    // Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.011;
    mbcCut_center=2.112; mbcCut_range=0.004;
    deltaMCut_center=0.1438; deltaMCut_range=0.006;
    diffD0Cut=-0.006;
    dPhiCutLess=0.1;
    
    branchingFr_mode=0.055;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
  }
  else if (decay=="KsK")
  {
    decayNumber=400;
    decay_tex="K_{S} K^{#pm}";
    
    // Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.008;
    mbcCut_center=2.112; mbcCut_range=0.007;
    deltaMCut_center=0.1438; deltaMCut_range=0.006;
    diffD0Cut=-0.004;
    dPhiCutLess=0.14;
    
    branchingFr_mode=0.0149;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.35; mbc_y2=0.6;
  }
  else if (decay=="pieta")
  {
    decayNumber=440;
    decay_tex="#pi^{#pm} #eta, #eta #rightarrow #gamma #gamma";
    
    //Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.016;
    mbcCut_center=2.112; mbcCut_range=0.008;
    deltaMCut_center=0.1438; deltaMCut_range=0.008;
    diffD0Cut=-0.004;
    dPhiCutLess=0.12;
    
    branchingFr_mode=0.0158*0.3931;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
  }
  else if (decay=="pietaprime")
  {
    decayNumber=460;
    decay_tex="#pi^{#pm} #eta', #eta' #rightarrow #pi^{+} #pi^{-} #eta, #eta #rightarrow #gamma #gamma";
    
    //Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.011;
    mbcCut_center=2.112; mbcCut_range=0.011;
    deltaMCut_center=0.1438; deltaMCut_range=0.013;
    diffD0Cut=-0.003;
    dPhiCutLess=0.11;
    
    branchingFr_mode=0.038*0.446*0.3931;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
  }
  else if (decay=="KKpipi0")
  {
    decayNumber=404;
    decay_tex="K^{+} K^{-} #pi^{#pm} #pi^{0}";
    
    //Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.01;
    mbcCut_center=2.112; mbcCut_range=0.004;
    deltaMCut_center=0.1438; deltaMCut_range=0.006;
    diffD0Cut=-0.006;
    dPhiCutLess=0.12;
    
    branchingFr_mode=0.056;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
  }
  else if (decay=="pipipi")
  {
    decayNumber=421;
    decay_tex="#pi^{+} #pi^{-} #pi^{#pm}";
    
    //Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.012;
    mbcCut_center=2.112; mbcCut_range=0.004;
    deltaMCut_center=0.1438; deltaMCut_range=0.006;
    diffD0Cut=-0.006;
    dPhiCutLess=0.1;
    
    branchingFr_mode=0.0111;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
  }
  else if (decay=="KsKmpipi")
  {
    decayNumber=406;
    decay_tex="K^{*+}K^{*0}";
    
    //Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.006;
    mbcCut_center=2.112; mbcCut_range=0.005;
    deltaMCut_center=0.1438; deltaMCut_range=0.008;
    diffD0Cut=-0.005;
    dPhiCutLess=0.13;
    
    branchingFr_mode=0.0164;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
  }
  else if (decay=="pipi0eta")
  {
    decayNumber=441;
    decay_tex="#eta #rho^{#pm}, #eta #rightarrow #gamma #gamma, #rho^{#pm} #rightarrow #pi^{#pm} #pi^{0}";
    
    //Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.015;
    mbcCut_center=2.112; mbcCut_range=0.004;
    deltaMCut_center=0.1438; deltaMCut_range=0.005;
    diffD0Cut=-0.007;
    dPhiCutLess=0.13;
    
    branchingFr_mode=0.130*1.*0.3931;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
  }
  else if (decay=="pietaprimerho")
  {
    decayNumber=480;
    decay_tex="#pi^{#pm} #eta', #eta' #rightarrow #rho^{0} #gamma";
    
    //Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.012;
    mbcCut_center=2.112; mbcCut_range=0.004;
    deltaMCut_center=0.1438; deltaMCut_range=0.007;
    diffD0Cut=-0.006;
    dPhiCutLess=0.11;
    
    branchingFr_mode=0.038*0.294;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
  }
  
  mbcCut_range*=1;
  diffD0Cut=-0.005;
  dPhiCutLess=0.12;
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

void setBounds()
{
  if (DeltaM_sideband)
  {
    lowerSide_xmin=0.1;
    lowerSide_xmax=deltaMCut_center-(1+2*sideband_dist)*deltaMCut_range;
    xmin=deltaMCut_center-deltaMCut_range;
    xmax=deltaMCut_center+deltaMCut_range;
    upperSide_xmin=deltaMCut_center+(1+2*sideband_dist)*deltaMCut_range;
    upperSide_xmax=0.25;
  }
  if (mBC_sideband)
  {
    lowerSide_xmin=2.06;
    lowerSide_xmax=mbcCut_center-(1+2*sideband_dist)*mbcCut_range;
    xmin=mbcCut_center-mbcCut_range;
    xmax=mbcCut_center+mbcCut_range;
    upperSide_xmin=mbcCut_center+(1+2*sideband_dist)*mbcCut_range;
    upperSide_xmax=2.155;
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

bool pi0MassCut(float pi0Mass)
{
  return (fabs(pi0Mass-pi0MassCut_center)>pi0MassCut_range);
}

bool electronEnergyCut(double electronEnergy, double positronEnergy)
{
  return (electronEnergy<electronEnergyThreshold && positronEnergy<electronEnergyThreshold);
}

//Chebyshev polynomials:
double T_1(double x) {return x;}
double T_2(double x) {return 2*x*x-1;}
double T_3(double x) {return 4*x*x*x-3*x;}

double a0_mc=-2.82368e+03;
double a1_mc=6.82327e+03;
double a2_mc=-2.75178e+03;
double a3_mc=1.86979e+03;
Double_t deltaM_MCFit(Double_t *x, Double_t *par)
{
  if (x[0]>(deltaMCut_center-(1+2*sideband_dist)*deltaMCut_range) && x[0]<(deltaMCut_center+(1+2*sideband_dist)*deltaMCut_range)) TF1::RejectPoint();
  return (a0_mc+a1_mc*T_1(x[0])+a2_mc*T_2(x[0])+a3_mc*T_3(x[0]))*par[0];
}

double a0_data=2.10895e+03;
double a1_data=-8.20634e+03;
double a2_data=2.08660e+03;
double a3_data=-2.57894e+03;
Double_t deltaM_DataFit(Double_t *x, Double_t *par)
{
  if (x[0]>(deltaMCut_center-(1+2*sideband_dist)*deltaMCut_range) && x[0]<(deltaMCut_center+(1+2*sideband_dist)*deltaMCut_range)) TF1::RejectPoint();
  return (a0_data+a1_data*T_1(x[0])+a2_data*T_2(x[0])+a3_data*T_3(x[0]))*par[0];
}

double b0_mc=-4.56878e+02;
double b1_mc=2.23446e+02;
Double_t mbc_MCFit(Double_t *x, Double_t *par)
{
  if (x[0]>(mbcCut_center-(1+2*sideband_dist)*mbcCut_range) && x[0]<(mbcCut_center+(1+2*sideband_dist)*mbcCut_range)) TF1::RejectPoint();
  return (b0_mc+b1_mc*x[0])*pow(2.155-x[0], 0.5)*par[0];
}

double b0_physics=-3.25895e+02;
double b1_physics=1.60872e+02;
Double_t mbc_DataFit(Double_t *x, Double_t *par)
{
  if (x[0]>(mbcCut_center-(1+2*sideband_dist)*mbcCut_range) && x[0]<(mbcCut_center+(1+2*sideband_dist)*mbcCut_range)) TF1::RejectPoint();
  return (b0_physics+b1_physics*x[0])*pow(2.155-x[0], 0.5)*par[0];
}

Double_t mbc_SignalFit(Double_t *x, Double_t *par)
{
  Double_t result;
  return result;
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
  float mbcCut_run;     float mbcCut_event;
  float deltaMCut_run;  float deltaMCut_event;
  float diffD0Cut_run;  float diffD0Cut_event;
  float dPhiCut_run;    float dPhiCut_event;
  float pi0MassCut_run; float pi0MassCut_event;
};

int DsTaggedAnalysis_GenericContinuumPhysics_electronFit()
{
  setFileNames();
  setValues();
  setBounds();
  std::cout<<"lowerSide_xmin="<<lowerSide_xmin<<std::endl;
  std::cout<<"lowerSide_xmax="<<lowerSide_xmax<<std::endl;
  std::cout<<"xmin="<<xmin<<std::endl;
  std::cout<<"xmax="<<xmax<<std::endl;
  std::cout<<"upperSide_xmin="<<upperSide_xmin<<std::endl;
  std::cout<<"upperSide_xmax="<<upperSide_xmax<<std::endl;
  
  TChain *signalTree=new TChain("DsTaggedDecaysProc/nt3");
  signalTree->Add(signalFileName_Dsp.c_str());
  signalTree->Add(signalFileName_Dsm.c_str());  // Make sure they're in the same ratio, or the ratio is at least known.
  
  TChain *converTree_Dsp=new TChain("DsTaggedDecaysProc/nt3");
  //converTree_Dsp->Add(converFileName_Dsp.c_str());

  TChain *converTree_Dsm=new TChain("DsTaggedDecaysProc/nt3");
  //converTree_Dsm->Add(converFileName_Dsm.c_str());

  TChain *genericTree = new TChain("DsTaggedDecaysProc/nt3");
  /*
  genericTree->Add("/nfs/cor/an3/souvik/Dataset39/DsTaggedProc_GenericMC_213586_214863.root");
  genericTree->Add("/nfs/cor/an3/souvik/Dataset40/DsTaggedProc_GenericMC_215307_217385.root");
  genericTree->Add("/nfs/cor/an3/souvik/Dataset41/DsTaggedProc_GenericMC_217687_219721.root");
  genericTree->Add("/nfs/cor/an3/souvik/Dataset47/DsTaggedProc_GenericMC_230474_232255.root");
  genericTree->Add("/nfs/cor/an3/souvik/Dataset48/DsTaggedProc_GenericMC_232264_234607.root");
  */
  TChain *continuTree = new TChain("DsTaggedDecaysProc/nt3");
  /*
  continuTree->Add("/nfs/cor/an3/souvik/Dataset39/DsTaggedProc_ContinuumMC_213586_214863.root");
  continuTree->Add("/nfs/cor/an3/souvik/Dataset40/DsTaggedProc_ContinuumMC_215307_217385.root");
  continuTree->Add("/nfs/cor/an3/souvik/Dataset41/DsTaggedProc_ContinuumMC_217687_219721.root");
  continuTree->Add("/nfs/cor/an3/souvik/Dataset47/DsTaggedProc_ContinuumMC_230474_232255.root");
  continuTree->Add("/nfs/cor/an3/souvik/Dataset48/DsTaggedProc_ContinuumMC_232264_234607.root");
  */
  TChain *physicsTree = new TChain("DsTaggedDecaysProc/nt3");
  /*
  physicsTree->Add("/nfs/cor/an3/souvik/Dataset39/DsTaggedDecaysProc_ReTaggedData_213586_214863.root");
  physicsTree->Add("/nfs/cor/an3/souvik/Dataset40/DsTaggedDecaysProc_ReTaggedData_215307_217385.root");
  physicsTree->Add("/nfs/cor/an3/souvik/Dataset41/DsTaggedDecaysProc_ReTaggedData_217687_219721.root");
  physicsTree->Add("/nfs/cor/an3/souvik/Dataset47/DsTaggedDecaysProc_ReTaggedData_230474_232255.root");
  physicsTree->Add("/nfs/cor/an3/souvik/Dataset48/DsTaggedDecaysProc_ReTaggedData_232264_234607.root");
  */
  NEvents signalEvents={0,0,0,0,0,0,0,0};
  NEvents converEvents_Dsp={0,0,0,0,0,0,0,0};
  NEvents converEvents_Dsm={0,0,0,0,0,0,0,0};
  NEvents genericEvents={0,0,0,0,0,0,0,0};
  NEvents continuEvents={0,0,0,0,0,0,0,0};
  NEvents physicsEvents={0,0,0,0,0,0,0,0};
  
  EventNumber signalNumber={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  EventNumber converNumber_Dsp={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  EventNumber converNumber_Dsm={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  EventNumber genericNumber={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  EventNumber continuNumber={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  EventNumber physicsNumber={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  float runNumber_signal, eventNumber_signal;
  float runNumber_conver_Dsp, eventNumber_conver_Dsp;
  float runNumber_conver_Dsm, eventNumber_conver_Dsm;
  float runNumber_generic, eventNumber_generic;
  float runNumber_continu, eventNumber_continu;
  float runNumber_physics, eventNumber_physics;
  float dsPlusM_signal, dsPlusCharge_signal, DeltaE_signal, MBC_signal, DeltaM_signal, decayMode_signal;
  float dsPlusM_conver_Dsp, dsPlusCharge_conver_Dsp, DeltaE_conver_Dsp, MBC_conver_Dsp, DeltaM_conver_Dsp, decayMode_conver_Dsp;
  float dsPlusM_conver_Dsm, dsPlusCharge_conver_Dsm, DeltaE_conver_Dsm, MBC_conver_Dsm, DeltaM_conver_Dsm, decayMode_conver_Dsm;
  float dsPlusM_generic, dsPlusCharge_generic, DeltaE_generic, MBC_generic, DeltaM_generic, decayMode_generic;
  float dsPlusM_continu, dsPlusCharge_continu, DeltaE_continu, MBC_continu, DeltaM_continu, decayMode_continu;
  float dsPlusM_physics, dsPlusCharge_physics, DeltaE_physics, MBC_physics, DeltaM_physics, decayMode_physics;
  float d0_e_signal, d0_p_signal, z0_e_signal, z0_p_signal, px_e_signal, py_e_signal, pz_e_signal, px_p_signal, py_p_signal, pz_p_signal, E_e_signal, E_p_signal, curv_e_signal, curv_p_signal;
  float d0_e_conver_Dsp, d0_p_conver_Dsp, z0_e_conver_Dsp, z0_p_conver_Dsp, px_e_conver_Dsp, py_e_conver_Dsp, pz_e_conver_Dsp, px_p_conver_Dsp, py_p_conver_Dsp, pz_p_conver_Dsp, E_e_conver_Dsp, E_p_conver_Dsp, curv_e_conver_Dsp, curv_p_conver_Dsp;
  float d0_e_conver_Dsm, d0_p_conver_Dsm, z0_e_conver_Dsm, z0_p_conver_Dsm, px_e_conver_Dsm, py_e_conver_Dsm, pz_e_conver_Dsm, px_p_conver_Dsm, py_p_conver_Dsm, pz_p_conver_Dsm, E_e_conver_Dsm, E_p_conver_Dsm, curv_e_conver_Dsm, curv_p_conver_Dsm;
  float d0_e_generic, d0_p_generic, z0_e_generic, z0_p_generic, px_e_generic, py_e_generic, pz_e_generic, px_p_generic, py_p_generic, pz_p_generic, E_e_generic, E_p_generic, curv_e_generic, curv_p_generic;
  float d0_e_continu, d0_p_continu, z0_e_continu, z0_p_continu, px_e_continu, py_e_continu, pz_e_continu, px_p_continu, py_p_continu, pz_p_continu, E_e_continu, E_p_continu, curv_e_continu, curv_p_continu;
  float d0_e_physics, d0_p_physics, z0_e_physics, z0_p_physics, px_e_physics, py_e_physics, pz_e_physics, px_p_physics, py_p_physics, pz_p_physics, E_e_physics, E_p_physics, curv_e_physics, curv_p_physics;
  float px_e_generic_MC, py_e_generic_MC, pz_e_generic_MC, E_e_generic_MC, px_p_generic_MC, py_p_generic_MC, pz_p_generic_MC, E_p_generic_MC;
  float px_e_continu_MC, py_e_continu_MC, pz_e_continu_MC, E_e_continu_MC, px_p_continu_MC, py_p_continu_MC, pz_p_continu_MC, E_p_continu_MC;
  float pi0Mass_signal;
  float pi0Mass_conver_Dsp;
  float pi0Mass_conver_Dsm;
  float pi0Mass_generic;
  float pi0Mass_continu;
  float pi0Mass_physics;
  float conversionBit_generic, nConversionEvents_generic_Dsp=0, nConversionEvents_generic_Dsm=0;
  float mee_signalMC;
  
  // Variables used to count sideband and signal entries
  double signal_beforeSignal=0, signal_afterSignal=0, signal_atSignal=0;
  double conver_Dsp_beforeSignal=0, conver_Dsp_afterSignal=0, conver_Dsp_atSignal=0;
  double conver_Dsm_beforeSignal=0, conver_Dsm_afterSignal=0, conver_Dsm_atSignal=0;
  double generic_beforeSignal=0, generic_afterSignal=0, generic_atSignal=0;
  double continu_beforeSignal=0, continu_afterSignal=0, continu_atSignal=0;
  double physics_beforeSignal=0, physics_afterSignal=0, physics_atSignal=0;
  
  TH1D *h_dsPlusM_signal = new TH1D("h_dsPlusM_signal", "m_{D_{S}^{+}} Signal Sample; GeV", 100, 1.9, 2.1); h_dsPlusM_signal->SetLineColor(kRed);
  TH1D *h_dsPlusM_conver_Dsp = new TH1D("h_dsPlusM_conver_Dsp", "m_{D_{S}^{+}} conver_Dsp Sample; GeV", 100, 1.9, 2.1); h_dsPlusM_conver_Dsp->SetLineColor(kRed);
  TH1D *h_dsPlusM_conver_Dsm = new TH1D("h_dsPlusM_conver_Dsm", "m_{D_{S}^{+}} conver_Dsm Sample; GeV", 100, 1.9, 2.1); h_dsPlusM_conver_Dsm->SetLineColor(kRed);
  TH1D *h_dsPlusM_generic = new TH1D("h_dsPlusM_generic", "m_{D_{S}^{+}} generic Sample; GeV", 100, 1.9, 2.1); h_dsPlusM_generic->SetLineColor(kRed);
  TH1D *h_dsPlusM_continu = new TH1D("h_dsPlusM_continu", "m_{D_{S}^{+}} continuum Background Sample; GeV", 100, 1.9, 2.1); h_dsPlusM_continu->SetLineColor(kBlue);
  TH1D *h_dsPlusM_physics = new TH1D("h_dsPlusM_physics", "m_{D_{S}^{+}} Data; GeV", 100, 1.9, 2.1); h_dsPlusM_physics->SetLineColor(kGreen);
  
  TH1D *h_MBC_signal = new TH1D("h_MBC_signal", "m_{BC} Signal Sample; GeV", 80, 2., 2.16); h_MBC_signal->SetLineColor(kCyan); h_MBC_signal->Sumw2();
  TH1D *h_MBC_conver_Dsp = new TH1D("h_MBC_conver_Dsp", "m_{BC} conver_Dsp Sample; GeV", 80, 2., 2.16); h_MBC_conver_Dsp->SetLineColor(kRed); h_MBC_conver_Dsp->Sumw2();
  TH1D *h_MBC_conver_Dsm = new TH1D("h_MBC_conver_Dsm", "m_{BC} conver_Dsm Sample; GeV", 80, 2., 2.16); h_MBC_conver_Dsm->SetLineColor(kRed); h_MBC_conver_Dsm->Sumw2();
  TH1D *h_MBC_generic = new TH1D("h_MBC_generic", "m_{BC} generic Sample; GeV", 80, 2., 2.16); h_MBC_generic->SetLineColor(kGreen); h_MBC_generic->Sumw2();
  TH1D *h_MBC_continu = new TH1D("h_MBC_continu", "m_{BC} continuum Background; GeV", 80, 2., 2.16); h_MBC_continu->SetLineColor(kBlue); h_MBC_continu->Sumw2();
  TH1D *h_MBC_physics = new TH1D("h_MBC_physics", "m_{BC} Data; GeV", 80, 2., 2.16); h_MBC_physics->SetLineColor(kMagenta); h_MBC_physics->Sumw2();
  
  TH1D *h_DeltaM_signal = new TH1D("h_DeltaM_signal", "#delta m Signal Sample; GeV", 80, 0.0, 0.4); h_DeltaM_signal->SetLineColor(kCyan); h_DeltaM_signal->Sumw2();
  TH1D *h_DeltaM_conver_Dsp = new TH1D("h_DeltaM_conver_Dsp", "#delta m conver_Dsp Sample; #deltam (GeV); Number of Events", 80, 0.0, 0.4); h_DeltaM_conver_Dsp->SetLineColor(kRed); h_DeltaM_conver_Dsp->Sumw2();
  TH1D *h_DeltaM_conver_Dsm = new TH1D("h_DeltaM_conver_Dsm", "#delta m conver_Dsm Sample; #deltam (GeV); Number of Events", 80, 0.0, 0.4); h_DeltaM_conver_Dsm->SetLineColor(kRed); h_DeltaM_conver_Dsm->Sumw2();
  TH1D *h_DeltaM_generic = new TH1D("h_DeltaM_generic", "#delta m generic Sample; #deltam (GeV); Number of Events", 80, 0.0, 0.4); h_DeltaM_generic->SetLineColor(kGreen); h_DeltaM_generic->Sumw2();
  TH1D *h_DeltaM_continu = new TH1D("h_DeltaM_continu", "#delta m continuum Background Sample; #deltam (GeV); Number of Events", 80, 0.0, 0.4); h_DeltaM_continu->SetLineColor(kBlue); h_DeltaM_continu->Sumw2();
  TH1D *h_DeltaM_physics = new TH1D("h_DeltaM_physics", "#delta m Data; #deltam (GeV); Number of Events", 80, 0.0, 0.4); h_DeltaM_physics->SetLineColor(kMagenta); h_DeltaM_physics->Sumw2();
  
  TH1D *h_diffD0_signal = new TH1D("h_diffD0_signal", "#Deltad_{0} Signal Sample; m", 50, -0.01, 0.01); h_diffD0_signal->SetLineColor(kRed);
  TH1D *h_diffD0_conver_Dsp = new TH1D("h_diffD0_conver_Dsp", "#Deltad_{0} conver_Dsp Sample; m", 50, -0.01, 0.01); h_diffD0_conver_Dsp->SetLineColor(kRed);
  TH1D *h_diffD0_conver_Dsm = new TH1D("h_diffD0_conver_Dsm", "#Deltad_{0} conver_Dsm Sample; m", 50, -0.01, 0.01); h_diffD0_conver_Dsm->SetLineColor(kRed);
  TH1D *h_diffD0_generic = new TH1D("h_diffD0_generic", "#Deltad_{0} generic Sample; m", 50, -0.01, 0.01); h_diffD0_generic->SetLineColor(kRed);
  TH1D *h_diffD0_continu = new TH1D("h_diffD0_continu", "#Deltad_{0} continuum Background Sample; m", 50, -0.01, 0.01); h_diffD0_continu->SetLineColor(kBlue);
  TH1D *h_diffD0_physics = new TH1D("h_diffD0_physics", "#Deltad_{0} Data; m", 50, -0.01, 0.01); h_diffD0_physics->SetLineColor(kGreen);
  
  TH1D *h_dPhi_signal = new TH1D("h_dPhi_signal", "#Delta#Phi Signal Sample", 50, -pi/2, pi/2); h_dPhi_signal->SetLineColor(kRed);
  TH1D *h_dPhi_conver_Dsp = new TH1D("h_dPhi_conver_Dsp", "#Delta#Phi conver_Dsp Sample", 50, -pi/2, pi/2); h_dPhi_conver_Dsp->SetLineColor(kRed);
  TH1D *h_dPhi_conver_Dsm = new TH1D("h_dPhi_conver_Dsm", "#Delta#Phi conver_Dsm Sample", 50, -pi/2, pi/2); h_dPhi_conver_Dsm->SetLineColor(kRed);
  TH1D *h_dPhi_generic = new TH1D("h_dPhi_generic", "#Delta#Phi generic Sample", 50, -pi/2, pi/2); h_dPhi_generic->SetLineColor(kRed);
  TH1D *h_dPhi_continu = new TH1D("h_dPhi_continu", "#Delta#Phi continuum Background Sample", 50, -pi/2, pi/2); h_dPhi_continu->SetLineColor(kBlue);
  TH1D *h_dPhi_physics = new TH1D("h_dPhi_physics", "#Delta#Phi Data", 50, -pi/2, pi/2); h_dPhi_physics->SetLineColor(kGreen);
  
  TH2D *h_dPhi_diffD0_conver_Dsp = new TH2D("h_dPhi_diffD0_conver_Dsp", "#Delta#Phi vs #Deltad_{0} Ds*+ conversion Sample; ; m", 50, -2., 2., 50, -0.01, 0.01); h_dPhi_diffD0_conver_Dsp->SetLineColor(kRed);
  TH2D *h_dPhi_diffD0_conver_Dsm = new TH2D("h_dPhi_diffD0_conver_Dsm", "#Delta#Phi vs #Deltad_{0} Ds*- conversion Sample; ; m", 50, -2., 2., 50, -0.01, 0.01); h_dPhi_diffD0_conver_Dsm->SetLineColor(kRed);
  TH2D *h_dPhi_diffD0_generic = new TH2D("h_dPhi_diffD0_generic", "#Delta#Phi vs #Deltad_{0} generic Sample; ; m", 50, -2., 2., 50, -0.01, 0.01); h_dPhi_diffD0_generic->SetLineColor(kRed);
  TH2D *h_dPhi_diffD0_continu = new TH2D("h_dPhi_diffD0_continu", "#Delta#Phi vs #Deltad_{0} continuum Background Sample; ; m", 50, -2., 2., 50, -0.01, 0.01); h_dPhi_diffD0_continu->SetLineColor(kBlue);
  TH2D *h_dPhi_diffD0_physics = new TH2D("h_dPhi_diffD0_physics", "#Delta#Phi vs #Deltad_{0} Data; ; m", 50, -2., 2., 50, -0.01, 0.01); h_dPhi_diffD0_physics->SetLineColor(kGreen);
  
  TH1D *h_electronE;
  if (DeltaM_sideband) h_electronE = new TH1D("h_electronE", "Electron Energy; Reco Electron Energy (GeV)", 100, 0.0, 0.5); 
  else h_electronE = new TH1D("h_electronE", "Electron Energy; Reco Electron Energy (GeV)", 100, 0.0, 0.2); 
  h_electronE->SetLineColor(kRed);
  
  TH1D *h_pi0_signal = new TH1D("h_pi0_signal", "Pion Mass; GeV", 50, 0.035, .235); h_pi0_signal->SetLineColor(kRed);
  TH1D *h_pi0_conver_Dsp = new TH1D("h_pi0_conver_Dsp", "Pion Mass; GeV", 50, 0.035, .235); h_pi0_conver_Dsp->SetLineColor(kRed);
  TH1D *h_pi0_conver_Dsm = new TH1D("h_pi0_conver_Dsm", "Pion Mass; GeV", 50, 0.035, .235); h_pi0_conver_Dsm->SetLineColor(kRed);
  TH1D *h_pi0_generic = new TH1D("h_pi0_generic", "Pion Mass", 50, 0.035, .235); h_pi0_generic->SetLineColor(kRed);
  TH1D *h_pi0_continu = new TH1D("h_pi0_continu", "Pion Mass", 50, 0.035, .235); h_pi0_continu->SetLineColor(kBlue);
  TH1D *h_pi0_physics = new TH1D("h_pi0_physics", "Pion Mass", 50, 0.035, .235); h_pi0_physics->SetLineColor(kGreen);
  
  TH1D *h_mee_signalMC = new TH1D("h_mee_signalMC", "h_mee_signalMC", 100, 0., .15);
  
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
  signalTree->SetBranchAddress("keeMass_MC", &(mee_signalMC));
  
  converTree_Dsp->SetBranchAddress("Run", &(runNumber_conver_Dsp));
  converTree_Dsp->SetBranchAddress("Event", &(eventNumber_conver_Dsp));
  converTree_Dsp->SetBranchAddress("dsPlusM", &(dsPlusM_conver_Dsp));
  converTree_Dsp->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_conver_Dsp));
  converTree_Dsp->SetBranchAddress("DecayMode", &(decayMode_conver_Dsp));
  converTree_Dsp->SetBranchAddress("DeltaE", &(DeltaE_conver_Dsp));
  converTree_Dsp->SetBranchAddress("MBC", &(MBC_conver_Dsp));
  converTree_Dsp->SetBranchAddress("DeltaM", &(DeltaM_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kElectron1D0_reco", &(d0_e_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kElectron2D0_reco", &(d0_p_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kElectron1Z0_reco", &(z0_e_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kElectron2Z0_reco", &(z0_p_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kElectron1Px_reco", &(px_e_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kElectron1Py_reco", &(py_e_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kElectron1Pz_reco", &(pz_e_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kElectron2Px_reco", &(px_p_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kElectron2Py_reco", &(py_p_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kElectron2Pz_reco", &(pz_p_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kElectron1E_reco", &(E_e_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kElectron2E_reco", &(E_p_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kElectron1Curv_reco", &(curv_e_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kElectron2Curv_reco", &(curv_p_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kPi0Mass_reco", &(pi0Mass_conver_Dsp));
  
  converTree_Dsm->SetBranchAddress("Run", &(runNumber_conver_Dsm));
  converTree_Dsm->SetBranchAddress("Event", &(eventNumber_conver_Dsm));
  converTree_Dsm->SetBranchAddress("dsPlusM", &(dsPlusM_conver_Dsm));
  converTree_Dsm->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_conver_Dsm));
  converTree_Dsm->SetBranchAddress("DecayMode", &(decayMode_conver_Dsm));
  converTree_Dsm->SetBranchAddress("DeltaE", &(DeltaE_conver_Dsm));
  converTree_Dsm->SetBranchAddress("MBC", &(MBC_conver_Dsm));
  converTree_Dsm->SetBranchAddress("DeltaM", &(DeltaM_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kElectron1D0_reco", &(d0_e_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kElectron2D0_reco", &(d0_p_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kElectron1Z0_reco", &(z0_e_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kElectron2Z0_reco", &(z0_p_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kElectron1Px_reco", &(px_e_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kElectron1Py_reco", &(py_e_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kElectron1Pz_reco", &(pz_e_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kElectron2Px_reco", &(px_p_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kElectron2Py_reco", &(py_p_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kElectron2Pz_reco", &(pz_p_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kElectron1E_reco", &(E_e_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kElectron2E_reco", &(E_p_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kElectron1Curv_reco", &(curv_e_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kElectron2Curv_reco", &(curv_p_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kPi0Mass_reco", &(pi0Mass_conver_Dsm));
  
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
  genericTree->SetBranchAddress("conversionBit", &(conversionBit_generic));
  
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
  
  physicsTree->SetBranchAddress("Run", &(runNumber_physics));
  physicsTree->SetBranchAddress("Event", &(eventNumber_physics));
  physicsTree->SetBranchAddress("dsPlusM", &(dsPlusM_physics));
  physicsTree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_physics));
  physicsTree->SetBranchAddress("DecayMode", &(decayMode_physics));
  physicsTree->SetBranchAddress("MBC", &(MBC_physics));
  physicsTree->SetBranchAddress("DeltaM", &(DeltaM_physics));
  physicsTree->SetBranchAddress("kElectron1D0_reco", &(d0_e_physics));
  physicsTree->SetBranchAddress("kElectron2D0_reco", &(d0_p_physics));
  physicsTree->SetBranchAddress("kElectron1Z0_reco", &(z0_e_physics));
  physicsTree->SetBranchAddress("kElectron2Z0_reco", &(z0_p_physics));
  physicsTree->SetBranchAddress("kElectron1Px_reco", &(px_e_physics));
  physicsTree->SetBranchAddress("kElectron1Py_reco", &(py_e_physics));
  physicsTree->SetBranchAddress("kElectron1Pz_reco", &(pz_e_physics));
  physicsTree->SetBranchAddress("kElectron2Px_reco", &(px_p_physics));
  physicsTree->SetBranchAddress("kElectron2Py_reco", &(py_p_physics));
  physicsTree->SetBranchAddress("kElectron2Pz_reco", &(pz_p_physics));
  physicsTree->SetBranchAddress("kElectron1E_reco", &(E_e_physics));
  physicsTree->SetBranchAddress("kElectron2E_reco", &(E_p_physics));
  physicsTree->SetBranchAddress("kElectron1Curv_reco", &(curv_e_physics));
  physicsTree->SetBranchAddress("kElectron2Curv_reco", &(curv_p_physics));
  physicsTree->SetBranchAddress("kPi0Mass_reco", &(pi0Mass_physics));
  
  int nRepetitions=0;
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
      h_mee_signalMC->Fill(mee_signalMC);
    }
    
    if (fabs(d0_e_signal)<0.005 && fabs(d0_p_signal)<0.005 &&
        fabs(z0_e_signal)<0.05 && fabs(z0_p_signal)<0.05 &&
        electronEnergyCut(E_e_signal, E_p_signal))
    {
      if (decayMode_signal==decayNumber)
      {
        if (!DsMass_sideband) h_dsPlusM_signal->Fill(dsPlusM_signal);
        if (DsMass_sideband || dsPlusMCut(dsPlusM_signal))
        {
          if (!mBC_sideband) h_MBC_signal->Fill(MBC_signal);
          if (mBC_sideband || MBCCut(MBC_signal))
          {
            if (!DeltaM_sideband) h_DeltaM_signal->Fill(DeltaM_signal);
            if (DeltaM_sideband || DeltaMCut(DeltaM_signal))
            {
              //h_dPhi_diffD0_signal->Fill(dPhi, d0_e_signal-d0_p_signal);
              h_diffD0_signal->Fill(d0_e_signal-d0_p_signal);
              if (diffD0_converSide ^ dD0(d0_e_signal-d0_p_signal))
              {
                if (!dPhi_sideband) h_dPhi_signal->Fill(dPhi);
                if (dPhi_sideband || dPhiCut(dPhi))
                {
                  h_pi0_signal->Fill(pi0Mass_signal);
                  if (pi0MassCut(pi0Mass_signal))
                  {
                    if (runNumber_signal!=signalNumber.pi0MassCut_run || eventNumber_signal!=signalNumber.pi0MassCut_event)
                    {
                      ++signalEvents.pi0MassCut;
                      signalNumber.pi0MassCut_run=runNumber_signal;
                      signalNumber.pi0MassCut_event=eventNumber_signal;
                    }
                    if (DsMass_sideband)
                    {
                      h_dsPlusM_signal->Fill(dsPlusM_signal);
                      if (dsPlusM_signal>lowerSide_xmin && dsPlusM_signal<=lowerSide_xmax) ++signal_beforeSignal;
                      if (dsPlusM_signal>xmin && dsPlusM_signal<xmax) ++signal_atSignal;
                      if (dsPlusM_signal>=upperSide_xmin && dsPlusM_signal<upperSide_xmax) ++signal_afterSignal;
                    }
                    if (mBC_sideband)
                    {
                      // Quick check to see if (run, event) repeats. Stealing some unused variables here.
                      if (runNumber_signal!=signalNumber.dPhiCut_run || eventNumber_signal!=signalNumber.dPhiCut_event)
                      {
                        signalNumber.dPhiCut_run=runNumber_signal;
                        signalNumber.dPhiCut_event=eventNumber_signal;
                        signalEvents.dsPlusMCut=dsPlusM_signal;
                        signalEvents.mbcCut=MBC_signal;
                        signalEvents.deltaMCut=DeltaM_signal;
                        signalEvents.diffD0Cut=d0_e_signal-d0_p_signal;
                        signalEvents.dPhiCut=dPhi;
                      }
                      else
                      {
                        std::cout<<"Repetition in (run, event) = ("<<runNumber_signal<<", "<<eventNumber_signal<<")";
                        std::cout<<" dsPlusM = "<<signalEvents.dsPlusMCut;
                        std::cout<<" mbcCut = "<<signalEvents.mbcCut;
                        std::cout<<" deltaMCut = "<<signalEvents.deltaMCut;
                        std::cout<<" diffD0Cut = "<<signalEvents.diffD0Cut;
                        std::cout<<" dPhiCut = "<<signalEvents.dPhiCut<<std::endl;
                        signalEvents.dsPlusMCut=dsPlusM_signal;
                        signalEvents.mbcCut=MBC_signal;
                        signalEvents.deltaMCut=DeltaM_signal;
                        signalEvents.diffD0Cut=d0_e_signal-d0_p_signal;
                        signalEvents.dPhiCut=dPhi;
                        ++nRepetitions;
                      }
                      
                      h_MBC_signal->Fill(MBC_signal);
                      if (MBC_signal>lowerSide_xmin && MBC_signal<=lowerSide_xmax) ++signal_beforeSignal;
                      if (MBC_signal>xmin && MBC_signal<xmax) ++signal_atSignal;
                      if (MBC_signal>=upperSide_xmin && MBC_signal<upperSide_xmax) ++signal_afterSignal;
                    }
                    if (DeltaM_sideband)
                    {
                      h_DeltaM_signal->Fill(DeltaM_signal);
                      if (DeltaM_signal>lowerSide_xmin && DeltaM_signal<=lowerSide_xmax) ++signal_beforeSignal;
                      if (DeltaM_signal>xmin && DeltaM_signal<xmax) ++signal_atSignal;
                      if (DeltaM_signal>=upperSide_xmin && DeltaM_signal<upperSide_xmax) ++signal_afterSignal;
                    }
                    if (dPhi_sideband)
                    {
                      h_dPhi_signal->Fill(dPhi);
                      if (dPhi>dPhiCutLess) ++signal_afterSignal;
                      if (dPhi<dPhiCutLess) ++signal_atSignal;
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
  
  std::cout<<"nRepetitions = "<<nRepetitions<<std::endl;
  std::cout<<"Signal processed"<<std::endl;
  
  int nConverEvents_Dsp=converTree_Dsp->GetEntries();
  for (int i=0; i<nConverEvents_Dsp; ++i)
  {
    converTree_Dsp->GetEvent(i);
    
    double phi_e=atan2(py_e_conver_Dsp, px_e_conver_Dsp);
    double phi_p=atan2(py_p_conver_Dsp, px_p_conver_Dsp);
    double dPhi=deltaPhi(phi_e, phi_p);
    
    if (runNumber_conver_Dsp!=converNumber_Dsp.noCut_run || eventNumber_conver_Dsp!=converNumber_Dsp.noCut_event)
    {
      ++converEvents_Dsp.noCut;
      converNumber_Dsp.noCut_run=runNumber_conver_Dsp;
      converNumber_Dsp.noCut_event=eventNumber_conver_Dsp;
    }
    
    if (fabs(d0_e_conver_Dsp)<0.005 && fabs(d0_p_conver_Dsp)<0.005 &&
        fabs(z0_e_conver_Dsp)<0.05 && fabs(z0_p_conver_Dsp)<0.05 &&
        electronEnergyCut(E_e_conver_Dsp, E_p_conver_Dsp))
    {
      if (decayMode_conver_Dsp==decayNumber)
      {
        if (!DsMass_sideband) h_dsPlusM_conver_Dsp->Fill(dsPlusM_conver_Dsp);
        if (DsMass_sideband || dsPlusMCut(dsPlusM_conver_Dsp))
        {
          if (!mBC_sideband) h_MBC_conver_Dsp->Fill(MBC_conver_Dsp);
          if (mBC_sideband || MBCCut(MBC_conver_Dsp))
          {
            if (!DeltaM_sideband) h_DeltaM_conver_Dsp->Fill(DeltaM_conver_Dsp);
            if (DeltaM_sideband || DeltaMCut(DeltaM_conver_Dsp))
            {
              h_diffD0_conver_Dsp->Fill(d0_e_conver_Dsp-d0_p_conver_Dsp);
              h_dPhi_diffD0_conver_Dsp->Fill(dPhi, d0_e_conver_Dsp-d0_p_conver_Dsp);
              if (diffD0_converSide ^ dD0(d0_e_conver_Dsp-d0_p_conver_Dsp))
              {
                //std::cout<<"dD0 = "<<d0_e_conver_Dsp-d0_p_conver_Dsp<<std::endl;
                if (!dPhi_sideband) h_dPhi_conver_Dsp->Fill(dPhi);
                if (dPhi_sideband || dPhiCut(dPhi))
                {
                  h_pi0_conver_Dsp->Fill(pi0Mass_conver_Dsp);
                  if (pi0MassCut(pi0Mass_conver_Dsp))
                  {
                    if (runNumber_conver_Dsp!=converNumber_Dsp.pi0MassCut_run || eventNumber_conver_Dsp!=converNumber_Dsp.pi0MassCut_event)
                    {
                      ++converEvents_Dsp.pi0MassCut;
                      converNumber_Dsp.pi0MassCut_run=runNumber_conver_Dsp;
                      converNumber_Dsp.pi0MassCut_event=eventNumber_conver_Dsp;
                    }
                    if (DsMass_sideband)
                    {
                      h_dsPlusM_conver_Dsp->Fill(dsPlusM_conver_Dsp);
                      if (dsPlusM_conver_Dsp>lowerSide_xmin && dsPlusM_conver_Dsp<=lowerSide_xmax) ++conver_Dsp_beforeSignal;
                      if (dsPlusM_conver_Dsp>xmin && dsPlusM_conver_Dsp<xmax) ++conver_Dsp_atSignal;
                      if (dsPlusM_conver_Dsp>=upperSide_xmin && dsPlusM_conver_Dsp<upperSide_xmax) ++conver_Dsp_afterSignal;
                    }
                    if (mBC_sideband)
                    {
                      h_MBC_conver_Dsp->Fill(MBC_conver_Dsp);
                      if (MBC_conver_Dsp>lowerSide_xmin && MBC_conver_Dsp<=lowerSide_xmax) ++conver_Dsp_beforeSignal;
                      if (MBC_conver_Dsp>xmin && MBC_conver_Dsp<xmax) ++conver_Dsp_atSignal;
                      if (MBC_conver_Dsp>=upperSide_xmin && MBC_conver_Dsp<upperSide_xmax) ++conver_Dsp_afterSignal;
                    }
                    if (DeltaM_sideband)
                    {
                      h_DeltaM_conver_Dsp->Fill(DeltaM_conver_Dsp);
                      if (DeltaM_conver_Dsp>lowerSide_xmin && DeltaM_conver_Dsp<=lowerSide_xmax) ++conver_Dsp_beforeSignal;
                      if (DeltaM_conver_Dsp>xmin && DeltaM_conver_Dsp<xmax) ++conver_Dsp_atSignal;
                      if (DeltaM_conver_Dsp>=upperSide_xmin && DeltaM_conver_Dsp<upperSide_xmax) ++conver_Dsp_afterSignal;
                    }
                    if (dPhi_sideband)
                    {
                      h_dPhi_conver_Dsp->Fill(dPhi);
                      if (dPhi>dPhiCutLess) ++conver_Dsp_afterSignal;
                      if (dPhi<dPhiCutLess) ++conver_Dsp_atSignal;
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
  
  std::cout<<"Ds+ Conversions processed"<<std::endl;
  
  int nConverEvents_Dsm=converTree_Dsm->GetEntries();
  for (int i=0; i<nConverEvents_Dsm; ++i)
  {
    converTree_Dsm->GetEvent(i);
    
    double phi_e=atan2(py_e_conver_Dsm, px_e_conver_Dsm);
    double phi_p=atan2(py_p_conver_Dsm, px_p_conver_Dsm);
    double dPhi=deltaPhi(phi_e, phi_p);
    
    if (runNumber_conver_Dsm!=converNumber_Dsm.noCut_run || eventNumber_conver_Dsm!=converNumber_Dsm.noCut_event)
    {
      ++converEvents_Dsm.noCut;
      converNumber_Dsm.noCut_run=runNumber_conver_Dsm;
      converNumber_Dsm.noCut_event=eventNumber_conver_Dsm;
    }
    
    if (fabs(d0_e_conver_Dsm)<0.005 && fabs(d0_p_conver_Dsm)<0.005 &&
        fabs(z0_e_conver_Dsm)<0.05 && fabs(z0_p_conver_Dsm)<0.05 &&
        electronEnergyCut(E_e_conver_Dsm, E_p_conver_Dsm))
    {
      if (decayMode_conver_Dsm==decayNumber)
      {
        if (!DsMass_sideband) h_dsPlusM_conver_Dsm->Fill(dsPlusM_conver_Dsm);
        if (DsMass_sideband || dsPlusMCut(dsPlusM_conver_Dsm))
        {
          if (!mBC_sideband) h_MBC_conver_Dsm->Fill(MBC_conver_Dsm);
          if (mBC_sideband || MBCCut(MBC_conver_Dsm))
          {
            if (!DeltaM_sideband) h_DeltaM_conver_Dsm->Fill(DeltaM_conver_Dsm);
            if (DeltaM_sideband || DeltaMCut(DeltaM_conver_Dsm))
            {
              h_diffD0_conver_Dsm->Fill(d0_e_conver_Dsm-d0_p_conver_Dsm);
              h_dPhi_diffD0_conver_Dsm->Fill(dPhi, d0_e_conver_Dsm-d0_p_conver_Dsm);
              if (diffD0_converSide ^ dD0(d0_e_conver_Dsm-d0_p_conver_Dsm))
              {
                if (!dPhi_sideband) h_dPhi_conver_Dsm->Fill(dPhi);
                if (dPhi_sideband || dPhiCut(dPhi))
                {
                  h_pi0_conver_Dsm->Fill(pi0Mass_conver_Dsm);
                  if (pi0MassCut(pi0Mass_conver_Dsm))
                  {
                    if (runNumber_conver_Dsm!=converNumber_Dsm.pi0MassCut_run || eventNumber_conver_Dsm!=converNumber_Dsm.pi0MassCut_event)
                    {
                      ++converEvents_Dsm.pi0MassCut;
                      converNumber_Dsm.pi0MassCut_run=runNumber_conver_Dsm;
                      converNumber_Dsm.pi0MassCut_event=eventNumber_conver_Dsm;
                    }
                    if (DsMass_sideband)
                    {
                      h_dsPlusM_conver_Dsm->Fill(dsPlusM_conver_Dsm);
                      if (dsPlusM_conver_Dsm>lowerSide_xmin && dsPlusM_conver_Dsm<=lowerSide_xmax) ++conver_Dsm_beforeSignal;
                      if (dsPlusM_conver_Dsm>xmin && dsPlusM_conver_Dsm<xmax) ++conver_Dsm_atSignal;
                      if (dsPlusM_conver_Dsm>=upperSide_xmin && dsPlusM_conver_Dsm<upperSide_xmax) ++conver_Dsm_afterSignal;
                    }
                    if (mBC_sideband)
                    {
                      h_MBC_conver_Dsm->Fill(MBC_conver_Dsm);
                      if (MBC_conver_Dsm>lowerSide_xmin && MBC_conver_Dsm<=lowerSide_xmax) ++conver_Dsm_beforeSignal;
                      if (MBC_conver_Dsm>xmin && MBC_conver_Dsm<xmax) ++conver_Dsm_atSignal;
                      if (MBC_conver_Dsm>=upperSide_xmin && MBC_conver_Dsm<upperSide_xmax) ++conver_Dsm_afterSignal;
                    }
                    if (DeltaM_sideband)
                    {
                      h_DeltaM_conver_Dsm->Fill(DeltaM_conver_Dsm);
                      if (DeltaM_conver_Dsm>lowerSide_xmin && DeltaM_conver_Dsm<=lowerSide_xmax) ++conver_Dsm_beforeSignal;
                      if (DeltaM_conver_Dsm>xmin && DeltaM_conver_Dsm<xmax) ++conver_Dsm_atSignal;
                      if (DeltaM_conver_Dsm>=upperSide_xmin && DeltaM_conver_Dsm<upperSide_xmax) ++conver_Dsm_afterSignal;
                    }
                    if (dPhi_sideband)
                    {
                      h_dPhi_conver_Dsm->Fill(dPhi);
                      if (dPhi>dPhiCutLess) ++conver_Dsm_afterSignal;
                      if (dPhi<dPhiCutLess) ++conver_Dsm_atSignal;
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
  
  std::cout<<"Ds- Conversions processed"<<std::endl;
  
  int ngenericEvents=genericTree->GetEntries();
  for (int i=0; i<ngenericEvents; ++i)
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
      if (conversionBit_generic==1) ++nConversionEvents_generic_Dsp;
      else if (conversionBit_generic==-1) ++nConversionEvents_generic_Dsm;
    }
    
    if (fabs(d0_e_generic)<0.005 && fabs(d0_p_generic)<0.005 &&
        fabs(z0_e_generic)<0.05 && fabs(z0_p_generic)<0.05 &&
        electronEnergyCut(E_e_generic, E_p_generic) &&
        conversionBit_generic==0)
    {
      if (decayMode_generic==decayNumber)
      {
        if (!DsMass_sideband) h_dsPlusM_generic->Fill(dsPlusM_generic);
        if (DsMass_sideband || dsPlusMCut(dsPlusM_generic))
        {
          if (!mBC_sideband) h_MBC_generic->Fill(MBC_generic);
          if (mBC_sideband || MBCCut(MBC_generic))
          {
            if (!DeltaM_sideband) h_DeltaM_generic->Fill(DeltaM_generic);
            if (DeltaM_sideband || DeltaMCut(DeltaM_generic))
            {
              h_diffD0_generic->Fill(d0_e_generic-d0_p_generic);
              h_dPhi_diffD0_generic->Fill(dPhi, d0_e_generic-d0_p_generic);
              if (diffD0_converSide ^ dD0(d0_e_generic-d0_p_generic))
              {
                if (!dPhi_sideband) h_dPhi_generic->Fill(dPhi);
                if (dPhi_sideband || dPhiCut(dPhi))
                {
                  h_pi0_generic->Fill(pi0Mass_generic);
                  if (pi0MassCut(pi0Mass_generic))
                  {
                    if (runNumber_generic!=genericNumber.pi0MassCut_run || eventNumber_generic!=genericNumber.pi0MassCut_event)
                    {
                      ++genericEvents.pi0MassCut;
                      genericNumber.pi0MassCut_run=runNumber_generic;
                      genericNumber.pi0MassCut_event=eventNumber_generic;
                    }
                    if (DsMass_sideband)
                    {
                      h_dsPlusM_generic->Fill(dsPlusM_generic);
                      if (dsPlusM_generic>lowerSide_xmin && dsPlusM_generic<=lowerSide_xmax) ++generic_beforeSignal;
                      if (dsPlusM_generic>xmin && dsPlusM_generic<xmax) ++generic_atSignal;
                      if (dsPlusM_generic>=upperSide_xmin && dsPlusM_generic<upperSide_xmax) ++generic_afterSignal;
                    }
                    if (mBC_sideband)
                    {
                      h_MBC_generic->Fill(MBC_generic);
                      if (MBC_generic>lowerSide_xmin && MBC_generic<=lowerSide_xmax) ++generic_beforeSignal;
                      if (MBC_generic>xmin && MBC_generic<xmax) ++generic_atSignal;
                      if (MBC_generic>=upperSide_xmin && MBC_generic<upperSide_xmax) ++generic_afterSignal;
                    }
                    if (DeltaM_sideband)
                    {
                      h_DeltaM_generic->Fill(DeltaM_generic);
                      if (DeltaM_generic>lowerSide_xmin && DeltaM_generic<=lowerSide_xmax) ++generic_beforeSignal;
                      if (DeltaM_generic>xmin && DeltaM_generic<xmax) ++generic_atSignal;
                      if (DeltaM_generic>=upperSide_xmin && DeltaM_generic<upperSide_xmax) ++generic_afterSignal;
                    }
                    if (dPhi_sideband)
                    {
                      h_dPhi_generic->Fill(dPhi);
                      if (dPhi>dPhiCutLess) ++generic_afterSignal;
                      if (dPhi<dPhiCutLess) ++generic_atSignal;
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
  
  std::cout<<"Generic processed"<<std::endl;
  
  int ncontinuEvents=continuTree->GetEntries();
  for (int i=0; i<ncontinuEvents; ++i)
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
        if (!DsMass_sideband) h_dsPlusM_continu->Fill(dsPlusM_continu);
        if (DsMass_sideband || dsPlusMCut(dsPlusM_continu))
        {
          if (!mBC_sideband) h_MBC_continu->Fill(MBC_continu);
          if (mBC_sideband || MBCCut(MBC_continu))
          {
            if (!DeltaM_sideband) h_DeltaM_continu->Fill(DeltaM_continu);
            if (DeltaM_sideband || DeltaMCut(DeltaM_continu))
            {
              h_diffD0_continu->Fill(d0_e_continu-d0_p_continu);
              h_dPhi_diffD0_continu->Fill(dPhi, d0_e_continu-d0_p_continu);
              if (diffD0_converSide ^ dD0(d0_e_continu-d0_p_continu))
              {
                if (!dPhi_sideband) h_dPhi_continu->Fill(dPhi);
                if (dPhi_sideband || dPhiCut(dPhi))
                {
                  h_pi0_continu->Fill(pi0Mass_continu);
                  if (pi0MassCut(pi0Mass_continu))
                  {
                    if (runNumber_continu!=continuNumber.pi0MassCut_run || eventNumber_continu!=continuNumber.pi0MassCut_event)
                    {
                      ++continuEvents.pi0MassCut;
                      continuNumber.pi0MassCut_run=runNumber_continu;
                      continuNumber.pi0MassCut_event=eventNumber_continu;
                    }
                    if (DsMass_sideband)
                    {
                      h_dsPlusM_continu->Fill(dsPlusM_continu);
                      if (dsPlusM_continu>lowerSide_xmin && dsPlusM_continu<=lowerSide_xmax) ++continu_beforeSignal;
                      if (dsPlusM_continu>xmin && dsPlusM_continu<xmax) ++continu_atSignal;
                      if (dsPlusM_continu>=upperSide_xmin && dsPlusM_continu<upperSide_xmax) ++continu_afterSignal;
                    }
                    if (mBC_sideband)
                    {
                      h_MBC_continu->Fill(MBC_continu);
                      if (MBC_continu>lowerSide_xmin && MBC_continu<=lowerSide_xmax) ++continu_beforeSignal;
                      if (MBC_continu>xmin && MBC_continu<xmax) ++continu_atSignal;
                      if (MBC_continu>=upperSide_xmin && MBC_continu<upperSide_xmax) ++continu_afterSignal;
                    }
                    if (DeltaM_sideband)
                    {
                      h_DeltaM_continu->Fill(DeltaM_continu);
                      if (DeltaM_continu>lowerSide_xmin && DeltaM_continu<=lowerSide_xmax) ++continu_beforeSignal;
                      if (DeltaM_continu>xmin && DeltaM_continu<xmax) ++continu_atSignal;
                      if (DeltaM_continu>=upperSide_xmin && DeltaM_continu<upperSide_xmax) ++continu_afterSignal;
                    }
                    if (dPhi_sideband)
                    {
                      h_dPhi_continu->Fill(dPhi);
                      if (dPhi>dPhiCutLess) ++continu_afterSignal;
                      if (dPhi<dPhiCutLess) ++continu_atSignal;
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
  
  std::cout<<"Continuum processed"<<std::endl;
  
  int nphysicsEvents=physicsTree->GetEntries();
  for (int i=0; i<nphysicsEvents; ++i)
  {
    physicsTree->GetEvent(i);
    
    double phi_e=atan2(py_e_physics, px_e_physics);
    double phi_p=atan2(py_p_physics, px_p_physics);
    double dPhi=deltaPhi(phi_e, phi_p);
    
    if (runNumber_physics!=physicsNumber.noCut_run || eventNumber_physics!=physicsNumber.noCut_event)
    {      
      ++physicsEvents.noCut;
      physicsNumber.noCut_run=runNumber_physics;
      physicsNumber.noCut_event=eventNumber_physics;
    }
    
    if (fabs(d0_e_physics)<0.005 && fabs(d0_p_physics)<0.005 && 
        fabs(z0_e_physics)<0.05 && fabs(z0_p_physics)<0.05 &&
        electronEnergyCut(E_e_physics, E_p_physics))
    {
      if (decayMode_physics==decayNumber)
      {
        if (!DsMass_sideband) h_dsPlusM_physics->Fill(dsPlusM_physics);
        if (DsMass_sideband || dsPlusMCut(dsPlusM_physics))
        {
          if (!mBC_sideband) h_MBC_physics->Fill(MBC_physics);
          if (mBC_sideband || MBCCut(MBC_physics))
          {
            if (!DeltaM_sideband) h_DeltaM_physics->Fill(DeltaM_physics);
            if (DeltaM_sideband || DeltaMCut(DeltaM_physics))
            {
              h_diffD0_physics->Fill(d0_e_physics-d0_p_physics);
              h_dPhi_diffD0_physics->Fill(dPhi, d0_e_physics-d0_p_physics);
              if (diffD0_converSide ^ dD0(d0_e_physics-d0_p_physics))
              {
                if (!dPhi_sideband) h_dPhi_physics->Fill(dPhi);
                if (dPhi_sideband || dPhiCut(dPhi))
                {
                  h_pi0_physics->Fill(pi0Mass_physics);
                  if (pi0MassCut(pi0Mass_physics))
                  {
                    if (runNumber_physics!=physicsNumber.pi0MassCut_run || eventNumber_physics!=physicsNumber.pi0MassCut_event)
                    {
                      ++physicsEvents.pi0MassCut;
                      physicsNumber.pi0MassCut_run=runNumber_physics;
                      physicsNumber.pi0MassCut_event=eventNumber_physics;
                    }
                    if (DsMass_sideband && !dsPlusMCut(dsPlusM_physics))
                    {
                      h_dsPlusM_physics->Fill(dsPlusM_physics);
                      if (dsPlusM_physics>lowerSide_xmin && dsPlusM_physics<=lowerSide_xmax) ++physics_beforeSignal;
                      if (dsPlusM_physics>xmin && dsPlusM_physics<xmax) ++physics_atSignal;
                      if (dsPlusM_physics>=upperSide_xmin && dsPlusM_physics<upperSide_xmax) ++physics_afterSignal;
                    }
                    if (mBC_sideband && !MBCCut(MBC_physics))
                    {
                      h_MBC_physics->Fill(MBC_physics);
                      if (MBC_physics>lowerSide_xmin && MBC_physics<=lowerSide_xmax) ++physics_beforeSignal;
                      if (MBC_physics>xmin && MBC_physics<xmax) ++physics_atSignal;
                      if (MBC_physics>=upperSide_xmin && MBC_physics<upperSide_xmax) ++physics_afterSignal;
                    }
                    if (DeltaM_sideband && !DeltaMCut(DeltaM_physics))
                    {
                      h_DeltaM_physics->Fill(DeltaM_physics);
                      if (DeltaM_physics>lowerSide_xmin && DeltaM_physics<=lowerSide_xmax) ++physics_beforeSignal;
                      if (DeltaM_physics>xmin && DeltaM_physics<xmax) ++physics_atSignal;
                      if (DeltaM_physics>=upperSide_xmin && DeltaM_physics<upperSide_xmax) ++physics_afterSignal;
                    }
                    if (dPhi_sideband && !dPhiCut(dPhi))
                    {
                      h_dPhi_physics->Fill(dPhi);
                      if (dPhi>dPhiCutLess) ++physics_afterSignal;
                      if (dPhi<dPhiCutLess) ++physics_atSignal;
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
  
  std::cout<<"Data processed"<<std::endl;
  
  std::cout<<"=== Conversions n-tuplizing Efficiency ==="<<std::endl;
  double ntuplizerEff_Dsp=converEvents_Dsp.noCut/nConversionSample_Dsp;
  double ntuplizerEff_Dsp_error=pow(ntuplizerEff_Dsp*(1-ntuplizerEff_Dsp)/nConversionSample_Dsp, 0.5);
  double ntuplizerEff_Dsm=converEvents_Dsm.noCut/nConversionSample_Dsm;
  double ntuplizerEff_Dsm_error=pow(ntuplizerEff_Dsm*(1-ntuplizerEff_Dsm)/nConversionSample_Dsm, 0.5);
  std::cout<<"Number of Dsp conversion events in sample = "<<nConversionSample_Dsp<<std::endl;
  std::cout<<"Number of Dsp conversion events after n-tuplizer = "<<converEvents_Dsp.noCut<<std::endl;
  std::cout<<" Efficiency of n-tuplizer on Dsp = "<<ntuplizerEff_Dsp<<" +- "<<ntuplizerEff_Dsp_error<<std::endl;
  std::cout<<"Number of Dsm conversion events in sample = "<<nConversionSample_Dsm<<std::endl;
  std::cout<<"Number of Dsm conversion events after n-tuplizer = "<<converEvents_Dsm.noCut<<std::endl;
  std::cout<<" Efficiency of n-tuplizer on Dsm = "<<ntuplizerEff_Dsm<<" +- "<<ntuplizerEff_Dsm_error<<std::endl;
  std::cout<<"=== Conversion Scales ==="<<std::endl;
  double expectedConversionEvents_Dsp=luminosity*prodCrossSection_DsDss*branchingFr_Dsstgamma*ntuplizerEff_Dsp/2;
  double expectedConversionEvents_Dsp_error=expectedConversionEvents_Dsp*pow(pow(luminosity_error/luminosity, 2)
                                                                            +pow(prodCrossSection_DsDss_error/prodCrossSection_DsDss, 2)
                                                                            +pow(branchingFr_Dsstgamma_error/branchingFr_Dsstgamma, 2)
                                                                            +pow(ntuplizerEff_Dsp_error/ntuplizerEff_Dsp, 2), 0.5);
  double expectedConversionEvents_Dsm=luminosity*prodCrossSection_DsDss*branchingFr_Dsstgamma*ntuplizerEff_Dsm/2;
  double expectedConversionEvents_Dsm_error=expectedConversionEvents_Dsm*pow(pow(luminosity_error/luminosity, 2)
                                                                            +pow(prodCrossSection_DsDss_error/prodCrossSection_DsDss, 2)
                                                                            +pow(branchingFr_Dsstgamma_error/branchingFr_Dsstgamma, 2)
                                                                            +pow(ntuplizerEff_Dsm_error/ntuplizerEff_Dsm, 2), 0.5);
  std::cout<<"Number of conversion Dsp events expected to be produced in "<<luminosity<<" /pb and survive n-tuplizing =  "<<expectedConversionEvents_Dsp<<"+-"<<expectedConversionEvents_Dsp_error<<std::endl;
  std::cout<<"Number of conversion Dsm events expected to be produced in "<<luminosity<<" /pb and survive n-tuplizing =  "<<expectedConversionEvents_Dsm<<"+-"<<expectedConversionEvents_Dsm_error<<std::endl;
  double seenConversionEvents_Dsp=nConversionEvents_generic_Dsp*genericScale;
  double seenConversionEvents_Dsp_error=nConversionEvents_generic_Dsp*genericScale_error;
  double seenConversionEvents_Dsm=nConversionEvents_generic_Dsm*genericScale;
  double seenConversionEvents_Dsm_error=nConversionEvents_generic_Dsm*genericScale_error;
  std::cout<<"Number of conversion events seen in Dsp = "<<seenConversionEvents_Dsp<<" +- "<<seenConversionEvents_Dsp_error<<std::endl;
  std::cout<<"Number of conversion events seen in Dsm = "<<seenConversionEvents_Dsm<<" +- "<<seenConversionEvents_Dsm_error<<std::endl;
  if (fabs(expectedConversionEvents_Dsp-seenConversionEvents_Dsp)>(expectedConversionEvents_Dsp_error+seenConversionEvents_Dsp_error)) std::cout<<"WARNING!! Dsp conversion numbers are off."<<std::endl;
  if (fabs(expectedConversionEvents_Dsm-seenConversionEvents_Dsm)>(expectedConversionEvents_Dsm_error+seenConversionEvents_Dsm_error)) std::cout<<"WARNING!! Dsm conversion numbers are off."<<std::endl;
  double converScale_Dsp=seenConversionEvents_Dsp/converEvents_Dsp.noCut;
  double converScale_Dsp_error=converScale_Dsp*pow(pow(seenConversionEvents_Dsp_error/seenConversionEvents_Dsp, 2)+1/converEvents_Dsp.noCut, 0.5);
  double converScale_Dsm=seenConversionEvents_Dsm/converEvents_Dsm.noCut;
  double converScale_Dsm_error=converScale_Dsm*pow(pow(seenConversionEvents_Dsm_error/seenConversionEvents_Dsm, 2)+1/converEvents_Dsm.noCut, 0.5);
  std::cout<<"Dsp Conversion scale = seenConversionEvents_Dsp/converEvents_Dsp.noCut = "<<converScale_Dsp<<"+-"<<converScale_Dsp_error<<std::endl;
  std::cout<<"Dsm Conversion scale = seenConversionEvents_Dsm/converEvents_Dsm.noCut = "<<converScale_Dsm<<"+-"<<converScale_Dsm_error<<std::endl;
  std::cout<<"=== Signal Scales ==="<<std::endl;
  std::cout<<"nSignalSample = "<<nSignalSample<<std::endl;
  std::cout<<"Produced signal events = "<<(luminosity*prodCrossSection_DsDss*branchingFr_signal*branchingFr_mode)<<std::endl;
  double signalScale=(luminosity*prodCrossSection_DsDss*branchingFr_signal*branchingFr_mode)/nSignalSample;
  double signalScale_error=signalScale*pow(pow(luminosity_error/luminosity, 2)+pow(prodCrossSection_DsDss_error/prodCrossSection_DsDss, 2), 0.5);
  std::cout<<"Signal Scale = "<<signalScale<<std::endl;
  std::cout<<"=== Signal Efficiency ==="<<std::endl;
  double signalEfficiency=signal_atSignal/nSignalSample;
  double signalEfficiency_error=signalEfficiency/sqrt(signal_atSignal);
  std::cout<<"Signal Efficiency = "<<signal_atSignal<<"/"<<nSignalSample<<"="<<signal_atSignal/nSignalSample<<" +- "<<signalEfficiency_error<<std::endl;
  std::cout<<"==="<<std::endl;
  
  // Scale plots appropriately
  h_dsPlusM_signal->Scale(signalScale);
  h_dsPlusM_conver_Dsp->Scale(converScale_Dsp);
  h_dsPlusM_conver_Dsm->Scale(converScale_Dsm);
  h_dsPlusM_generic->Scale(genericScale);
  h_dsPlusM_continu->Scale(continuScale);
  h_MBC_signal->Scale(signalScale);
  h_MBC_conver_Dsp->Scale(converScale_Dsp);
  h_MBC_conver_Dsm->Scale(converScale_Dsm);
  h_MBC_generic->Scale(genericScale);
  h_MBC_continu->Scale(continuScale);
  h_DeltaM_signal->Scale(signalScale);
  h_DeltaM_conver_Dsp->Scale(converScale_Dsp);
  h_DeltaM_conver_Dsm->Scale(converScale_Dsm);
  h_DeltaM_generic->Scale(genericScale);
  h_DeltaM_continu->Scale(continuScale);
  h_diffD0_signal->Scale(signalScale);
  h_diffD0_conver_Dsp->Scale(converScale_Dsp);
  h_diffD0_conver_Dsm->Scale(converScale_Dsm);
  h_diffD0_generic->Scale(genericScale);
  h_diffD0_continu->Scale(continuScale);
  h_dPhi_signal->Scale(signalScale);
  h_dPhi_conver_Dsp->Scale(converScale_Dsp);
  h_dPhi_conver_Dsm->Scale(converScale_Dsm);
  h_dPhi_generic->Scale(genericScale);
  h_dPhi_continu->Scale(continuScale);
  h_pi0_signal->Scale(signalScale);
  h_pi0_conver_Dsp->Scale(converScale_Dsp);
  h_pi0_conver_Dsm->Scale(converScale_Dsm);
  h_pi0_generic->Scale(genericScale);
  h_pi0_continu->Scale(continuScale);
  
  // And scale numbers too
  signal_beforeSignal*=signalScale; signal_afterSignal*=signalScale; signal_atSignal*=signalScale;
  conver_Dsp_beforeSignal*=converScale_Dsp; conver_Dsp_afterSignal*=converScale_Dsp; conver_Dsp_atSignal*=converScale_Dsp;
  conver_Dsm_beforeSignal*=converScale_Dsm; conver_Dsm_afterSignal*=converScale_Dsm; conver_Dsm_atSignal*=converScale_Dsm;
  generic_beforeSignal*=genericScale; generic_afterSignal*=genericScale; generic_atSignal*=genericScale;
  continu_beforeSignal*=continuScale; continu_afterSignal*=continuScale; continu_atSignal*=continuScale;
  
  // Now calculate errors
  double signal_beforeSignal_error;
  if (signal_beforeSignal>0) signal_beforeSignal_error=signal_beforeSignal*pow(pow(signalScale_error/signalScale, 2)+signalScale/signal_beforeSignal, 0.5);
  else signal_beforeSignal_error=0;
  
  double signal_afterSignal_error;
  if (signal_afterSignal>0) signal_afterSignal_error=signal_afterSignal*pow(pow(signalScale_error/signalScale, 2)+signalScale/signal_afterSignal, 0.5);
  else signal_afterSignal_error=0;
  
  double signal_atSignal_error;
  if (signal_atSignal>0) signal_atSignal_error=signal_atSignal*pow(pow(signalScale_error/signalScale, 2)+signalScale/signal_atSignal, 0.5);
  else signal_atSignal_error=0;
  
  double conver_Dsp_beforeSignal_error;
  if (conver_Dsp_beforeSignal>0) conver_Dsp_beforeSignal_error=conver_Dsp_beforeSignal*pow(pow(converScale_Dsp_error/converScale_Dsp, 2)+converScale_Dsp/conver_Dsp_beforeSignal, 0.5);
  else conver_Dsp_beforeSignal_error=0;
  
  double conver_Dsp_afterSignal_error;
  if (conver_Dsp_afterSignal>0) conver_Dsp_afterSignal_error=conver_Dsp_afterSignal*pow(pow(converScale_Dsp_error/converScale_Dsp, 2)+converScale_Dsp/conver_Dsp_afterSignal, 0.5);
  else conver_Dsp_afterSignal_error=0;
  
  double conver_Dsp_atSignal_error;
  if (conver_Dsp_atSignal>0) conver_Dsp_atSignal_error=conver_Dsp_atSignal*pow(pow(converScale_Dsp_error/converScale_Dsp, 2)+converScale_Dsp/conver_Dsp_atSignal, 0.5);
  else conver_Dsp_atSignal_error=0;
  
  double conver_Dsm_beforeSignal_error;
  if (conver_Dsm_beforeSignal>0) conver_Dsm_beforeSignal_error=conver_Dsm_beforeSignal*pow(pow(converScale_Dsm_error/converScale_Dsm, 2)+converScale_Dsm/conver_Dsm_beforeSignal, 0.5);
  else conver_Dsm_beforeSignal_error=0;
  
  double conver_Dsm_afterSignal_error;
  if (conver_Dsm_afterSignal>0) conver_Dsm_afterSignal_error=conver_Dsm_afterSignal*pow(pow(converScale_Dsm_error/converScale_Dsm, 2)+converScale_Dsm/conver_Dsm_afterSignal, 0.5);
  else conver_Dsm_afterSignal_error=0;
  
  double conver_Dsm_atSignal_error;
  if (conver_Dsm_atSignal>0) conver_Dsm_atSignal_error=conver_Dsm_atSignal*pow(pow(converScale_Dsm_error/converScale_Dsm, 2)+converScale_Dsm/conver_Dsm_atSignal, 0.5);
  else conver_Dsm_atSignal_error=0;
  
  double generic_beforeSignal_error;
  if (generic_beforeSignal>0) generic_beforeSignal_error=generic_beforeSignal*pow(pow(genericScale_error/genericScale, 2)+genericScale/generic_beforeSignal, 0.5);
  else generic_beforeSignal_error=0;
  
  double generic_afterSignal_error;
  if (generic_afterSignal>0) generic_afterSignal_error=generic_afterSignal*pow(pow(genericScale_error/genericScale, 2)+genericScale/generic_afterSignal, 0.5);
  else generic_afterSignal_error=0;
  
  double generic_atSignal_error;
  if (generic_atSignal>0) generic_atSignal_error=generic_atSignal*pow(pow(genericScale_error/genericScale, 2)+genericScale/generic_atSignal, 0.5);
  else generic_atSignal_error=0;
  
  double continu_beforeSignal_error;
  if (continu_beforeSignal>0) continu_beforeSignal_error=continu_beforeSignal*pow(pow(continuScale_error/continuScale, 2)+continuScale/continu_beforeSignal, 0.5);
  else continu_beforeSignal_error=0;
  
  double continu_afterSignal_error;
  if (continu_afterSignal>0) continu_afterSignal_error=continu_afterSignal*pow(pow(continuScale_error/continuScale, 2)+continuScale/continu_afterSignal, 0.5);
  else continu_afterSignal_error=0;
  
  double continu_atSignal_error;
  if (continu_atSignal>0) continu_atSignal_error=continu_atSignal*pow(pow(continuScale_error/continuScale, 2)+continuScale/continu_atSignal, 0.5);
  else continu_atSignal_error=0;
  
  double mc_beforeSignal=conver_Dsp_beforeSignal+conver_Dsm_beforeSignal+generic_beforeSignal+continu_beforeSignal;
  double mc_beforeSignal_error=pow(pow(conver_Dsp_beforeSignal_error, 2)+pow(conver_Dsm_beforeSignal_error, 2)+pow(generic_beforeSignal_error, 2)+pow(continu_beforeSignal_error, 2), 0.5);
  double mc_afterSignal=conver_Dsp_afterSignal+conver_Dsm_afterSignal+generic_afterSignal+continu_afterSignal;
  double mc_afterSignal_error=pow(pow(conver_Dsp_afterSignal_error, 2)+pow(conver_Dsm_afterSignal_error, 2)+pow(generic_afterSignal_error, 2)+pow(continu_afterSignal_error, 2), 0.5);
  double mc_atSignal=conver_Dsp_atSignal+conver_Dsm_atSignal+generic_atSignal+continu_atSignal;
  double mc_atSignal_error=pow(pow(conver_Dsp_atSignal_error, 2)+pow(conver_Dsm_atSignal_error, 2)+pow(generic_atSignal_error, 2)+pow(continu_atSignal_error, 2), 0.5);
  double physics_beforeSignal_error=pow(physics_beforeSignal, 0.5);
  double physics_afterSignal_error=pow(physics_afterSignal, 0.5);
  double physics_atSignal_error=pow(physics_atSignal, 0.5);
  double dataMC_beforeSignal=physics_beforeSignal/mc_beforeSignal;
  double dataMC_beforeSignal_error=dataMC_beforeSignal*pow(1/physics_beforeSignal+pow(mc_beforeSignal_error/mc_beforeSignal, 2), 0.5);
  double dataMC_afterSignal=physics_afterSignal/mc_afterSignal;
  double dataMC_afterSignal_error=dataMC_afterSignal*pow(1/physics_afterSignal+pow(mc_afterSignal_error/mc_afterSignal, 2), 0.5);
  double dataMC_atSignal=physics_atSignal/mc_atSignal;
  double dataMC_atSignal_error=dataMC_atSignal*pow(1/physics_atSignal+pow(mc_atSignal_error/mc_atSignal, 2), 0.5);
  
  std::cout<<"Sideband summary for "<<decay<<std::endl;
  
  gROOT->SetStyle("Plain");
  double ymin, ymax;
  
  if (mBC_sideband)
  {
    std::cout<<"=== mBC sideband ==="<<std::endl;
    std::cout<<"Signal MC: ";
    std::cout<<" before signal = "<<signal_beforeSignal<<"+-"<<signal_beforeSignal_error;
    std::cout<<", after signal = "<<signal_afterSignal<<"+-"<<signal_afterSignal_error;
    std::cout<<", at signal = "<<signal_atSignal<<"+-"<<signal_atSignal_error<<std::endl;
    std::cout<<"Conversion Dsp: ";
    std::cout<<" before signal = "<<conver_Dsp_beforeSignal<<"+-"<<conver_Dsp_beforeSignal_error;
    std::cout<<", after signal = "<<conver_Dsp_afterSignal<<"+-"<<conver_Dsp_afterSignal_error;
    std::cout<<", at signal = "<<conver_Dsp_atSignal<<"+-"<<conver_Dsp_atSignal_error<<std::endl;
    std::cout<<"Conversion Dsm: ";
    std::cout<<" before signal = "<<conver_Dsm_beforeSignal<<"+-"<<conver_Dsm_beforeSignal_error;
    std::cout<<", after signal = "<<conver_Dsm_afterSignal<<"+-"<<conver_Dsm_afterSignal_error;
    std::cout<<", at signal = "<<conver_Dsm_atSignal<<"+-"<<conver_Dsm_atSignal_error<<std::endl;
    std::cout<<"Generic MC: ";
    std::cout<<" before signal = "<<generic_beforeSignal<<"+-"<<generic_beforeSignal_error;
    std::cout<<", after signal = "<<generic_afterSignal<<"+-"<<generic_afterSignal_error;
    std::cout<<", at signal = "<<generic_atSignal<<"+-"<<generic_atSignal_error<<std::endl;
    std::cout<<"Continuum MC: ";
    std::cout<<" before signal = "<<continu_beforeSignal<<"+-"<<continu_beforeSignal_error;
    std::cout<<", after signal = "<<continu_afterSignal<<"+-"<<continu_afterSignal_error;
    std::cout<<", at signal = "<<continu_atSignal<<"+-"<<continu_atSignal_error<<std::endl;
    std::cout<<"Sum MC: ";
    std::cout<<" before signal = "<<mc_beforeSignal<<"+-"<<mc_beforeSignal_error;
    std::cout<<", after signal = "<<mc_afterSignal<<"+-"<<mc_afterSignal_error;
    std::cout<<", at signal = "<<mc_atSignal<<"+-"<<mc_atSignal_error<<std::endl;
    std::cout<<"Data (electron-fitted): "<<std::endl;
    std::cout<<" before signal = "<<physics_beforeSignal<<"+-"<<physics_beforeSignal_error<<" Data/MC = "<<dataMC_beforeSignal<<"+-"<<dataMC_beforeSignal_error<<std::endl;
    std::cout<<" after signal = "<<physics_afterSignal<<"+-"<<physics_afterSignal_error<<" Data/MC = "<<dataMC_afterSignal<<"+-"<<dataMC_afterSignal_error<<std::endl;
    std::cout<<" at signal = "<<physics_atSignal<<"+-"<<physics_atSignal_error<<" Data/MC = "<<dataMC_atSignal<<"+-"<<dataMC_atSignal_error<<std::endl;
    std::cout<<"===================="<<std::endl;
    std::cout<<" TeX output "<<std::endl;
    std::cout<<setprecision(2);
    std::cout<<std::fixed;
    std::cout<<" & "<<conver_Dsp_beforeSignal+conver_Dsm_beforeSignal<<" $\\pm$ "<<pow(pow(conver_Dsp_beforeSignal_error, 2)+pow(conver_Dsm_beforeSignal_error, 2), 0.5);
    std::cout<<" & "<<generic_beforeSignal<<" $\\pm$ "<<generic_beforeSignal_error;
    std::cout<<" & "<<continu_beforeSignal<<" $\\pm$ "<<continu_beforeSignal_error;
    std::cout<<" & "<<mc_beforeSignal<<" $\\pm$ "<<mc_beforeSignal_error;
    std::cout<<" & "<<physics_beforeSignal<<" $\\pm$ "<<physics_beforeSignal_error;
    std::cout<<" & "<<dataMC_beforeSignal<<" $\\pm$ "<<dataMC_beforeSignal_error;
    std::cout<<" & "<<conver_Dsp_afterSignal+conver_Dsm_afterSignal<<" $\\pm$ "<<pow(pow(conver_Dsp_afterSignal_error, 2)+pow(conver_Dsm_afterSignal_error, 2), 0.5);    
    std::cout<<" & "<<generic_afterSignal<<" $\\pm$ "<<generic_afterSignal_error;
    std::cout<<" & "<<continu_afterSignal<<" $\\pm$ "<<continu_afterSignal_error;
    std::cout<<" & "<<mc_afterSignal<<" $\\pm$ "<<mc_afterSignal_error;
    std::cout<<" & "<<physics_afterSignal<<" $\\pm$ "<<physics_afterSignal_error;
    std::cout<<" & "<<dataMC_afterSignal<<" $\\pm$ "<<dataMC_afterSignal_error<<" \\\\  \\hline"<<std::endl;
    std::cout<<setprecision(3);
    std::cout<<fixed;
    
    TLine *line;
    /*
    TCanvas *MBC = new TCanvas("MBC", "", 400, 600);
    MBC->Divide(2,3);
    MBC->cd(1);
    h_MBC_signal->Draw();
    ymax=(h_MBC_signal->GetMaximum())*0.75;
    line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
    line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
    line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
    MBC->cd(2);
    h_MBC_conver_Dsp->Draw();
    ymax=(h_MBC_conver_Dsp->GetMaximum())*0.75;
    line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
    line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
    line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
    MBC->cd(3);
    h_MBC_conver_Dsm->Draw();
    ymax=(h_MBC_conver_Dsm->GetMaximum())*0.75;
    line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
    line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
    line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
    MBC->cd(4);
    h_MBC_generic->Draw();
    ymax=(h_MBC_generic->GetMaximum())*0.75;
    line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
    line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
    line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
    MBC->cd(5);
    h_MBC_continu->Draw();
    ymax=(h_MBC_continu->GetMaximum())*0.75;
    line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
    line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
    line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
    MBC->cd(6);
    h_MBC_physics->Draw();
    ymax=(h_MBC_physics->GetMaximum())*0.75;
    line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
    line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
    line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
    */
    
    h_MBC_conver_Dsp->SetFillColor(kRed);
    h_MBC_conver_Dsm->SetFillColor(kRed);
    h_MBC_generic->SetFillColor(kGreen);
    h_MBC_continu->SetFillColor(kBlue);
    h_MBC_signal->SetFillColor(kCyan);
    h_MBC_signal->SetLineColor(kBlack);
    
    /*
    if (decay=="pietaprime") // artificially add a data point
    {
      h_MBC_physics->Fill(2.08);
    }
    */
    
    /*
    // Fit the MC and data histograms and calculate expected background in signal region.
    // Fit MC
    TF1 *f_Sum_background=new TF1("f_Sum_background", mbc_MCFit, lowerSide_xmin, upperSide_xmax, 1);
    f_Sum_background->SetLineColor(kBlack);
    h_MBC_physics->Fit("f_Sum_background", "RLLEFM");
    double expBg_Sides=h_MBC_physics->Integral(h_MBC_physics->FindBin(lowerSide_xmin), h_MBC_physics->FindBin(lowerSide_xmax))
                      +h_MBC_physics->Integral(h_MBC_physics->FindBin(upperSide_xmin), h_MBC_physics->FindBin(upperSide_xmax));
    double expBg_MC=f_Sum_background->Integral(xmin, xmax)*500;
    double expBg_MC_error=expBg_MC/pow(expBg_Sides, 0.5);
    std::cout<<"Expected background from MC = "<<expBg_MC<<" +- "<<expBg_MC_error<<std::endl;
    // Fit Data
    TF1* f_physics=new TF1("f_physics", mbc_DataFit, lowerSide_xmin, upperSide_xmax, 1);
    f_physics->SetLineColor(kMagenta);
    h_MBC_physics->Fit(f_physics, "RLLEFM+");
    double expBg_data=f_physics->Integral(xmin, xmax)*500;
    double expBg_data_error=expBg_data/pow(expBg_Sides, 0.5);
    std::cout<<"Expected background from data = "<<expBg_data<<" +- "<<expBg_data_error<<std::endl;
    */
    
    TCanvas *c_MBC_Stacked = new TCanvas("c_MBC_Stacked");
    THStack *s_MBC_sideband=new THStack("s_MBC_sideband", ("m_{BC} Sidebands in "+decay).c_str());
    s_MBC_sideband->Add(h_MBC_continu, "hist");
    s_MBC_sideband->Add(h_MBC_generic, "hist");
    s_MBC_sideband->Add(h_MBC_conver_Dsp, "hist");
    s_MBC_sideband->Add(h_MBC_conver_Dsm, "hist");
    s_MBC_sideband->Add(h_MBC_signal, "hist");
    double stack_ymax=h_MBC_signal->GetMaximum();
    double physics_ymax=h_MBC_physics->GetMaximum();
    if (physics_ymax>stack_ymax) 
    {
      s_MBC_sideband->SetMaximum(physics_ymax*1.5);
      ymax=physics_ymax;
    }
    else
    {
      s_MBC_sideband->SetMaximum(stack_ymax*1.5);
      ymax=stack_ymax;
    }
    s_MBC_sideband->Draw();
    h_MBC_physics->Draw("SAME");
    s_MBC_sideband->GetXaxis()->SetTitle("m_{BC} (GeV)");
    s_MBC_sideband->GetYaxis()->SetTitle("Number of Events / 2 MeV");
    line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
    line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
    line=new TLine(xmin, ymax, xmin, 0); line->Draw();
    line=new TLine(xmax, ymax, xmax, 0); line->Draw();
    line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
    line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
    TLegend *legendMBC=new TLegend(mbc_x1, mbc_y1, mbc_x2, mbc_y2);
    std::string signal_string="Signal MC: "; signal_string+=itoa(signalScale*h_MBC_signal->GetEntries()); signal_string+=" Events";
    std::string continu_string="Continuum MC: "; continu_string+=itoa(continuScale*h_MBC_continu->GetEntries()); continu_string+=" Events";
    std::string generic_string="Generic MC: "; generic_string+=itoa(genericScale*h_MBC_generic->GetEntries()); generic_string+=" Events";
    std::string conver_string="Conversion MC: "; conver_string+=itoa(converScale_Dsp*h_MBC_conver_Dsp->GetEntries()+converScale_Dsm*h_MBC_conver_Dsm->GetEntries()); conver_string+=" Events";
    std::string physics_string="Data: "; physics_string+=itoa(h_MBC_physics->GetEntries()); physics_string+=" Events";
    legendMBC->AddEntry(h_MBC_signal, signal_string.c_str());
    legendMBC->AddEntry(h_MBC_continu, continu_string.c_str());
    legendMBC->AddEntry(h_MBC_generic, generic_string.c_str());
    legendMBC->AddEntry(h_MBC_conver_Dsp, continu_string.c_str());
    legendMBC->AddEntry(h_MBC_physics, physics_string.c_str());
    legendMBC->SetFillColor(kWhite);
    legendMBC->Draw();
    /*
    std::string filename_eps=decay+"_Sideband_MBC_002.eps";
    c_MBC_Stacked->SaveAs(filename_eps.c_str());
    */
    
    // Calculate Signal Efficiency
    // ---------------------------
    std::string title="m_{BC} Distribution in D_{s}^{*#pm} #rightarrow D_{s}^{#pm} e^{+} e^{-}, D_s^{#pm} #rightarrow ";
    title+=decay_tex;
    TCanvas *c_MBC_Signal = new TCanvas("c_MBC_Signal", "c_MBC_Signal");
    h_MBC_signal->SetTitle(title.c_str());
    h_MBC_signal->Scale(1./(nSignalEvents*signalScale));
    h_MBC_signal->GetYaxis()->SetTitle("Efficiency / 2 MeV");
    h_MBC_signal->GetXaxis()->SetTitle("m_{BC} (GeV)");
    h_MBC_signal->GetXaxis()->CenterTitle();
    h_MBC_signal->GetYaxis()->CenterTitle();
    h_MBC_signal->GetYaxis()->SetTitleOffset(1.3);
    h_MBC_signal->Draw();
    ymax=h_MBC_signal->GetMaximum()*0.9;
    line=new TLine(xmin, ymax, xmin, 0); line->Draw();
    line=new TLine(xmax, ymax, xmax, 0); line->Draw();
    std::string filename_eps=decay+"_SignalEfficiency_MBC.eps";
    c_MBC_Signal->SaveAs(filename_eps.c_str());
    
    // Save histograms in a ROOT file
    /*
    std::string filename_root=decay+"_Sideband_MBC_002.root";
    TFile *tfile=new TFile(filename_root.c_str(), "RECREATE");
    h_MBC_signal->Write();
    h_MBC_conver_Dsp->Write();
    h_MBC_conver_Dsm->Write();
    h_MBC_generic->Write();
    h_MBC_continu->Write();
    h_MBC_physics->Write();
    tfile->Close();
    */
  }
  
  if (DeltaM_sideband)
  {
    std::cout<<"=== DeltaM sideband ==="<<std::endl;
    std::cout<<"Signal MC: ";
    std::cout<<" before signal = "<<signal_beforeSignal<<"+-"<<signal_beforeSignal_error;
    std::cout<<", after signal = "<<signal_afterSignal<<"+-"<<signal_afterSignal_error;
    std::cout<<", at signal = "<<signal_atSignal<<"+-"<<signal_atSignal_error<<std::endl;
    std::cout<<"Conversion Dsp: ";
    std::cout<<" before signal = "<<conver_Dsp_beforeSignal<<"+-"<<conver_Dsp_beforeSignal_error;
    std::cout<<", after signal = "<<conver_Dsp_afterSignal<<"+-"<<conver_Dsp_afterSignal_error;
    std::cout<<", at signal = "<<conver_Dsp_atSignal<<"+-"<<conver_Dsp_atSignal_error<<std::endl;
    std::cout<<"Conversion Dsm: ";
    std::cout<<" before signal = "<<conver_Dsm_beforeSignal<<"+-"<<conver_Dsm_beforeSignal_error;
    std::cout<<", after signal = "<<conver_Dsm_afterSignal<<"+-"<<conver_Dsm_afterSignal_error;
    std::cout<<", at signal = "<<conver_Dsm_atSignal<<"+-"<<conver_Dsm_atSignal_error<<std::endl;
    std::cout<<"Generic MC: ";
    std::cout<<" before signal = "<<generic_beforeSignal<<"+-"<<generic_beforeSignal_error;
    std::cout<<", after signal = "<<generic_afterSignal<<"+-"<<generic_afterSignal_error;
    std::cout<<", at signal = "<<generic_atSignal<<"+-"<<generic_atSignal_error<<std::endl;
    std::cout<<"Continuum MC: ";
    std::cout<<" before signal = "<<continu_beforeSignal<<"+-"<<continu_beforeSignal_error;
    std::cout<<", after signal = "<<continu_afterSignal<<"+-"<<continu_afterSignal_error;
    std::cout<<", at signal = "<<continu_atSignal<<"+-"<<continu_atSignal_error<<std::endl;
    std::cout<<"Sum MC: ";
    std::cout<<" before signal = "<<mc_beforeSignal<<"+-"<<mc_beforeSignal_error;
    std::cout<<", after signal = "<<mc_afterSignal<<"+-"<<mc_afterSignal_error;
    std::cout<<", at signal = "<<mc_atSignal<<"+-"<<mc_atSignal_error<<std::endl;
    std::cout<<"Data (electron-fitted): "<<std::endl;
    std::cout<<" before signal = "<<physics_beforeSignal<<"+-"<<physics_beforeSignal_error<<" Data/MC = "<<dataMC_beforeSignal<<"+-"<<dataMC_beforeSignal_error<<std::endl;
    std::cout<<" after signal = "<<physics_afterSignal<<"+-"<<physics_afterSignal_error<<" Data/MC = "<<dataMC_afterSignal<<"+-"<<dataMC_afterSignal_error<<std::endl;
    std::cout<<" at signal = "<<physics_atSignal<<"+-"<<physics_atSignal_error<<" Data/MC = "<<dataMC_atSignal<<"+-"<<dataMC_atSignal_error<<std::endl;
    std::cout<<"===================="<<std::endl;
    std::cout<<" TeX output "<<std::endl;
    std::cout<<setprecision(2);
    std::cout<<std::fixed;
    std::cout<<" & "<<conver_Dsp_beforeSignal+conver_Dsm_beforeSignal<<" $\\pm$ "<<pow(pow(conver_Dsp_beforeSignal_error, 2)+pow(conver_Dsm_beforeSignal_error, 2), 0.5);
    std::cout<<" & "<<generic_beforeSignal<<" $\\pm$ "<<generic_beforeSignal_error;
    std::cout<<" & "<<continu_beforeSignal<<" $\\pm$ "<<continu_beforeSignal_error;
    std::cout<<" & "<<mc_beforeSignal<<" $\\pm$ "<<mc_beforeSignal_error;
    std::cout<<" & "<<physics_beforeSignal<<" $\\pm$ "<<physics_beforeSignal_error;
    std::cout<<" & "<<dataMC_beforeSignal<<" $\\pm$ "<<dataMC_beforeSignal_error;
    std::cout<<" & "<<conver_Dsp_afterSignal+conver_Dsm_afterSignal<<" $\\pm$ "<<pow(pow(conver_Dsp_afterSignal_error, 2)+pow(conver_Dsm_afterSignal_error, 2), 0.5);    
    std::cout<<" & "<<generic_afterSignal<<" $\\pm$ "<<generic_afterSignal_error;
    std::cout<<" & "<<continu_afterSignal<<" $\\pm$ "<<continu_afterSignal_error;
    std::cout<<" & "<<mc_afterSignal<<" $\\pm$ "<<mc_afterSignal_error;
    std::cout<<" & "<<physics_afterSignal<<" $\\pm$ "<<physics_afterSignal_error;
    std::cout<<" & "<<dataMC_afterSignal<<" $\\pm$ "<<dataMC_afterSignal_error<<" \\\\  \\hline"<<std::endl;
    std::cout<<setprecision(3);
    std::cout<<std::fixed;
    
    TLine *line;
    /*
    TCanvas *DeltaM = new TCanvas("DeltaM", "", 400, 600);
    DeltaM->Divide(2,3);
    DeltaM->cd(1);
    h_DeltaM_signal->Draw();
    ymax=(h_DeltaM_signal->GetMaximum())*0.75;
    line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
    line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
    line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
    DeltaM->cd(2);
    h_DeltaM_conver_Dsp->Draw();
    ymax=(h_DeltaM_conver_Dsp->GetMaximum())*0.75;
    line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
    line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
    line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
    DeltaM->cd(3);
    h_DeltaM_conver_Dsm->Draw();
    ymax=(h_DeltaM_conver_Dsm->GetMaximum())*0.75;
    line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
    line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
    line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
    DeltaM->cd(4);
    h_DeltaM_generic->Draw();
    ymax=(h_DeltaM_generic->GetMaximum())*0.75;
    line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
    line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
    line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
    DeltaM->cd(5);
    h_DeltaM_continu->Draw();
    ymax=(h_DeltaM_continu->GetMaximum())*0.75;
    line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
    line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
    line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
    DeltaM->cd(6);
    h_DeltaM_physics->Draw();
    ymax=(h_DeltaM_physics->GetMaximum())*0.75;
    line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
    line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
    line = new TLine(xmin,ymax,xmin,0); line->Draw();
    line = new TLine(xmax,ymax,xmax,0); line->Draw();
    line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
    line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
    */
    
    h_DeltaM_conver_Dsp->SetFillColor(kRed); 
    h_DeltaM_conver_Dsm->SetFillColor(kRed); 
    h_DeltaM_generic->SetFillColor(kGreen); 
    h_DeltaM_continu->SetFillColor(kBlue); 
    h_DeltaM_signal->SetFillColor(kCyan); 
    h_DeltaM_signal->SetLineColor(kBlack);
    
    /*
    if (decay=="pietaprime") // artificially add a data point
    {
      h_DeltaM_physics->Fill(0.2);
    }
    */
    
    
    // Fit the MC and data histograms and calculate expected background in signal region.
    // Fit MC
    TF1 *f_Sum_background=new TF1("f_Sum_background", deltaM_MCFit, 0.1, 0.25, 1);
    f_Sum_background->SetLineColor(kBlack);
    h_DeltaM_physics->Fit("f_Sum_background", "RLLEFM");
    double expBg_Sides=h_DeltaM_physics->Integral(h_DeltaM_physics->FindBin(lowerSide_xmin), h_DeltaM_physics->FindBin(lowerSide_xmax))
                      +h_DeltaM_physics->Integral(h_DeltaM_physics->FindBin(upperSide_xmin), h_DeltaM_physics->FindBin(upperSide_xmax));
    double expBg_MC=f_Sum_background->Integral(xmin, xmax)*200;
    double expBg_MC_error=expBg_MC/pow(expBg_Sides, 0.5);
    std::cout<<"Expected background from MC = "<<expBg_MC<<" +- "<<expBg_MC_error<<std::endl;
    // Fit Data
    TF1* f_physics=new TF1("f_physics", deltaM_DataFit, 0.1, 0.25, 1);
    f_physics->SetLineColor(kMagenta);
    h_DeltaM_physics->Fit(f_physics, "RLLEFM+");
    double expBg_data=f_physics->Integral(xmin, xmax)*200;
    double expBg_data_error=expBg_data/pow(expBg_Sides, 0.5);
    std::cout<<"Expected background from data = "<<expBg_data<<" +- "<<expBg_data_error<<std::endl;
    
    
    TCanvas *c_DeltaM_Stacked = new TCanvas("c_DeltaM_Stacked");
    THStack *s_DeltaM_sideband=new THStack("s_DeltaM_sideband", ("#deltam Sidebands in "+decay).c_str());
    s_DeltaM_sideband->Add(h_DeltaM_continu, "hist");
    s_DeltaM_sideband->Add(h_DeltaM_generic, "hist");
    s_DeltaM_sideband->Add(h_DeltaM_conver_Dsp, "hist");
    s_DeltaM_sideband->Add(h_DeltaM_conver_Dsm, "hist");
    s_DeltaM_sideband->Add(h_DeltaM_signal, "hist");
    double stack_ymax=h_DeltaM_signal->GetMaximum();
    double physics_ymax=h_DeltaM_physics->GetMaximum();
    if (physics_ymax>stack_ymax)
    {
      s_DeltaM_sideband->SetMaximum(physics_ymax*1.5);
      ymax=physics_ymax;
    }
    else
    {
      ymax=stack_ymax;
    }
    s_DeltaM_sideband->Draw();
    h_DeltaM_physics->Draw("SAME");
    s_DeltaM_sideband->GetXaxis()->SetTitle("#deltam (GeV)");
    s_DeltaM_sideband->GetYaxis()->SetTitle("Number of Events / 5 MeV");
    line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
    line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
    line=new TLine(xmin, ymax, xmin, 0); line->Draw();
    line=new TLine(xmax, ymax, xmax, 0); line->Draw();
    line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
    line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
    TLegend *legendDeltaM=new TLegend(deltam_x1, deltam_y1, deltam_x2, deltam_y2);
    std::string signal_string="Signal MC: "; signal_string+=itoa(signalScale*h_DeltaM_signal->GetEntries()); signal_string+=" Events";
    std::string continu_string="Continuum MC: "; continu_string+=itoa(continuScale*h_DeltaM_continu->GetEntries()); continu_string+=" Events";
    std::string generic_string="Generic MC: "; generic_string+=itoa(genericScale*h_DeltaM_generic->GetEntries()); generic_string+=" Events";
    std::string conver_string="Conversion MC: "; conver_string+=itoa(converScale_Dsp*h_DeltaM_conver_Dsp->GetEntries()+converScale_Dsm*h_DeltaM_conver_Dsm->GetEntries()); conver_string+=" Events";
    std::string physics_string="Data: "; physics_string+=itoa(h_DeltaM_physics->GetEntries()); physics_string+=" Events";
    legendDeltaM->AddEntry(h_DeltaM_signal, signal_string.c_str());
    legendDeltaM->AddEntry(h_DeltaM_continu, continu_string.c_str());
    legendDeltaM->AddEntry(h_DeltaM_generic, generic_string.c_str());
    legendDeltaM->AddEntry(h_DeltaM_conver_Dsp, conver_string.c_str());
    legendDeltaM->AddEntry(h_DeltaM_physics, physics_string.c_str());
    legendDeltaM->SetFillColor(kWhite);
    legendDeltaM->Draw();
    
    std::string filename_eps=decay+"_Sideband_DeltaM_005.eps";
    c_DeltaM_Stacked->SaveAs(filename_eps.c_str());
    
    // Save histograms in a ROOT file
    /*
    std::string filename_root=decay+"_Sideband_DeltaM_005.root";
    TFile *tfile = new TFile(filename_root.c_str(), "RECREATE");
    h_DeltaM_signal->Write();
    h_DeltaM_conver_Dsp->Write();
    h_DeltaM_conver_Dsm->Write();
    h_DeltaM_generic->Write();
    h_DeltaM_continu->Write();
    h_DeltaM_physics->Write();
    tfile->Close();
    */
  }
  
  if (dPhi_sideband)
  {
    std::cout<<"=== dPhi sideband ==="<<std::endl;
    std::cout<<"Signal: ";
    std::cout<<" after signal = "<<signal_afterSignal<<"+-"<<signal_afterSignal_error<<std::endl;
    std::cout<<"Conversion Dsp: ";
    std::cout<<" after signal = "<<conver_Dsp_afterSignal<<"+-"<<conver_Dsp_afterSignal_error<<std::endl;
    std::cout<<"Conversion Dsm: ";
    std::cout<<" after signal = "<<conver_Dsm_afterSignal<<"+-"<<conver_Dsm_afterSignal_error<<std::endl;
    std::cout<<"Generic MC: ";
    std::cout<<" after signal = "<<generic_afterSignal<<"+-"<<generic_afterSignal_error<<std::endl;
    std::cout<<"Continuum MC: ";
    std::cout<<" after signal = "<<continu_afterSignal<<"+-"<<continu_afterSignal_error<<std::endl;
    std::cout<<"Sum MC: ";
    std::cout<<" after signal = "<<mc_afterSignal<<"+-"<<mc_afterSignal_error<<std::endl;
    std::cout<<"Data (electron-fitted): ";
    std::cout<<" after signal = "<<physics_afterSignal<<"+-"<<physics_afterSignal_error<<", Data/MC = "<<dataMC_afterSignal<<"+-"<<dataMC_afterSignal_error<<std::endl;
    std::cout<<"==================="<<std::endl;
    std::cout<<" TeX output "<<std::endl;
    std::cout.precision(2);
    std::cout<<" & "<<conver_Dsp_afterSignal+conver_Dsm_afterSignal<<" $\\pm$ "<<pow(pow(conver_Dsp_afterSignal_error, 2)+pow(conver_Dsm_afterSignal_error, 2), 0.5);    
    std::cout<<" & "<<generic_afterSignal<<" $\\pm$ "<<generic_afterSignal_error;
    std::cout<<" & "<<continu_afterSignal<<" $\\pm$ "<<continu_afterSignal_error;
    std::cout<<" & "<<mc_afterSignal<<" $\\pm$ "<<mc_afterSignal_error;
    std::cout<<" & "<<physics_afterSignal<<" $\\pm$ "<<physics_afterSignal_error;
    std::cout<<" & "<<dataMC_afterSignal<<" $\\pm$ "<<dataMC_afterSignal_error<<" \\\\  \\hline"<<std::endl;
    
    TLine *line;
    TCanvas *dPhi = new TCanvas("dPhi", "", 300, 1200);
    dPhi->Divide(2,3);
    dPhi->cd(1);
    h_dPhi_signal->Draw();
    ymax=(h_dPhi_signal->GetMaximum())*0.75;
    line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
    dPhi->cd(2);
    h_dPhi_conver_Dsp->Draw();
    ymax=(h_dPhi_conver_Dsp->GetMaximum())*0.75;
    line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
    dPhi->cd(3);
    h_dPhi_conver_Dsm->Draw();
    ymax=(h_dPhi_conver_Dsm->GetMaximum())*0.75;
    line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
    dPhi->cd(4);
    h_dPhi_generic->Draw();
    ymax=(h_dPhi_generic->GetMaximum())*0.75;
    line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
    dPhi->cd(5);
    h_dPhi_continu->Draw();
    ymax=(h_dPhi_continu->GetMaximum())*0.75;
    line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
    dPhi->cd(6);
    h_dPhi_physics->Draw();
    ymax=(h_dPhi_physics->GetMaximum())*0.75;
    line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
    
    TCanvas *dPhi_diffD0 = new TCanvas("dPhi_diffD0", "", 300, 1200);
    ymin=diffD0Cut; ymax=0.01;
    dPhi_diffD0->Divide(1,5);
    dPhi_diffD0->cd(1);
    h_dPhi_diffD0_conver_Dsp->Draw("box");
    line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
    line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
    line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
    line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
    dPhi_diffD0->cd(2);
    h_dPhi_diffD0_conver_Dsm->Draw("box");
    line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
    line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
    line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
    line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
    dPhi_diffD0->cd(3);
    h_dPhi_diffD0_generic->Draw("box");
    line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
    line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
    line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
    line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
    dPhi_diffD0->cd(4);
    h_dPhi_diffD0_continu->Draw("box");
    line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
    line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
    line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
    line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
    dPhi_diffD0->cd(5);
    h_dPhi_diffD0_physics->Draw("box");
    line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
    line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
    line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
    line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  }
  
  if (!mBC_sideband && !DeltaM_sideband && !DsMass_sideband && !dPhi_sideband)
  {
    
  }
  
  TCanvas *c_mee_signalMC = new TCanvas("c_mee_signalMC", "c_mee_signalMC");
  h_mee_signalMC->Draw();
  
  return 0;
}
