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
#include "TVector3.h"
#include <set>
#include <map>
#include "time.h"

double pi=3.14159265358979;

bool DeltaM_sideband=false;
bool mBC_sideband=true;
double sideband_size=4;
double sideband_dist=0.5;
double photonMatchingAngle=0.14; // 0.25
double DsMatchingAngle=0.25;

std::string decay="KsKmpipi";
// options: KKpi, KsK, pieta, pietaprime, KKpipi0, pipipi, KsKmpipi, pipi0eta, pietaprimerho
int decayNumber;

double dsPlusMCut_center, dsPlusMCut_range;
double mbcCut_center, mbcCut_range;
double deltaMCut_center, deltaMCut_range;
double branchingFr_mode;

double nConversionSample_Dsp;
double nConversionSample_Dsm;
double genericScale=1/19.2;
double genericScale_error=0.8/(19.2*19.2); // 1/(19.2 +- 1.5) = 0.0521 +- 0.0041
double continuScale=1/5.;
double continuScale_error=0.0;
double nWrongConversionSample_Dsp=4511222;
double nWrongConversionSample_Dsm=4896941;

double luminosity=586; // /pb
double luminosity_error=6;
double prodCrossSection_DsDss=948;
double prodCrossSection_DsDss_error=36;
double branchingFr_Dsstgamma=0.942;

double lowerSide_xmin;
double lowerSide_xmax;
double xmin;
double xmax;
double upperSide_xmin;
double upperSide_xmax;

double deltam_x1, deltam_x2, deltam_y1, deltam_y2;
double mbc_x1, mbc_x2, mbc_y1, mbc_y2;

std::string itoa (int value)
{
	char buffer[50];
	sprintf(buffer,"%d",value);
	std::string str(buffer);
	return str;
}

std::string converFileName_Dsp="/nfs/cor/an2/souvik/MC_gamma_Dsp_";
std::string converFileName_Dsm="/nfs/cor/an2/souvik/MC_gamma_Dsm_";

void setValues()
{

  if (decay=="KKpi")
  {
    decayNumber=401;
    
    // Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.011;
    mbcCut_center=2.112; mbcCut_range=0.004;
    //deltaMCut_center=0.1438; deltaMCut_range=0.006;
    deltaMCut_center=0.14; deltaMCut_range=0.02; // Widened for Ds gamma
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
  }
  else if (decay=="KsK")
  {
    decayNumber=400;
    
    // Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.008;
    mbcCut_center=2.112; mbcCut_range=0.007;
    // deltaMCut_center=0.1438; deltaMCut_range=0.006;
    deltaMCut_center=0.14; deltaMCut_range=0.02; // Widened for Ds gamma
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.35; mbc_y2=0.6;
  }
  else if (decay=="pieta")
  {
    decayNumber=440;
    
    //Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.016;
    mbcCut_center=2.112; mbcCut_range=0.008;
    // deltaMCut_center=0.1438; deltaMCut_range=0.008;
    deltaMCut_center=0.14; deltaMCut_range=0.02; // Widened for Ds gamma
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
  }
  else if (decay=="pietaprime")
  {
    decayNumber=460;
    
    //Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.011;
    mbcCut_center=2.112; mbcCut_range=0.011;
    // deltaMCut_center=0.1438; deltaMCut_range=0.013;
    deltaMCut_center=0.14; deltaMCut_range=0.02; // Widened for Ds gamma
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
  }
  else if (decay=="KKpipi0")
  {
    decayNumber=404;
    
    //Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.01;
    mbcCut_center=2.112; mbcCut_range=0.004;
    //deltaMCut_center=0.1438; deltaMCut_range=0.006;
    deltaMCut_center=0.14; deltaMCut_range=0.02; // Widened for Ds gamma
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
  }
  else if (decay=="pipipi")
  {
    decayNumber=421;
    
    //Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.012;
    mbcCut_center=2.112; mbcCut_range=0.004;
    //deltaMCut_center=0.1438; deltaMCut_range=0.006;
    deltaMCut_center=0.14; deltaMCut_range=0.02; // Widened for Ds gamma
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
  }
  else if (decay=="KsKmpipi")
  {
    decayNumber=406;
    
    //Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.006;
    mbcCut_center=2.112; mbcCut_range=0.005;
    //deltaMCut_center=0.1438; deltaMCut_range=0.008;
    deltaMCut_center=0.14; deltaMCut_range=0.02; // Widened for Ds gamma
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
  }
  else if (decay=="pipi0eta")
  {
    decayNumber=441;
    
    //Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.015;
    mbcCut_center=2.112; mbcCut_range=0.004;
    //deltaMCut_center=0.1438; deltaMCut_range=0.005;
    deltaMCut_center=0.14; deltaMCut_range=0.02; // Widened for Ds gamma
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
  }
  else if (decay=="pietaprimerho")
  {
    decayNumber=480;
    
    //Optimized Cuts
    dsPlusMCut_center=1.96849; dsPlusMCut_range=0.012;
    mbcCut_center=2.112; mbcCut_range=0.004;
    //deltaMCut_center=0.1438; deltaMCut_range=0.007;
    deltaMCut_center=0.14; deltaMCut_range=0.02; // Widened for Ds gamma
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
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
    upperSide_xmax=2.150;
  }
}

void setFileNames()
{
  if (decay=="KKpi")
  {
    converFileName_Dsp=converFileName_Dsp+"KKpi/DsTaggedDecaysProc_MC_gamma_Dsp_KKpi.root";
    converFileName_Dsm=converFileName_Dsm+"KKpi/DsTaggedDecaysProc_MC_gamma_Dsm_KKpi.root";
    nConversionSample_Dsp=99880;
    nConversionSample_Dsm=99880;
  }
  else if (decay=="KsK")
  {
    converFileName_Dsp=converFileName_Dsp+"KsK/DsTaggedDecaysProc_MC_gamma_Dsp_KsK.root";
    converFileName_Dsm=converFileName_Dsm+"KsK/DsTaggedDecaysProc_MC_gamma_Dsm_KsK.root";
    nConversionSample_Dsp=938239;
    nConversionSample_Dsm=974094;
  }
  else if (decay=="pieta")
  {
    converFileName_Dsp=converFileName_Dsp+"pieta/DsTaggedDecaysProc_MC_gamma_Dsp_pieta.root";
    converFileName_Dsm=converFileName_Dsm+"pieta/DsTaggedDecaysProc_MC_gamma_Dsm_pieta.root";
    nConversionSample_Dsp=99880;
    nConversionSample_Dsm=99880;
  }
  else if (decay=="pietaprime")
  {
    converFileName_Dsp=converFileName_Dsp+"pietaprime/DsTaggedDecaysProc_MC_gamma_Dsp_pietaprime.root";
    converFileName_Dsm=converFileName_Dsm+"pietaprime/DsTaggedDecaysProc_MC_gamma_Dsm_pietaprime.root";
    nConversionSample_Dsp=99880;
    nConversionSample_Dsm=79904;
  }
  else if (decay=="KKpipi0")
  {
    converFileName_Dsp=converFileName_Dsp+"KKpipi0/DsTaggedDecaysProc_MC_gamma_Dsp_KKpipi0.root";
    converFileName_Dsm=converFileName_Dsm+"KKpipi0/DsTaggedDecaysProc_MC_gamma_Dsm_KKpipi0.root";
    nConversionSample_Dsp=96355;
    nConversionSample_Dsm=99880;
  }
  else if (decay=="pipipi")
  {
    converFileName_Dsp=converFileName_Dsp+"pipipi/DsTaggedDecaysProc_MC_gamma_Dsp_pipipi.root";
    converFileName_Dsm=converFileName_Dsm+"pipipi/DsTaggedDecaysProc_MC_gamma_Dsm_pipipi.root";
    nConversionSample_Dsp=99880;
    nConversionSample_Dsm=99880;
  }
  else if (decay=="KsKmpipi")
  {
    converFileName_Dsp=converFileName_Dsp+"KsKmpipi/DsTaggedDecaysProc_MC_gamma_Dsp_KsKmpipi.root";
    converFileName_Dsm=converFileName_Dsm+"KsKmpipi/DsTaggedDecaysProc_MC_gamma_Dsm_KsKmpipi.root";
    nConversionSample_Dsp=99880;
    nConversionSample_Dsm=99880;
  }
  else if (decay=="pipi0eta")
  {
    converFileName_Dsp=converFileName_Dsp+"pipi0eta/DsTaggedDecaysProc_MC_gamma_Dsp_pipi0eta.root";
    converFileName_Dsm=converFileName_Dsm+"pipi0eta/DsTaggedDecaysProc_MC_gamma_Dsm_pipi0eta.root";
    nConversionSample_Dsp=99880;
    nConversionSample_Dsm=99880;
  }
  else if (decay=="pietaprimerho")
  {
    converFileName_Dsp=converFileName_Dsp+"pietaprimerho/DsTaggedDecaysProc_MC_gamma_Dsp_pietaprimerho.root";
    converFileName_Dsm=converFileName_Dsm+"pietaprimerho/DsTaggedDecaysProc_MC_gamma_Dsm_pietaprimerho.root";
    nConversionSample_Dsp=99880;
    nConversionSample_Dsm=99321;
  }
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

bool x925Cut(float x925)
{
  return (x925>1.0);
}

struct NEvents {
  float noCut;
  float tagCut;
  float dsPlusMCut;
  float mbcCut;
  float deltaMCut;
};

struct EventNumber {
  float noCut_run;      float noCut_event;
  float tagCut_run;     float tagCut_event;
  float dsPlusMCut_run; float dsPlusMCut_event;
  float mbcCut_run;     float mbcCut_event;
  float deltaMCut_run;  float deltaMCut_event;
};

int DsGammaAnalysis()
{
  
  setValues();
  setFileNames();
  setBounds();
  
  TChain *converTree_Dsp=new TChain("DsTaggedDecaysProc/nt6");
  converTree_Dsp->Add(converFileName_Dsp.c_str());
  
  TChain *converTree_Dsm=new TChain("DsTaggedDecaysProc/nt6");
  converTree_Dsm->Add(converFileName_Dsm.c_str());
  
  TChain *genericTree = new TChain("DsTaggedDecaysProc/nt6");
  genericTree->Add("/nfs/cor/an3/souvik/Dataset39/DsTaggedProc_GenericMC_213586_214863_backup.root");
  genericTree->Add("/nfs/cor/an3/souvik/Dataset40/DsTaggedProc_GenericMC_215307_217385_backup.root");
  genericTree->Add("/nfs/cor/an3/souvik/Dataset41/DsTaggedProc_GenericMC_217687_219721_backup.root");
  genericTree->Add("/nfs/cor/an3/souvik/Dataset47/DsTaggedProc_GenericMC_230474_232255_backup.root");
  genericTree->Add("/nfs/cor/an3/souvik/Dataset48/DsTaggedProc_GenericMC_232264_234607_backup.root");
  
  TChain *continuTree = new TChain("DsTaggedDecaysProc/nt6");
  continuTree->Add("/nfs/cor/an3/souvik/Dataset39/DsTaggedProc_ContinuumMC_213586_214863.root");
  continuTree->Add("/nfs/cor/an3/souvik/Dataset40/DsTaggedProc_ContinuumMC_215307_217385.root");
  continuTree->Add("/nfs/cor/an3/souvik/Dataset41/DsTaggedProc_ContinuumMC_217687_219721.root");
  continuTree->Add("/nfs/cor/an3/souvik/Dataset47/DsTaggedProc_ContinuumMC_230474_232255.root");
  continuTree->Add("/nfs/cor/an3/souvik/Dataset48/DsTaggedProc_ContinuumMC_232264_234607.root");
  
  TChain *physicsTree = new TChain("DsTaggedDecaysProc/nt6");
  physicsTree->Add("/nfs/cor/an3/souvik/Dataset39/DsTaggedDecaysProc_ReTaggedData_213586_214863.root");
  physicsTree->Add("/nfs/cor/an3/souvik/Dataset40/DsTaggedDecaysProc_ReTaggedData_215307_217385.root");
  physicsTree->Add("/nfs/cor/an3/souvik/Dataset41/DsTaggedDecaysProc_ReTaggedData_217687_219721.root");
  physicsTree->Add("/nfs/cor/an3/souvik/Dataset47/DsTaggedDecaysProc_ReTaggedData_230474_232255.root");
  physicsTree->Add("/nfs/cor/an3/souvik/Dataset48/DsTaggedDecaysProc_ReTaggedData_232264_234607.root");
  
  TChain *wrongConverTree_Dsp=new TChain("DsTaggedDecaysProc/nt6");
  wrongConverTree_Dsp->Add("/nfs/cor/an3/souvik/MC_gamma_Dsp_generic/DsTaggedDecaysProc_MC_gamma_Dsp_generic.root");
  
  TChain *wrongConverTree_Dsm=new TChain("DsTaggedDecaysProc/nt6");
  wrongConverTree_Dsm->Add("/nfs/cor/an3/souvik/MC_gamma_Dsm_generic/DsTaggedDecaysProc_MC_gamma_Dsm_generic.root");
  
  NEvents converEvents_Dsp={0,0,0,0,0};
  NEvents converEvents_Dsm={0,0,0,0,0};
  NEvents genericEvents={0,0,0,0,0};
  NEvents continuEvents={0,0,0,0,0};
  NEvents physicsEvents={0,0,0,0,0};
  
  EventNumber converNumber_Dsp={0,0,0,0,0,0,0,0,0,0};
  EventNumber converNumber_Dsm={0,0,0,0,0,0,0,0,0,0};
  EventNumber genericNumber={0,0,0,0,0,0,0,0,0,0};
  EventNumber continuNumber={0,0,0,0,0,0,0,0,0,0};
  EventNumber physicsNumber={0,0,0,0,0,0,0,0,0,0};
  
  float runNumber_conver_Dsp, eventNumber_conver_Dsp;
  float runNumber_conver_Dsm, eventNumber_conver_Dsm;
  float runNumber_generic, eventNumber_generic;
  float runNumber_continu, eventNumber_continu;
  float runNumber_physics, eventNumber_physics;
  float runNumber_wrongConver_Dsp, eventNumber_wrongConver_Dsp;
  float runNumber_wrongConver_Dsm, eventNumber_wrongConver_Dsm;
  float beamEnergy_wrongConver_Dsp, beamEnergy_wrongConver_Dsm;
  float dsPlusM_conver_Dsp, dsPlusCharge_conver_Dsp, MBC_conver_Dsp, DeltaM_conver_Dsp, decayMode_conver_Dsp, e9oe25Unf_conver_Dsp;
  float dsPlusM_conver_Dsm, dsPlusCharge_conver_Dsm, MBC_conver_Dsm, DeltaM_conver_Dsm, decayMode_conver_Dsm, e9oe25Unf_conver_Dsm;
  float dsPlusM_generic, dsPlusCharge_generic, MBC_generic, DeltaM_generic, decayMode_generic, e9oe25Unf_generic;
  float dsPlusM_continu, dsPlusCharge_continu, MBC_continu, DeltaM_continu, decayMode_continu, e9oe25Unf_continu;
  float dsPlusM_physics, dsPlusCharge_physics, MBC_physics, DeltaM_physics, decayMode_physics, e9oe25Unf_physics;
  float dsPlusM_wrongConver_Dsp, dsPlusCharge_wrongConver_Dsp, MBC_wrongConver_Dsp, DeltaM_wrongConver_Dsp, decayMode_wrongConver_Dsp, e9oe25Unf_wrongConver_Dsp;
  float dsPlusM_wrongConver_Dsm, dsPlusCharge_wrongConver_Dsm, MBC_wrongConver_Dsm, DeltaM_wrongConver_Dsm, decayMode_wrongConver_Dsm, e9oe25Unf_wrongConver_Dsm;
  float photonE_conver_Dsp, photonE_conver_Dsm;
  float photonE_conver_Dsp_MC, photonE_conver_Dsm_MC;
  float dsPlusPx_reco_conver_Dsp, dsPlusPy_reco_conver_Dsp, dsPlusPz_reco_conver_Dsp;
  float dsPlusPx_MC_conver_Dsp, dsPlusPy_MC_conver_Dsp, dsPlusPz_MC_conver_Dsp;
  float photonPx_reco_conver_Dsp, photonPy_reco_conver_Dsp, photonPz_reco_conver_Dsp;
  float photonPx_MC_conver_Dsp, photonPy_MC_conver_Dsp, photonPz_MC_conver_Dsp;
  float dsPlusPx_reco_conver_Dsm, dsPlusPy_reco_conver_Dsm, dsPlusPz_reco_conver_Dsm;
  float dsPlusPx_MC_conver_Dsm, dsPlusPy_MC_conver_Dsm, dsPlusPz_MC_conver_Dsm;
  float photonPx_reco_conver_Dsm, photonPy_reco_conver_Dsm, photonPz_reco_conver_Dsm;
  float photonPx_MC_conver_Dsm, photonPy_MC_conver_Dsm, photonPz_MC_conver_Dsm;
  float dsPlusPx_reco_wrongConver_Dsp, dsPlusPy_reco_wrongConver_Dsp, dsPlusPz_reco_wrongConver_Dsp;
  float dsPlusPx_MC_wrongConver_Dsp, dsPlusPy_MC_wrongConver_Dsp, dsPlusPz_MC_wrongConver_Dsp;
  float photonE_reco_wrongConver_Dsp, photonPx_reco_wrongConver_Dsp, photonPy_reco_wrongConver_Dsp, photonPz_reco_wrongConver_Dsp;
  float photonPx_MC_wrongConver_Dsp, photonPy_MC_wrongConver_Dsp, photonPz_MC_wrongConver_Dsp;
  float dsPlusPx_reco_wrongConver_Dsm, dsPlusPy_reco_wrongConver_Dsm, dsPlusPz_reco_wrongConver_Dsm;
  float dsPlusPx_MC_wrongConver_Dsm, dsPlusPy_MC_wrongConver_Dsm, dsPlusPz_MC_wrongConver_Dsm;
  float photonE_reco_wrongConver_Dsm, photonPx_reco_wrongConver_Dsm, photonPy_reco_wrongConver_Dsm, photonPz_reco_wrongConver_Dsm;
  float photonPx_MC_wrongConver_Dsm, photonPy_MC_wrongConver_Dsm, photonPz_MC_wrongConver_Dsm;
  float conversionBit_generic;
  
  // Variables used to count sideband and signal entries
  double conver_Dsp_beforeSignal=0, conver_Dsp_afterSignal=0, conver_Dsp_atSignal=0;
  double conver_Dsm_beforeSignal=0, conver_Dsm_afterSignal=0, conver_Dsm_atSignal=0;
  double generic_beforeSignal=0, generic_afterSignal=0, generic_atSignal=0;
  double continu_beforeSignal=0, continu_afterSignal=0, continu_atSignal=0;
  double physics_beforeSignal=0, physics_afterSignal=0, physics_atSignal=0;
  
  TH1D *h_dsPlusM_conver_Dsp = new TH1D("h_dsPlusM_conver_Dsp", "m_{D_{S}^{+}} conver_Dsp Sample; GeV", 100, 1.9, 2.05); h_dsPlusM_conver_Dsp->SetLineColor(kRed);
  TH1D *h_dsPlusM_conver_Dsm = new TH1D("h_dsPlusM_conver_Dsm", "m_{D_{S}^{-}} conver_Dsm Sample; GeV", 100, 1.9, 2.05); h_dsPlusM_conver_Dsm->SetLineColor(kRed);
  TH1D *h_dsPlusM_conver = new TH1D("h_dsPlusM_conver", "m_{D_{S}^{#pm}} Conversion Sample; m_{D_{S}^{#pm}} (GeV); # Events", 100, 1.9, 2.05); h_dsPlusM_conver->SetLineColor(kRed);
  TH1D *h_dsPlusM_generic = new TH1D("h_dsPlusM_generic", "m_{D_{S}^{+}} generic Sample; GeV", 100, 1.9, 2.1); h_dsPlusM_generic->SetLineColor(kRed);
  TH1D *h_dsPlusM_continu = new TH1D("h_dsPlusM_continu", "m_{D_{S}^{+}} continuum Background Sample; GeV", 100, 1.9, 2.1); h_dsPlusM_continu->SetLineColor(kBlue);
  TH1D *h_dsPlusM_physics = new TH1D("h_dsPlusM_physics", "m_{D_{S}^{+}} Data; GeV", 100, 1.9, 2.1); h_dsPlusM_physics->SetLineColor(kGreen);
  
  TH1D *h_MBC_conver_Dsp = new TH1D("h_MBC_conver_Dsp", "m_{BC} conver_Dsp Sample; GeV; # Events", 120, 2.04, 2.16); h_MBC_conver_Dsp->SetLineColor(kRed); h_MBC_conver_Dsp->Sumw2();
  TH1D *h_MBC_conver_Dsm = new TH1D("h_MBC_conver_Dsm", "m_{BC} conver_Dsm Sample; GeV; # Events", 120, 2.04, 2.16); h_MBC_conver_Dsm->SetLineColor(kRed); h_MBC_conver_Dsm->Sumw2();
  TH1D *h_MBC_conver = new TH1D("h_MBC_conver", "m_{BC} Conversion MC Sample; m_{BC} (GeV); # Events", 120, 2.04, 2.16); h_MBC_conver->SetLineColor(kRed); h_MBC_conver->Sumw2();
  TH1D *h_MBC_generic = new TH1D("h_MBC_generic", "m_{BC} generic Sample; GeV", 120, 2.04, 2.16); h_MBC_generic->SetLineColor(kGreen); h_MBC_generic->Sumw2();
  TH1D *h_MBC_generic_veto = new TH1D("h_MBC_generic_veto", "m_{BC} generic Sample excluding Ds gamma Events; GeV", 120, 2.04, 2.16); h_MBC_generic_veto->SetLineColor(kGreen); h_MBC_generic_veto->Sumw2();
  TH1D *h_MBC_generic_wrong = new TH1D("h_MBC_generic_wrong", "m_{BC} generic Sample with wrong-sign Ds reconstructed; GeV", 120, 2.04, 2.16); h_MBC_generic_wrong->Sumw2();
  TH1D *h_MBC_continu = new TH1D("h_MBC_continu", "m_{BC} continuum Background; GeV", 120, 2.04, 2.16); h_MBC_continu->SetLineColor(kBlue); h_MBC_continu->Sumw2();
  TH1D *h_MBC_physics = new TH1D("h_MBC_physics", "m_{BC} Data; GeV", 120, 2.04, 2.16); h_MBC_physics->SetLineColor(kMagenta); h_MBC_physics->Sumw2();
  TH1D *h_MBC_wrongConver_Dsp = new TH1D("h_MBC_wrongConver_Dsp", "m_{BC} wrongConver_Dsp Sample; GeV; # Events", 120, 2.04, 2.16); h_MBC_wrongConver_Dsp->SetLineColor(kRed); h_MBC_wrongConver_Dsp->Sumw2();
  TH1D *h_MBC_wrongConver_Dsm = new TH1D("h_MBC_wrongConver_Dsm", "m_{BC} wrongConver_Dsm Sample; GeV; # Events", 120, 2.04, 2.16); h_MBC_wrongConver_Dsm->SetLineColor(kRed); h_MBC_wrongConver_Dsm->Sumw2();
  TH1D *h_MBC_wrongConver = new TH1D("h_MBC_wrongConver", "m_{BC} wrongConversion MC Sample; m_{BC} (GeV); # Events", 120, 2.04, 2.16); h_MBC_wrongConver->SetLineColor(kRed); h_MBC_wrongConver->Sumw2();
  
  // Strictly matched, i.e. Ds and photon matched
  TH1D *h_MBC_conver_Dsp_matched = new TH1D("h_MBC_conver_Dsp_matched", "m_{BC} conver_Dsp Sample (Matched); GeV; # Events", 120, 2.04, 2.16); h_MBC_conver_Dsp_matched->SetLineColor(kBlack); h_MBC_conver_Dsp_matched->Sumw2();
  TH1D *h_MBC_conver_Dsm_matched = new TH1D("h_MBC_conver_Dsm_matched", "m_{BC} conver_Dsm Sample (Matched); GeV; # Events", 120, 2.04, 2.16); h_MBC_conver_Dsm_matched->SetLineColor(kBlack); h_MBC_conver_Dsm_matched->Sumw2();
  TH1D *h_MBC_conver_matched = new TH1D("h_MBC_conver_matched", "m_{BC} Conversion MC Sample (Matched); m_{BC} (GeV); # Events", 120, 2.04, 2.16); h_MBC_conver_matched->SetLineColor(kBlack); h_MBC_conver_matched->Sumw2();
  // Ds matched, Photon unmatched
  TH1D *h_MBC_conver_Dsp_unmatched = new TH1D("h_MBC_conver_Dsp_unmatched", "m_{BC} conver_Dsp Sample (Unmatched); GeV; # Events", 120, 2.04, 2.16); h_MBC_conver_Dsp_unmatched->SetLineColor(kBlack); h_MBC_conver_Dsp_unmatched->Sumw2();
  TH1D *h_MBC_conver_Dsm_unmatched = new TH1D("h_MBC_conver_Dsm_unmatched", "m_{BC} conver_Dsm Sample (Unmatched); GeV; # Events", 120, 2.04, 2.16); h_MBC_conver_Dsm_unmatched->SetLineColor(kBlack); h_MBC_conver_Dsm_unmatched->Sumw2();
  TH1D *h_MBC_conver_unmatched = new TH1D("h_MBC_conver_unmatched", "m_{BC} Conversion MC Sample (Unmatched); GeV; # Events", 120, 2.04, 2.16); h_MBC_conver_unmatched->SetLineColor(kBlack); h_MBC_conver_unmatched->Sumw2();
  
  // Strictly matched, i.e. Ds matched & photon matched
  TH1D *h_MBC_wrongConver_Dsp_matched = new TH1D("h_MBC_wrongConver_Dsp_matched", "m_{BC} conver_Dsp Sample (Matched); GeV; # Events", 120, 2.04, 2.16); h_MBC_wrongConver_Dsp_matched->SetLineColor(kBlack); h_MBC_wrongConver_Dsp_matched->Sumw2();
  TH1D *h_MBC_wrongConver_Dsm_matched = new TH1D("h_MBC_wrongConver_Dsm_matched", "m_{BC} conver_Dsm Sample (Matched); GeV; # Events", 120, 2.04, 2.16); h_MBC_wrongConver_Dsm_matched->SetLineColor(kBlack); h_MBC_wrongConver_Dsm_matched->Sumw2();
  TH1D *h_MBC_wrongConver_matched = new TH1D("h_MBC_wrongConver_matched", "m_{BC} Conversion MC Sample (Matched); m_{BC} (GeV); # Events", 120, 2.04, 2.16); h_MBC_wrongConver_matched->SetLineColor(kBlack); h_MBC_wrongConver_matched->Sumw2();
  // Photon unmatched
  TH1D *h_MBC_wrongConver_Dsp_unmatched = new TH1D("h_MBC_wrongConver_Dsp_unmatched", "m_{BC} conver_Dsp Sample (Unmatched); GeV; # Events", 120, 2.04, 2.16); h_MBC_wrongConver_Dsp_unmatched->SetLineColor(kBlack); h_MBC_wrongConver_Dsp_unmatched->Sumw2();
  TH1D *h_MBC_wrongConver_Dsm_unmatched = new TH1D("h_MBC_wrongConver_Dsm_unmatched", "m_{BC} conver_Dsm Sample (Unmatched); GeV; # Events", 120, 2.04, 2.16); h_MBC_wrongConver_Dsm_unmatched->SetLineColor(kBlack); h_MBC_wrongConver_Dsm_unmatched->Sumw2();
  TH1D *h_MBC_wrongConver_unmatched = new TH1D("h_MBC_wrongConver_unmatched", "m_{BC} Conversion MC Sample (Unmatched); GeV; # Events", 120, 2.04, 2.16); h_MBC_wrongConver_unmatched->SetLineColor(kBlack); h_MBC_wrongConver_unmatched->Sumw2();
  // Photon and Ds unmatched
  TH1D *h_MBC_wrongConver_Dsp_strictUnmatched = new TH1D("h_MBC_wrongConver_Dsp_strictUnmatched", "h_MBC_wrongConver_Dsp_strictUnmatched", 120, 2.04, 2.16); h_MBC_wrongConver_Dsp_strictUnmatched->Sumw2();
  TH1D *h_MBC_wrongConver_Dsm_strictUnmatched = new TH1D("h_MBC_wrongConver_Dsm_strictUnmatched", "h_MBC_wrongConver_Dsm_strictUnmatched", 120, 2.04, 2.16); h_MBC_wrongConver_Dsm_strictUnmatched->Sumw2();
  TH1D *h_MBC_wrongConver_strictUnmatched = new TH1D("h_MBC_wrongConver_strictUnmatched", "h_MBC_wrongConver_Dsm_strictUnmatched", 120, 2.04, 2.16); h_MBC_wrongConver_strictUnmatched->Sumw2();
  
  TH2D *h_photonE_MBC_conver = new TH2D("h_photonE_MBC_conver", "h_photonE_MBC_conver", 50, 2.04, 2.16, 50, -0.1, 0.1);
  TH1D *h_photonE_wrongConver_Dsp_unmatched = new TH1D("h_photonE_wrongConver_Dsp_unmatched", "h_photonE_wrongConver_Dsp_unmatched", 100, 0., 0.5); h_photonE_wrongConver_Dsp_unmatched->Sumw2();
  TH1D *h_photonE_wrongConver_Dsm_unmatched = new TH1D("h_photonE_wrongConver_Dsm_unmatched", "h_photonE_wrongConver_Dsm_unmatched", 100, 0., 0.5); h_photonE_wrongConver_Dsm_unmatched->Sumw2();
  TH1D *h_photonE_wrongConver_unmatched = new TH1D("h_photonE_wrongConver_unmatched", "h_photonE_wrongConver_unmatched", 100, 0., 0.5); h_photonE_wrongConver_unmatched->Sumw2();
  TH1D *h_photonE_wrongConver_Dsp_matched = new TH1D("h_photonE_wrongConver_Dsp_matched", "h_photonE_wrongConver_Dsp_matched", 100, 0., 0.5); h_photonE_wrongConver_Dsp_matched->Sumw2();
  TH1D *h_photonE_wrongConver_Dsm_matched = new TH1D("h_photonE_wrongConver_Dsm_matched", "h_photonE_wrongConver_Dsm_matched", 100, 0., 0.5); h_photonE_wrongConver_Dsm_matched->Sumw2();
  TH1D *h_photonE_wrongConver_matched = new TH1D("h_photonE_wrongConver_matched", "h_photonE_wrongConver_matched", 100, 0., 0.5); h_photonE_wrongConver_matched->Sumw2();
  
  TH1D *h_DeltaM_conver_Dsp = new TH1D("h_DeltaM_conver_Dsp", "#deltaM conver_Dsp Sample; GeV; # Events", 150, 0.05, 0.2); h_DeltaM_conver_Dsp->SetLineColor(kRed);
  TH1D *h_DeltaM_conver_Dsm = new TH1D("h_DeltaM_conver_Dsm", "#deltaM conver_Dsm Sample; GeV; # Events", 150, 0.05, 0.2); h_DeltaM_conver_Dsm->SetLineColor(kRed);
  TH1D *h_DeltaM_conver = new TH1D("h_DeltaM_conver", "#deltam Conversion MC Sample; #deltam (GeV); # Events", 150, 0.05, 0.2); h_DeltaM_conver->SetLineColor(kRed);
  TH1D *h_DeltaM_generic = new TH1D("h_DeltaM_generic", "#delta m generic Sample; #deltam (GeV); Number of Events", 150, 0.05, 0.2); h_DeltaM_generic->SetLineColor(kGreen); h_DeltaM_generic->Sumw2();
  TH1D *h_DeltaM_generic_veto = new TH1D("h_DeltaM_generic_veto", "#delta m generic Sample excluding Ds gamma Events; #deltam (GeV); Number of Events", 150, 0.05, 0.2); h_DeltaM_generic_veto->SetLineColor(kGreen); h_DeltaM_generic_veto->Sumw2();
  TH1D *h_DeltaM_generic_wrong = new TH1D("h_DeltaM_generic_wrong", "#delta m generic Sample with wrong-sign Ds reconstructed; #delta m (GeV)", 150, 0.05, 0.2); h_DeltaM_generic_wrong->Sumw2();
  TH1D *h_DeltaM_continu = new TH1D("h_DeltaM_continu", "#delta m continuum Background Sample; #deltam (GeV); Number of Events", 150, 0.05, 0.2); h_DeltaM_continu->SetLineColor(kBlue); h_DeltaM_continu->Sumw2();
  TH1D *h_DeltaM_physics = new TH1D("h_DeltaM_physics", "#delta m Data; #deltam (GeV); Number of Events", 150, 0.05, 0.2); h_DeltaM_physics->SetLineColor(kMagenta); h_DeltaM_physics->Sumw2();
  TH1D *h_DeltaM_wrongConver_Dsp = new TH1D("h_DeltaM_wrongConver_Dsp", "m_{BC} wrongConver_Dsp Sample; GeV; # Events", 150, 0.05, 0.2); h_DeltaM_wrongConver_Dsp->SetLineColor(kRed); h_DeltaM_wrongConver_Dsp->Sumw2();
  TH1D *h_DeltaM_wrongConver_Dsm = new TH1D("h_DeltaM_wrongConver_Dsm", "m_{BC} wrongConver_Dsm Sample; GeV; # Events", 150, 0.05, 0.2); h_DeltaM_wrongConver_Dsm->SetLineColor(kRed); h_DeltaM_wrongConver_Dsm->Sumw2();
  TH1D *h_DeltaM_wrongConver = new TH1D("h_DeltaM_wrongConver", "m_{BC} wrongConversion MC Sample; m_{BC} (GeV); # Events", 150, 0.05, 0.2); h_DeltaM_wrongConver->SetLineColor(kRed); h_DeltaM_wrongConver->Sumw2();
  
  TH1D *h_e9oe25Unf_conver_Dsp = new TH1D("h_e9oe25Unf_conver_Dsp", "E9/E25 Unfolded Data", 100, 0.5, 2.); h_e9oe25Unf_conver_Dsp->SetLineColor(kRed);
  TH1D *h_e9oe25Unf_conver_Dsm = new TH1D("h_e9oe25Unf_conver_Dsm", "E9/E25 Unfolded Data", 100, 0.5, 2.); h_e9oe25Unf_conver_Dsm->SetLineColor(kRed);
  TH1D *h_e9oe25Unf_generic = new TH1D("h_e9oe25Unf_generic", "E9/E25 Unfolded Data", 100, 0.5, 2.); h_e9oe25Unf_generic->SetLineColor(kGreen);
  TH1D *h_e9oe25Unf_continu = new TH1D("h_e9oe25Unf_continu", "E9/E25 Unfolded Data", 100, 0.5, 2.); h_e9oe25Unf_continu->SetLineColor(kBlue);
  TH1D *h_e9oe25Unf_physics = new TH1D("h_e9oe25Unf_physics", "E9/E25 Unfolded Data", 100, 0.5, 2.); h_e9oe25Unf_physics->SetLineColor(kMagenta);
  
  TH1D *h_angle_photon = new TH1D("h_angle_photon", "h_angle_photon", 100, 0., pi);
  TH1D *h_angle_wrongDs = new TH1D("h_angle_wrongDs", "h_angle_wrongDs", 100, 0., pi);
  
  converTree_Dsp->SetBranchAddress("Run", &(runNumber_conver_Dsp));
  converTree_Dsp->SetBranchAddress("Event", &(eventNumber_conver_Dsp));
  converTree_Dsp->SetBranchAddress("dsPlusM", &(dsPlusM_conver_Dsp));
  converTree_Dsp->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_conver_Dsp));
  converTree_Dsp->SetBranchAddress("DecayMode", &(decayMode_conver_Dsp));
  converTree_Dsp->SetBranchAddress("MBC", &(MBC_conver_Dsp));
  converTree_Dsp->SetBranchAddress("DeltaM", &(DeltaM_conver_Dsp));
  converTree_Dsp->SetBranchAddress("e9oe25Unf", &(e9oe25Unf_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kPhotonE_reco", &(photonE_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kPhotonE_MC", &(photonE_conver_Dsp_MC));
  // Need these for MC truth matching
  converTree_Dsp->SetBranchAddress("kDsPlusPx_reco", &(dsPlusPx_reco_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kDsPlusPy_reco", &(dsPlusPy_reco_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kDsPlusPz_reco", &(dsPlusPz_reco_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kDsPlusPx_MC", &(dsPlusPx_MC_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kDsPlusPy_MC", &(dsPlusPy_MC_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kDsPlusPz_MC", &(dsPlusPz_MC_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kPhotonPx_reco", &(photonPx_reco_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kPhotonPy_reco", &(photonPy_reco_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kPhotonPz_reco", &(photonPz_reco_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kPhotonPx_MC", &(photonPx_MC_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kPhotonPy_MC", &(photonPy_MC_conver_Dsp));
  converTree_Dsp->SetBranchAddress("kPhotonPz_MC", &(photonPz_MC_conver_Dsp));
  
  converTree_Dsm->SetBranchAddress("Run", &(runNumber_conver_Dsm));
  converTree_Dsm->SetBranchAddress("Event", &(eventNumber_conver_Dsm));
  converTree_Dsm->SetBranchAddress("dsPlusM", &(dsPlusM_conver_Dsm));
  converTree_Dsm->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_conver_Dsm));
  converTree_Dsm->SetBranchAddress("DecayMode", &(decayMode_conver_Dsm));
  converTree_Dsm->SetBranchAddress("MBC", &(MBC_conver_Dsm));
  converTree_Dsm->SetBranchAddress("DeltaM", &(DeltaM_conver_Dsm));
  converTree_Dsm->SetBranchAddress("e9oe25Unf", &(e9oe25Unf_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kPhotonE_reco", &(photonE_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kPhotonE_MC", &(photonE_conver_Dsm_MC));
  // Need these for MC truth matching
  converTree_Dsm->SetBranchAddress("kDsPlusPx_reco", &(dsPlusPx_reco_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kDsPlusPy_reco", &(dsPlusPy_reco_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kDsPlusPz_reco", &(dsPlusPz_reco_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kDsPlusPx_MC", &(dsPlusPx_MC_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kDsPlusPy_MC", &(dsPlusPy_MC_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kDsPlusPz_MC", &(dsPlusPz_MC_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kPhotonPx_reco", &(photonPx_reco_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kPhotonPy_reco", &(photonPy_reco_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kPhotonPz_reco", &(photonPz_reco_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kPhotonPx_MC", &(photonPx_MC_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kPhotonPy_MC", &(photonPy_MC_conver_Dsm));
  converTree_Dsm->SetBranchAddress("kPhotonPz_MC", &(photonPz_MC_conver_Dsm));
  
  genericTree->SetBranchAddress("Run", &(runNumber_generic));
  genericTree->SetBranchAddress("Event", &(eventNumber_generic));
  genericTree->SetBranchAddress("dsPlusM", &(dsPlusM_generic));
  genericTree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_generic));
  genericTree->SetBranchAddress("DecayMode", &(decayMode_generic));
  genericTree->SetBranchAddress("MBC", &(MBC_generic));
  genericTree->SetBranchAddress("DeltaM", &(DeltaM_generic));
  genericTree->SetBranchAddress("e9oe25Unf", &(e9oe25Unf_generic));
  genericTree->SetBranchAddress("conversionBit", &(conversionBit_generic));
  
  continuTree->SetBranchAddress("Run", &(runNumber_continu));
  continuTree->SetBranchAddress("Event", &(eventNumber_continu));
  continuTree->SetBranchAddress("dsPlusM", &(dsPlusM_continu));
  continuTree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_continu));
  continuTree->SetBranchAddress("DecayMode", &(decayMode_continu));
  continuTree->SetBranchAddress("MBC", &(MBC_continu));
  continuTree->SetBranchAddress("DeltaM", &(DeltaM_continu));
  continuTree->SetBranchAddress("e9oe25Unf", &(e9oe25Unf_continu));
  
  physicsTree->SetBranchAddress("Run", &(runNumber_physics));
  physicsTree->SetBranchAddress("Event", &(eventNumber_physics));
  physicsTree->SetBranchAddress("dsPlusM", &(dsPlusM_physics));
  physicsTree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_physics));
  physicsTree->SetBranchAddress("DecayMode", &(decayMode_physics));
  physicsTree->SetBranchAddress("MBC", &(MBC_physics));
  physicsTree->SetBranchAddress("DeltaM", &(DeltaM_physics));
  physicsTree->SetBranchAddress("e9oe25Unf", &(e9oe25Unf_physics));
  
  wrongConverTree_Dsp->SetBranchAddress("Run", &(runNumber_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("Event", &(eventNumber_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("Ebeam", &(beamEnergy_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("dsPlusM", &(dsPlusM_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("DecayMode", &(decayMode_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("MBC", &(MBC_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("DeltaM", &(DeltaM_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("e9oe25Unf", &(e9oe25Unf_wrongConver_Dsp));
  // Need these for MC truth matching
  wrongConverTree_Dsp->SetBranchAddress("kDsPlusPx_reco", &(dsPlusPx_reco_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("kDsPlusPy_reco", &(dsPlusPy_reco_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("kDsPlusPz_reco", &(dsPlusPz_reco_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("kDsPlusPx_MC", &(dsPlusPx_MC_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("kDsPlusPy_MC", &(dsPlusPy_MC_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("kDsPlusPz_MC", &(dsPlusPz_MC_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("kPhotonPx_reco", &(photonPx_reco_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("kPhotonPy_reco", &(photonPy_reco_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("kPhotonPz_reco", &(photonPz_reco_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("kPhotonE_reco", &(photonE_reco_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("kPhotonPx_MC", &(photonPx_MC_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("kPhotonPy_MC", &(photonPy_MC_wrongConver_Dsp));
  wrongConverTree_Dsp->SetBranchAddress("kPhotonPz_MC", &(photonPz_MC_wrongConver_Dsp));
  
  wrongConverTree_Dsm->SetBranchAddress("Run", &(runNumber_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("Event", &(eventNumber_wrongConver_Dsm));
  wrongConverTree_Dsp->SetBranchAddress("Ebeam", &(beamEnergy_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("dsPlusM", &(dsPlusM_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("DecayMode", &(decayMode_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("MBC", &(MBC_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("DeltaM", &(DeltaM_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("e9oe25Unf", &(e9oe25Unf_wrongConver_Dsm));
  // Need these for MC truth matching
  wrongConverTree_Dsm->SetBranchAddress("kDsPlusPx_reco", &(dsPlusPx_reco_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("kDsPlusPy_reco", &(dsPlusPy_reco_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("kDsPlusPz_reco", &(dsPlusPz_reco_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("kDsPlusPx_MC", &(dsPlusPx_MC_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("kDsPlusPy_MC", &(dsPlusPy_MC_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("kDsPlusPz_MC", &(dsPlusPz_MC_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("kPhotonPx_reco", &(photonPx_reco_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("kPhotonPy_reco", &(photonPy_reco_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("kPhotonPz_reco", &(photonPz_reco_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("kPhotonE_reco", &(photonE_reco_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("kPhotonPx_MC", &(photonPx_MC_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("kPhotonPy_MC", &(photonPy_MC_wrongConver_Dsm));
  wrongConverTree_Dsm->SetBranchAddress("kPhotonPz_MC", &(photonPz_MC_wrongConver_Dsm));
  
  int nConverEvents_Dsp=converTree_Dsp->GetEntries();
  std::cout<<"Processing Conver_Dsp MC Sample of "<<nConverEvents_Dsp<<" candidates ... "; std::cout.flush();
  int startTime=time(0);
  for (int i=0; i<nConverEvents_Dsp; ++i)
  // for (int i=0; i<1; ++i)
  {
    converTree_Dsp->GetEvent(i);
    
    if (runNumber_conver_Dsp!=converNumber_Dsp.noCut_run || eventNumber_conver_Dsp!=converNumber_Dsp.noCut_event)
    {
      ++converEvents_Dsp.noCut;
      converNumber_Dsp.noCut_run=runNumber_conver_Dsp;
      converNumber_Dsp.noCut_event=eventNumber_conver_Dsp;
    }
    
    if (decayMode_conver_Dsp==decayNumber && dsPlusCharge_conver_Dsp==1 && x925Cut(e9oe25Unf_conver_Dsp))
    {
      h_dsPlusM_conver_Dsp->Fill(dsPlusM_conver_Dsp);
      if (dsPlusMCut(dsPlusM_conver_Dsp))
      {
        if (!mBC_sideband) h_MBC_conver_Dsp->Fill(MBC_conver_Dsp);
        if (mBC_sideband || MBCCut(MBC_conver_Dsp))
        {
          if (!DeltaM_sideband) h_DeltaM_conver_Dsp->Fill(DeltaM_conver_Dsp);
          if (DeltaM_sideband || DeltaMCut(DeltaM_conver_Dsp))
          {
            if (runNumber_conver_Dsp!=converNumber_Dsp.deltaMCut_run || eventNumber_conver_Dsp!=converNumber_Dsp.deltaMCut_event)
            {
              ++converEvents_Dsp.deltaMCut;
              converNumber_Dsp.deltaMCut_run=runNumber_conver_Dsp;
              converNumber_Dsp.deltaMCut_event=eventNumber_conver_Dsp;
            }
            if (mBC_sideband)
            {
              h_MBC_conver_Dsp->Fill(MBC_conver_Dsp);
              h_photonE_MBC_conver->Fill(MBC_conver_Dsp, photonE_conver_Dsp-photonE_conver_Dsp_MC);
              if (MBC_conver_Dsp>lowerSide_xmin && MBC_conver_Dsp<=lowerSide_xmax) ++conver_Dsp_beforeSignal;
              if (MBC_conver_Dsp>xmin && MBC_conver_Dsp<xmax) ++conver_Dsp_atSignal;
              if (MBC_conver_Dsp>=upperSide_xmin && MBC_conver_Dsp<upperSide_xmax) ++conver_Dsp_afterSignal;
              TVector3 dsPlus_reco_conver_Dsp(dsPlusPx_reco_conver_Dsp, dsPlusPy_reco_conver_Dsp, dsPlusPz_reco_conver_Dsp);
              TVector3 dsPlus_MC_conver_Dsp(dsPlusPx_MC_conver_Dsp, dsPlusPy_MC_conver_Dsp, dsPlusPz_MC_conver_Dsp);
              TVector3 photon_reco_conver_Dsp(photonPx_reco_conver_Dsp, photonPy_reco_conver_Dsp, photonPz_reco_conver_Dsp);
              TVector3 photon_MC_conver_Dsp(photonPx_MC_conver_Dsp, photonPy_MC_conver_Dsp, photonPz_MC_conver_Dsp);
              double dsPlus_angle=dsPlus_reco_conver_Dsp.Angle(dsPlus_MC_conver_Dsp);
              double photon_angle=photon_reco_conver_Dsp.Angle(photon_MC_conver_Dsp);
              if (dsPlus_angle<DsMatchingAngle)
              {
                 if (photon_angle<photonMatchingAngle) h_MBC_conver_Dsp_matched->Fill(MBC_conver_Dsp);
                 else h_MBC_conver_Dsp_unmatched->Fill(MBC_conver_Dsp);
              }
            }
            if (DeltaM_sideband)
            {
              h_DeltaM_conver_Dsp->Fill(DeltaM_conver_Dsp);
              if (DeltaM_conver_Dsp>lowerSide_xmin && DeltaM_conver_Dsp<=lowerSide_xmax) ++conver_Dsp_beforeSignal;
              if (DeltaM_conver_Dsp>xmin && DeltaM_conver_Dsp<xmax) ++conver_Dsp_atSignal;
              if (DeltaM_conver_Dsp>=upperSide_xmin && DeltaM_conver_Dsp<upperSide_xmax) ++conver_Dsp_afterSignal;
            }
            
            h_e9oe25Unf_conver_Dsp->Fill(e9oe25Unf_conver_Dsp);
            
          }
        }
      }
    }
    if ((time(0)-startTime)>1)
    {
      std::cout<<(float)i*100/(float)nConverEvents_Dsp<<"%, "; std::cout.flush();
      startTime=time(0);
    }
  }
  
  std::cout<<"done."<<std::endl;
  
  int nConverEvents_Dsm=converTree_Dsm->GetEntries();
  std::cout<<"Processing Conver_Dsm MC Sample of "<<nConverEvents_Dsm<<" candidates ... "; std::cout.flush();
  startTime=time(0);
  for (int i=0; i<nConverEvents_Dsm; ++i)
  // for (int i=0; i<1; ++i)
  {
    converTree_Dsm->GetEvent(i);
    
    if (runNumber_conver_Dsm!=converNumber_Dsm.noCut_run || eventNumber_conver_Dsm!=converNumber_Dsm.noCut_event)
    {
      ++converEvents_Dsm.noCut;
      converNumber_Dsm.noCut_run=runNumber_conver_Dsm;
      converNumber_Dsm.noCut_event=eventNumber_conver_Dsm;
    }
    
    if (decayMode_conver_Dsm==decayNumber && dsPlusCharge_conver_Dsm==-1 && x925Cut(e9oe25Unf_conver_Dsm))
    {
      h_dsPlusM_conver_Dsm->Fill(dsPlusM_conver_Dsm);
      if (dsPlusMCut(dsPlusM_conver_Dsm))
      {
        if (!mBC_sideband) h_MBC_conver_Dsm->Fill(MBC_conver_Dsm);
        if (mBC_sideband || MBCCut(MBC_conver_Dsm))
        {
          if (!DeltaM_sideband) h_DeltaM_conver_Dsm->Fill(DeltaM_conver_Dsm);
          if (DeltaM_sideband || DeltaMCut(DeltaM_conver_Dsm))
          {
            if (runNumber_conver_Dsm!=converNumber_Dsm.deltaMCut_run || eventNumber_conver_Dsm!=converNumber_Dsm.deltaMCut_event)
            {
              ++converEvents_Dsm.deltaMCut;
              converNumber_Dsm.deltaMCut_run=runNumber_conver_Dsm;
              converNumber_Dsm.deltaMCut_event=eventNumber_conver_Dsm;
            }
            if (mBC_sideband)
            {
              h_MBC_conver_Dsm->Fill(MBC_conver_Dsm);
              h_photonE_MBC_conver->Fill(MBC_conver_Dsm, photonE_conver_Dsm-photonE_conver_Dsm_MC);
              if (MBC_conver_Dsm>lowerSide_xmin && MBC_conver_Dsm<=lowerSide_xmax) ++conver_Dsm_beforeSignal;
              if (MBC_conver_Dsm>xmin && MBC_conver_Dsm<xmax) ++conver_Dsm_atSignal;
              if (MBC_conver_Dsm>=upperSide_xmin && MBC_conver_Dsm<upperSide_xmax) ++conver_Dsm_afterSignal;
              TVector3 dsPlus_reco_conver_Dsm(dsPlusPx_reco_conver_Dsm, dsPlusPy_reco_conver_Dsm, dsPlusPz_reco_conver_Dsm);
              TVector3 dsPlus_MC_conver_Dsm(dsPlusPx_MC_conver_Dsm, dsPlusPy_MC_conver_Dsm, dsPlusPz_MC_conver_Dsm);
              TVector3 photon_reco_conver_Dsm(photonPx_reco_conver_Dsm, photonPy_reco_conver_Dsm, photonPz_reco_conver_Dsm);
              TVector3 photon_MC_conver_Dsm(photonPx_MC_conver_Dsm, photonPy_MC_conver_Dsm, photonPz_MC_conver_Dsm);
              double dsPlus_angle=dsPlus_reco_conver_Dsm.Angle(dsPlus_MC_conver_Dsm);
              double photon_angle=photon_reco_conver_Dsm.Angle(photon_MC_conver_Dsm);
              if (dsPlus_angle<DsMatchingAngle)
              {
                 if (photon_angle<photonMatchingAngle) h_MBC_conver_Dsm_matched->Fill(MBC_conver_Dsm);
                 else h_MBC_conver_Dsm_unmatched->Fill(MBC_conver_Dsm);
              }
            }
            if (DeltaM_sideband)
            {
              h_DeltaM_conver_Dsm->Fill(DeltaM_conver_Dsm);
              if (DeltaM_conver_Dsm>lowerSide_xmin && DeltaM_conver_Dsm<=lowerSide_xmax) ++conver_Dsm_beforeSignal;
              if (DeltaM_conver_Dsm>xmin && DeltaM_conver_Dsm<xmax) ++conver_Dsm_atSignal;
              if (DeltaM_conver_Dsm>=upperSide_xmin && DeltaM_conver_Dsm<upperSide_xmax) ++conver_Dsm_afterSignal;
            }
            
            h_e9oe25Unf_conver_Dsm->Fill(e9oe25Unf_conver_Dsm);
            
          }
        }
      }
    }
    if ((time(0)-startTime)>1)
    {
      std::cout<<(float)i*100/(float)nConverEvents_Dsm<<"%, "; std::cout.flush();
      startTime=time(0);
    }
  }
  
  std::cout<<"done."<<std::endl;
  
  int nGenericEvents=genericTree->GetEntries();
  std::cout<<"Processing Generic MC Sample of "<<nGenericEvents<<" candidates ... "; std::cout.flush();
  startTime=time(0);
  for (int i=0; i<nGenericEvents; ++i)
  // for (int i=0; i<1; ++i)
  {
    genericTree->GetEvent(i);
    
    if (runNumber_generic!=genericNumber.noCut_run || eventNumber_generic!=genericNumber.noCut_event)
    {
      ++genericEvents.noCut;
      genericNumber.noCut_run=runNumber_generic;
      genericNumber.noCut_event=eventNumber_generic;
    }
    
    if (decayMode_generic==decayNumber && x925Cut(e9oe25Unf_generic))
    {
      h_dsPlusM_generic->Fill(dsPlusM_generic);
      if (dsPlusMCut(dsPlusM_generic))
      {
        if (!mBC_sideband) h_MBC_generic->Fill(MBC_generic);
        if (mBC_sideband || MBCCut(MBC_generic))
        {
          if (!DeltaM_sideband) h_DeltaM_generic->Fill(DeltaM_generic);
          if (DeltaM_sideband || DeltaMCut(DeltaM_generic))
          {
            if (runNumber_generic!=genericNumber.deltaMCut_run || eventNumber_generic!=genericNumber.deltaMCut_event)
            {
              ++genericEvents.deltaMCut;
              genericNumber.deltaMCut_run=runNumber_generic;
              genericNumber.deltaMCut_event=eventNumber_generic;
            }
            if (mBC_sideband)
            {
              if (conversionBit_generic*dsPlusCharge_generic==1) h_MBC_generic->Fill(MBC_generic);
              if (conversionBit_generic==0) h_MBC_generic_veto->Fill(MBC_generic);
              if (conversionBit_generic*dsPlusCharge_generic==-1) h_MBC_generic_wrong->Fill(MBC_generic);
              if (MBC_generic>lowerSide_xmin && MBC_generic<=lowerSide_xmax) ++generic_beforeSignal;
              if (MBC_generic>xmin && MBC_generic<xmax) ++generic_atSignal;
              if (MBC_generic>=upperSide_xmin && MBC_generic<upperSide_xmax) ++generic_afterSignal;
            }
            if (DeltaM_sideband)
            {
              h_DeltaM_generic->Fill(DeltaM_generic);
              if (conversionBit_generic==0) h_DeltaM_generic_veto->Fill(DeltaM_generic);
              if (conversionBit_generic*dsPlusCharge_generic==-1) h_DeltaM_generic_wrong->Fill(DeltaM_generic);
              if (DeltaM_generic>lowerSide_xmin && DeltaM_generic<=lowerSide_xmax) ++generic_beforeSignal;
              if (DeltaM_generic>xmin && DeltaM_generic<xmax) ++generic_atSignal;
              if (DeltaM_generic>=upperSide_xmin && DeltaM_generic<upperSide_xmax) ++generic_afterSignal;
            }
            
            h_e9oe25Unf_generic->Fill(e9oe25Unf_generic);
            
          }
        }
      }
    }
    if ((time(0)-startTime)>10)
    {
      std::cout<<(float)i*100/(float)nGenericEvents<<"%, "; std::cout.flush();
      startTime=time(0);
    }
  }
  
  std::cout<<"done."<<std::endl;
  
  int nContinuEvents=continuTree->GetEntries();
  std::cout<<"Processing Continuum MC Sample... "; std::cout.flush();
  startTime=time(0);
  for (int i=0; i<nContinuEvents; ++i)
  // for (int i=0; i<1; ++i)
  {
    continuTree->GetEvent(i);
    
    if (runNumber_continu!=continuNumber.noCut_run || eventNumber_continu!=continuNumber.noCut_event)
    {
      ++continuEvents.noCut;
      continuNumber.noCut_run=runNumber_continu;
      continuNumber.noCut_event=eventNumber_continu;
    }
    
    if (decayMode_continu==decayNumber && x925Cut(e9oe25Unf_continu))
    {
      h_dsPlusM_continu->Fill(dsPlusM_continu);
      if (dsPlusMCut(dsPlusM_continu))
      {
        if (!mBC_sideband) h_MBC_continu->Fill(MBC_continu);
        if (mBC_sideband || MBCCut(MBC_continu))
        {
          if (!DeltaM_sideband) h_DeltaM_continu->Fill(DeltaM_continu);
          if (DeltaM_sideband || DeltaMCut(DeltaM_continu))
          {
            if (runNumber_continu!=continuNumber.deltaMCut_run || eventNumber_continu!=continuNumber.deltaMCut_event)
            {
              ++continuEvents.deltaMCut;
              continuNumber.deltaMCut_run=runNumber_continu;
              continuNumber.deltaMCut_event=eventNumber_continu;
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
            
            h_e9oe25Unf_continu->Fill(e9oe25Unf_continu);
            
          }
        }
      }
    }
    if ((time(0)-startTime)>10)
    {
      std::cout<<(float)i*100/(float)nContinuEvents<<"%, "; std::cout.flush();
      startTime=time(0);
    }
  }
  
  std::cout<<"done."<<std::endl;
  
  int nPhysicsEvents=physicsTree->GetEntries();
  std::cout<<"Processing Physics Data... "; std::cout.flush();
  startTime=time(0);
  for (int i=0; i<nPhysicsEvents; ++i)
  // for (int i=0; i<1; ++i)
  {
    physicsTree->GetEvent(i);
    
    if (runNumber_physics!=physicsNumber.noCut_run || eventNumber_physics!=physicsNumber.noCut_event)
    {
      ++physicsEvents.noCut;
      physicsNumber.noCut_run=runNumber_physics;
      physicsNumber.noCut_event=eventNumber_physics;
    }
    
    if (decayMode_physics==decayNumber && x925Cut(e9oe25Unf_physics))
    {
      h_dsPlusM_physics->Fill(dsPlusM_physics);
      if (dsPlusMCut(dsPlusM_physics))
      {
        if (!mBC_sideband) h_MBC_physics->Fill(MBC_physics);
        if (mBC_sideband || MBCCut(MBC_physics))
        {
          if (!DeltaM_sideband) h_DeltaM_physics->Fill(DeltaM_physics);
          if (DeltaM_sideband || DeltaMCut(DeltaM_physics))
          {
            if (runNumber_physics!=physicsNumber.deltaMCut_run || eventNumber_physics!=physicsNumber.deltaMCut_event)
            {
              ++physicsEvents.deltaMCut;
              physicsNumber.deltaMCut_run=runNumber_physics;
              physicsNumber.deltaMCut_event=eventNumber_physics;
            }
            if (mBC_sideband)
            {
              h_MBC_physics->Fill(MBC_physics);
              if (MBC_physics>lowerSide_xmin && MBC_physics<=lowerSide_xmax) ++physics_beforeSignal;
              if (MBC_physics>xmin && MBC_physics<xmax) ++physics_atSignal;
              if (MBC_physics>=upperSide_xmin && MBC_physics<upperSide_xmax) ++physics_afterSignal;
            }
            if (DeltaM_sideband)
            {
              h_DeltaM_physics->Fill(DeltaM_physics);
              if (DeltaM_physics>lowerSide_xmin && DeltaM_physics<=lowerSide_xmax) ++physics_beforeSignal;
              if (DeltaM_physics>xmin && DeltaM_physics<xmax) ++physics_atSignal;
              if (DeltaM_physics>=upperSide_xmin && DeltaM_physics<upperSide_xmax) ++physics_afterSignal;
            }
            
            h_e9oe25Unf_physics->Fill(e9oe25Unf_physics);
            
          }
        }
      }
    }
    if ((time(0)-startTime)>10)
    {
      std::cout<<(float)i*100/(float)nPhysicsEvents<<"%, "; std::cout.flush();
      startTime=time(0);
    }
  }
  
  std::cout<<"done."<<std::endl;
  
  int nWrongConverEvents_Dsp=wrongConverTree_Dsp->GetEntries();
  std::cout<<"Processing Wrong Conver_Dsp MC Sample of "<<nWrongConverEvents_Dsp<<" candidates ... "; std::cout.flush();
  startTime=time(0);
  for (int i=0; i<nWrongConverEvents_Dsp; ++i)
  // for (int i=0; i<1; ++i)
  {
    wrongConverTree_Dsp->GetEvent(i);
    
    if (decayMode_wrongConver_Dsp==decayNumber && dsPlusCharge_wrongConver_Dsp==-1 && x925Cut(e9oe25Unf_wrongConver_Dsp))
    {
      if (dsPlusMCut(dsPlusM_wrongConver_Dsp))
      {
        if (mBC_sideband || MBCCut(MBC_wrongConver_Dsp))
        {
          if (DeltaM_sideband || DeltaMCut(DeltaM_wrongConver_Dsp))
          {
            if (mBC_sideband) 
            {
              h_MBC_wrongConver_Dsp->Fill(MBC_wrongConver_Dsp);
              TVector3 photon_reco_wrongConver_Dsp(photonPx_reco_wrongConver_Dsp, photonPy_reco_wrongConver_Dsp, photonPz_reco_wrongConver_Dsp);
              TVector3 photon_MC_wrongConver_Dsp(photonPx_MC_wrongConver_Dsp, photonPy_MC_wrongConver_Dsp, photonPz_MC_wrongConver_Dsp);
              double photon_angle=photon_reco_wrongConver_Dsp.Angle(photon_MC_wrongConver_Dsp);
              TVector3 Dsm_reco(dsPlusPx_reco_wrongConver_Dsp, dsPlusPy_reco_wrongConver_Dsp, dsPlusPz_reco_wrongConver_Dsp);
              TVector3 Dsp_MC(dsPlusPx_MC_wrongConver_Dsp, dsPlusPy_MC_wrongConver_Dsp, dsPlusPz_MC_wrongConver_Dsp);
              TVector3 Dsm_MC=-(Dsp_MC+photon_MC_wrongConver_Dsp);
              double Ds_angle=Dsm_reco.Angle(Dsm_MC);
              h_angle_photon->Fill(photon_angle);
              h_angle_wrongDs->Fill(Ds_angle);
              if (Ds_angle<DsMatchingAngle)
              {
                if (photon_angle<photonMatchingAngle)
                {
                  h_MBC_wrongConver_Dsp_matched->Fill(MBC_wrongConver_Dsp);
                  h_photonE_wrongConver_Dsp_matched->Fill(photonE_reco_wrongConver_Dsp);
                }
                else
                {
                  h_MBC_wrongConver_Dsp_unmatched->Fill(MBC_wrongConver_Dsp);
                  h_photonE_wrongConver_Dsp_unmatched->Fill(photonE_reco_wrongConver_Dsp);
                }
              }
              else
              {
                if (photon_angle>photonMatchingAngle)
                {
                  h_MBC_wrongConver_Dsp_strictUnmatched->Fill(MBC_wrongConver_Dsp);
                }
              }
            }
            if (DeltaM_sideband) h_DeltaM_wrongConver_Dsp->Fill(DeltaM_wrongConver_Dsp);
          }
        }
      }
    }
    if ((time(0)-startTime)>10)
    {
      std::cout<<(float)i*100/(float)nWrongConverEvents_Dsp<<"%, "; std::cout.flush();
      startTime=time(0);
    }
  }
  
  std::cout<<"done."<<std::endl;
  
  int nWrongConverEvents_Dsm=wrongConverTree_Dsm->GetEntries();
  std::cout<<"Processing Wrong Conver_Dsm MC Sample of "<<nWrongConverEvents_Dsm<<" candidates ... "; std::cout.flush();
  startTime=time(0);
  for (int i=0; i<nWrongConverEvents_Dsm; ++i)
  // for (int i=0; i<1; ++i)
  {
    wrongConverTree_Dsm->GetEvent(i);
    
    if (decayMode_wrongConver_Dsm==decayNumber && dsPlusCharge_wrongConver_Dsm==+1 && x925Cut(e9oe25Unf_wrongConver_Dsm))
    {
      if (dsPlusMCut(dsPlusM_wrongConver_Dsm))
      {
        if (mBC_sideband || MBCCut(MBC_wrongConver_Dsm))
        {
          if (DeltaM_sideband || DeltaMCut(DeltaM_wrongConver_Dsm))
          {
            if (mBC_sideband) 
            {
              h_MBC_wrongConver_Dsm->Fill(MBC_wrongConver_Dsm);
              TVector3 photon_reco_wrongConver_Dsm(photonPx_reco_wrongConver_Dsm, photonPy_reco_wrongConver_Dsm, photonPz_reco_wrongConver_Dsm);
              TVector3 photon_MC_wrongConver_Dsm(photonPx_MC_wrongConver_Dsm, photonPy_MC_wrongConver_Dsm, photonPz_MC_wrongConver_Dsm);
              double photon_angle=photon_reco_wrongConver_Dsm.Angle(photon_MC_wrongConver_Dsm);
              TVector3 Dsp_reco(dsPlusPx_reco_wrongConver_Dsm, dsPlusPy_reco_wrongConver_Dsm, dsPlusPz_reco_wrongConver_Dsm);
              TVector3 Dsm_MC(dsPlusPx_MC_wrongConver_Dsm, dsPlusPy_MC_wrongConver_Dsm, dsPlusPz_MC_wrongConver_Dsm);
              TVector3 Dsp_MC=-(Dsm_MC+photon_MC_wrongConver_Dsm);
              double Ds_angle=Dsp_reco.Angle(Dsp_MC);
              if (Ds_angle<DsMatchingAngle)
              {
                if (photon_angle<photonMatchingAngle)
                {
                  h_MBC_wrongConver_Dsm_matched->Fill(MBC_wrongConver_Dsm);
                  h_photonE_wrongConver_Dsm_matched->Fill(photonE_reco_wrongConver_Dsm);
                }
                else
                {
                  h_MBC_wrongConver_Dsm_unmatched->Fill(MBC_wrongConver_Dsm);
                  h_photonE_wrongConver_Dsm_unmatched->Fill(photonE_reco_wrongConver_Dsm);
                }
              }
              else
              {
                if (photon_angle>photonMatchingAngle)
                {
                  h_MBC_wrongConver_Dsm_strictUnmatched->Fill(MBC_wrongConver_Dsm);
                }
              }
            }
            if (DeltaM_sideband) h_DeltaM_wrongConver_Dsm->Fill(DeltaM_wrongConver_Dsm);
          }
        }
      }
    }
    if ((time(0)-startTime)>10)
    {
      std::cout<<(float)i*100/(float)nWrongConverEvents_Dsm<<"%, "; std::cout.flush();
      startTime=time(0);
    }
  }
  
  std::cout<<"done."<<std::endl;
  
  // Report the efficiencies of the cuts on the conver samples
  double eff_Dsp=conver_Dsp_atSignal/nConversionSample_Dsp;
  double eff_Dsm=conver_Dsm_atSignal/nConversionSample_Dsm;
  double eff_Dsp_error=eff_Dsp/sqrt(conver_Dsp_atSignal);
  double eff_Dsm_error=eff_Dsm/sqrt(conver_Dsm_atSignal);
  std::cout<<"--- Cut Efficiencies --- "<<std::endl;
  std::cout<<" We started with "<<nConversionSample_Dsp<<" Ds+ gamma events and ended up with "<<conver_Dsp_atSignal
           <<", giving us a Ds+ efficiency = "<<eff_Dsp<<"+-"<<eff_Dsp_error<<std::endl;
  std::cout<<" We started with "<<nConversionSample_Dsm<<" Ds- gamma events and ended up with "<<conver_Dsm_atSignal
           <<", giving us a Ds- efficiency = "<<eff_Dsm<<"+-"<<eff_Dsm_error<<std::endl;
  std::cout<<" Average = "<<(eff_Dsp+eff_Dsm)/2<<"+-"<<(eff_Dsp_error+eff_Dsm_error)/2<<std::endl;
  std::cout<<"---"<<std::endl;
  
  // Scale plots appropriately
  h_MBC_conver_Dsp->Scale(1./(2.*nConversionSample_Dsp));
  h_MBC_conver_Dsm->Scale(1./(2.*nConversionSample_Dsm));
  h_MBC_generic->Scale(genericScale);
  h_MBC_generic_veto->Scale(genericScale);
  h_MBC_generic_wrong->Scale(genericScale);
  h_MBC_continu->Scale(continuScale);
  h_MBC_wrongConver_Dsp->Scale(1./(2.*nWrongConversionSample_Dsp));
  h_MBC_wrongConver_Dsm->Scale(1./(2.*nWrongConversionSample_Dsm));
  h_DeltaM_conver_Dsp->Scale(1./(2.*nConversionSample_Dsp));
  h_DeltaM_conver_Dsm->Scale(1./(2.*nConversionSample_Dsp));
  h_DeltaM_generic->Scale(genericScale);
  h_DeltaM_generic_veto->Scale(genericScale);
  h_DeltaM_generic_wrong->Scale(genericScale);
  h_DeltaM_continu->Scale(continuScale);
  h_DeltaM_wrongConver_Dsp->Scale(1./(2.*nWrongConversionSample_Dsp));
  h_DeltaM_wrongConver_Dsm->Scale(1./(2.*nWrongConversionSample_Dsm));
  h_MBC_conver_Dsp_matched->Scale(1./(2.*nConversionSample_Dsp));
  h_MBC_conver_Dsm_matched->Scale(1./(2.*nConversionSample_Dsm));
  h_MBC_conver_Dsp_unmatched->Scale(1./(2.*nConversionSample_Dsp));
  h_MBC_conver_Dsm_unmatched->Scale(1./(2.*nConversionSample_Dsm));
  h_MBC_wrongConver_Dsp_matched->Scale(1./(2.*nWrongConversionSample_Dsp));
  h_MBC_wrongConver_Dsm_matched->Scale(1./(2.*nWrongConversionSample_Dsm));
  h_MBC_wrongConver_Dsp_unmatched->Scale(1./(2.*nWrongConversionSample_Dsp));
  h_MBC_wrongConver_Dsm_unmatched->Scale(1./(2.*nWrongConversionSample_Dsm));
  h_MBC_wrongConver_Dsp_strictUnmatched->Scale(1./(2.*nWrongConversionSample_Dsp));
  h_MBC_wrongConver_Dsm_strictUnmatched->Scale(1./(2.*nWrongConversionSample_Dsm));
  h_photonE_wrongConver_Dsp_unmatched->Scale(1./(2.*nWrongConversionSample_Dsp));
  h_photonE_wrongConver_Dsm_unmatched->Scale(1./(2.*nWrongConversionSample_Dsm));
  h_photonE_wrongConver_Dsp_matched->Scale(1./(2.*nWrongConversionSample_Dsp));
  h_photonE_wrongConver_Dsm_matched->Scale(1./(2.*nWrongConversionSample_Dsm));
  
  // And scale numbers too
  generic_beforeSignal*=genericScale; generic_afterSignal*=genericScale; generic_atSignal*=genericScale;
  continu_beforeSignal*=continuScale; continu_afterSignal*=continuScale; continu_atSignal*=continuScale;
  
  // Add conversion plots
  h_MBC_conver->Add(h_MBC_conver_Dsp, h_MBC_conver_Dsm);
  h_MBC_conver_matched->Add(h_MBC_conver_Dsp_matched, h_MBC_conver_Dsm_matched);
  h_MBC_conver_unmatched->Add(h_MBC_conver_Dsp_unmatched, h_MBC_conver_Dsm_unmatched);
  h_MBC_wrongConver_matched->Add(h_MBC_wrongConver_Dsp_matched, h_MBC_wrongConver_Dsm_matched);
  h_MBC_wrongConver_unmatched->Add(h_MBC_wrongConver_Dsp_unmatched, h_MBC_wrongConver_Dsm_unmatched);
  h_MBC_wrongConver_strictUnmatched->Add(h_MBC_wrongConver_Dsp_strictUnmatched, h_MBC_wrongConver_Dsm_strictUnmatched);
  h_DeltaM_conver->Add(h_DeltaM_conver_Dsp, h_DeltaM_conver_Dsm);
  h_MBC_wrongConver->Add(h_MBC_wrongConver_Dsp, h_MBC_wrongConver_Dsm);
  h_DeltaM_wrongConver->Add(h_DeltaM_wrongConver_Dsp, h_DeltaM_wrongConver_Dsm);
  h_photonE_wrongConver_unmatched->Add(h_photonE_wrongConver_Dsp_unmatched, h_photonE_wrongConver_Dsm_unmatched);
  h_photonE_wrongConver_matched->Add(h_photonE_wrongConver_Dsp_matched, h_photonE_wrongConver_Dsm_matched);
  
  // Now calculate errors
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
  
  gROOT->SetStyle("Plain");
  double ymin, ymax;
  TLine *line;
  
  if (mBC_sideband)
  {
    /*
    h_MBC_generic->SetFillColor(kGreen);
    h_MBC_continu->SetFillColor(kBlue);
    
    TCanvas *c_MBC_Stacked = new TCanvas("c_MBC_Stacked");
    THStack *s_MBC_sideband=new THStack("s_MBC_sideband", ("m_{BC} Sidebands in "+decay).c_str());
    s_MBC_sideband->Add(h_MBC_continu, "hist");
    s_MBC_sideband->Add(h_MBC_generic, "hist");
    double stack_ymax=h_MBC_generic->GetMaximum();
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
    std::string continu_string="Continuum MC: "; continu_string+=itoa(continuScale*h_MBC_continu->GetEntries()); continu_string+=" Events";
    std::string generic_string="Generic MC: "; generic_string+=itoa(genericScale*h_MBC_generic->GetEntries()); generic_string+=" Events";
    std::string physics_string="Data: "; physics_string+=itoa(h_MBC_physics->GetEntries()); physics_string+=" Events";
    legendMBC->AddEntry(h_MBC_continu, continu_string.c_str());
    legendMBC->AddEntry(h_MBC_generic, generic_string.c_str());
    legendMBC->AddEntry(h_MBC_physics, physics_string.c_str());
    legendMBC->SetFillColor(kWhite);
    legendMBC->Draw();
    
    // std::string filename_eps=decay+"_Sideband_MBC_002.eps";
    // c_MBC_Stacked->SaveAs(filename_eps.c_str());
    TCanvas *c_MBC_conver = new TCanvas("c_MBC_conver");
    c_MBC_conver->Divide(1,3);
    c_MBC_conver->cd(1);
    h_MBC_conver_Dsp->Draw();
    c_MBC_conver->cd(2);
    h_MBC_conver_Dsm->Draw();
    c_MBC_conver->cd(3);
    h_MBC_conver->Draw();
    
    TCanvas *c_MBC_conver_matched = new TCanvas("c_MBC_conver_matched");
    c_MBC_conver_matched->Divide(1,3);
    c_MBC_conver_matched->cd(1);
    h_MBC_conver_Dsp_matched->Draw();
    c_MBC_conver_matched->cd(2);
    h_MBC_conver_Dsm_matched->Draw();
    c_MBC_conver_matched->cd(3);
    h_MBC_conver_matched->Draw();
    
    TCanvas *c_MBC_conver_unmatched=new TCanvas("c_MBC_conver_unmatched");
    c_MBC_conver_unmatched->Divide(1,3);
    c_MBC_conver_unmatched->cd(1);
    h_MBC_conver_Dsp_unmatched->Draw();
    c_MBC_conver_unmatched->cd(2);
    h_MBC_conver_Dsm_unmatched->Draw();
    c_MBC_conver_unmatched->cd(3);
    h_MBC_conver_unmatched->Draw();
    
    TCanvas *c_MBC_wrongConver = new TCanvas("c_MBC_wrongConver");
    c_MBC_wrongConver->Divide(1,3);
    c_MBC_wrongConver->cd(1);
    h_MBC_wrongConver_Dsp->Draw();
    c_MBC_wrongConver->cd(2);
    h_MBC_wrongConver_Dsm->Draw();
    c_MBC_wrongConver->cd(3);
    h_MBC_wrongConver->Draw();
    
    TCanvas *c_MBC_wrongConver_matched = new TCanvas("c_MBC_wrongConver_matched");
    c_MBC_wrongConver_matched->Divide(1,3);
    c_MBC_wrongConver_matched->cd(1);
    h_MBC_wrongConver_Dsp_matched->Draw();
    c_MBC_wrongConver_matched->cd(2);
    h_MBC_wrongConver_Dsm_matched->Draw();
    c_MBC_wrongConver_matched->cd(3);
    h_MBC_wrongConver_matched->Draw();
    
    TCanvas *c_MBC_wrongConver_unmatched=new TCanvas("c_MBC_wrongConver_unmatched");
    c_MBC_wrongConver_unmatched->Divide(1,3);
    c_MBC_wrongConver_unmatched->cd(1);
    h_MBC_wrongConver_Dsp_unmatched->Draw();
    c_MBC_wrongConver_unmatched->cd(2);
    h_MBC_wrongConver_Dsm_unmatched->Draw();
    c_MBC_wrongConver_unmatched->cd(3);
    h_MBC_wrongConver_unmatched->Draw();
    
    TCanvas *c_MBC_wrongConver_strictUnmatched = new TCanvas("c_MBC_wrongConver_strictUnmatched");
    h_MBC_wrongConver_strictUnmatched->Draw();
    
    TCanvas *c_angle_photon = new TCanvas("c_angle_photon", "c_angle_photon");
    h_angle_photon->Draw();
    ymax=0.75*h_angle_photon->GetMaximum();
    line=new TLine(photonMatchingAngle, ymax, photonMatchingAngle, 0); line->Draw();
    
    TCanvas *c_angle_wrongDs = new TCanvas("c_angle_wrongDs", "c_angle_wrongDs");
    h_angle_wrongDs->Draw();
    ymax=0.75*h_angle_wrongDs->GetMaximum();
    line=new TLine(DsMatchingAngle, ymax, DsMatchingAngle, 0); line->Draw();
    
    TCanvas *c_photonE_MBC = new TCanvas("c_photonE_MBC");
    h_photonE_MBC_conver->Draw("box");
    
    TCanvas *c_photonE_wrongSign = new TCanvas("c_photonE_wrongSign");
    c_photonE_wrongSign->Divide(1,3);
    c_photonE_wrongSign->cd(1);
    h_photonE_wrongConver_Dsp_unmatched->Draw();
    c_photonE_wrongSign->cd(2);
    h_photonE_wrongConver_Dsm_unmatched->Draw();
    c_photonE_wrongSign->cd(3);
    h_photonE_wrongConver_unmatched->Draw();
    
    TCanvas *c_photonE_wrongSign_matched = new TCanvas("c_photonE_wrongSign_matched");
    c_photonE_wrongSign_matched->Divide(1,3);
    c_photonE_wrongSign_matched->cd(1);
    h_photonE_wrongConver_Dsp_matched->Draw();
    c_photonE_wrongSign_matched->cd(2);
    h_photonE_wrongConver_Dsm_matched->Draw();
    c_photonE_wrongSign_matched->cd(3);
    h_photonE_wrongConver_matched->Draw();
    */
    // Save histograms in a ROOT file
    std::string filename_root=decay;
    filename_root+="_DsGamma_MBC.root";
    TFile *tfile=new TFile(filename_root.c_str(), "RECREATE");
    h_MBC_conver->Write();
    h_MBC_conver_matched->Write();
    h_MBC_conver_unmatched->Write();
    h_MBC_wrongConver->Write();
    h_MBC_wrongConver_matched->Write();
    h_MBC_wrongConver_unmatched->Write();
    h_MBC_wrongConver_strictUnmatched->Write();
    h_MBC_generic->Write();
    h_MBC_generic_veto->Write();
    h_MBC_generic_wrong->Write();
    h_MBC_continu->Write();
    h_MBC_physics->Write();
    h_photonE_wrongConver_matched->Write();
    h_photonE_wrongConver_unmatched->Write();
    h_angle_photon->Write();
    h_angle_wrongDs->Write();
    tfile->Close();
    
  }
  
  if (DeltaM_sideband)
  {
    /*
    h_DeltaM_generic->SetFillColor(kGreen); 
    h_DeltaM_continu->SetFillColor(kBlue);
    
    TCanvas *c_DeltaM_Stacked = new TCanvas("c_DeltaM_Stacked");
    THStack *s_DeltaM_sideband=new THStack("s_DeltaM_sideband", ("#deltam Sidebands in "+decay).c_str());
    s_DeltaM_sideband->Add(h_DeltaM_continu, "hist");
    s_DeltaM_sideband->Add(h_DeltaM_generic, "hist");
    double stack_ymax=h_DeltaM_generic->GetMaximum();
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
    std::string continu_string="Continuum MC: "; continu_string+=itoa(continuScale*h_DeltaM_continu->GetEntries()); continu_string+=" Events";
    std::string generic_string="Generic MC: "; generic_string+=itoa(genericScale*h_DeltaM_generic->GetEntries()); generic_string+=" Events";
    std::string physics_string="Data: "; physics_string+=itoa(h_DeltaM_physics->GetEntries()); physics_string+=" Events";
    legendDeltaM->AddEntry(h_DeltaM_continu, continu_string.c_str());
    legendDeltaM->AddEntry(h_DeltaM_generic, generic_string.c_str());
    legendDeltaM->AddEntry(h_DeltaM_physics, physics_string.c_str());
    legendDeltaM->SetFillColor(kWhite);
    legendDeltaM->Draw();
    
    //std::string filename_eps=decay+"_Sideband_DeltaM_005.eps";
    //c_DeltaM_Stacked->SaveAs(filename_eps.c_str());
    
    TCanvas *c_DeltaM_conver = new TCanvas("c_DeltaM_conver");
    c_DeltaM_conver->Divide(1,3);
    c_DeltaM_conver->cd(1);
    h_DeltaM_conver_Dsp->Draw();
    c_DeltaM_conver->cd(2);
    h_DeltaM_conver_Dsm->Draw();
    c_DeltaM_conver->cd(3);
    h_DeltaM_conver->Draw();
    ymax=h_DeltaM_conver->GetMaximum()*0.95;
    line=new TLine(xmin, ymax, xmin, 0); line->Draw();
    line=new TLine(xmax, ymax, xmax, 0); line->Draw();
    
    // Save histograms in a ROOT file
    std::string filename_root=decay;
    filename_root+="_DsGamma_DeltaM.root";
    TFile *tfile = new TFile(filename_root.c_str(), "RECREATE");
    h_DeltaM_conver_Dsp->Write();
    h_DeltaM_conver_Dsm->Write();
    h_DeltaM_conver->Write();
    h_DeltaM_wrongConver_Dsp->Write();
    h_DeltaM_wrongConver_Dsm->Write();
    h_DeltaM_wrongConver->Write();
    h_DeltaM_generic->Write();
    h_DeltaM_generic_veto->Write();
    h_DeltaM_generic_wrong->Write();
    h_DeltaM_continu->Write();
    h_DeltaM_physics->Write();
    tfile->Close();
    */
  }
  
  TCanvas *c_e9oe25Unf = new TCanvas("e9oe25Unf", "e9oe25Unf", 500, 1000);
  c_e9oe25Unf->Divide(1,5);
  c_e9oe25Unf->cd(1);
  h_e9oe25Unf_conver_Dsp->Draw();
  c_e9oe25Unf->cd(2);
  h_e9oe25Unf_conver_Dsm->Draw();
  c_e9oe25Unf->cd(3);
  h_e9oe25Unf_generic->Draw();
  c_e9oe25Unf->cd(4);
  h_e9oe25Unf_continu->Draw();
  c_e9oe25Unf->cd(5);
  h_e9oe25Unf_physics->Draw();
  
  return 0;
}
