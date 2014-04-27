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
#include "TPaveStats.h"
#include <set>
#include <map>

bool DeltaM_sideband=false;
bool mBC_sideband=false;

double pi=3.14159265358979;
double dsPlusMCut_center=1.96849, dsPlusMCut_range=0.01;
double deltaECut_center=0.012, deltaECut_range=0.02;
double mbcCut_center=2.112, mbcCut_range=0.008;
double deltaMCut_center=0.158, deltaMCut_range=0.01;
double diffD0Cut=-0.004;
double dPhiCutLess=0.07;

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
};

int DsTaggedAnalysis_GenericContinuumPhysics_pipi0eta()
{

  TFile *genericFile = new TFile("/nfs/cor/an2/souvik/Dataset47/DsTaggedProc_GenericMC_230474_232255.root");
  genericFile->cd("DsTaggedDecaysProc");
  TTree *genericTree = (TTree*)gDirectory->Get("nt1"); // pionFit
  
  TFile *continuFile = new TFile("/nfs/cor/an2/souvik/Dataset47/DsTaggedProc_ContinuumMC_230474_232255.root");
  continuFile->cd("DsTaggedDecaysProc");
  TTree *continuTree = (TTree*)gDirectory->Get("nt1"); // pionFit
  
  TFile *physicsFile = new TFile("/nfs/cor/an2/souvik/Dataset47/DsTaggedProc_Data_230474_232255.root");
  physicsFile->cd("DsTaggedDecaysProc");
  TTree *physicsTree = (TTree*)gDirectory->Get("nt1"); // pionFit
  
  NEvents genericEvents={0,0,0,0,0,0,0,0};
  NEvents continuEvents={0,0,0,0,0,0,0,0};
  NEvents physicsEvents={0,0,0,0,0,0,0,0};
  
  EventNumber genericNumber={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  EventNumber continuNumber={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  EventNumber physicsNumber={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  float runNumber_generic, eventNumber_generic;
  float runNumber_continu, eventNumber_continu;
  float runNumber_physics, eventNumber_physics;
  float dsPlusM_generic, dsPlusCharge_generic, DeltaE_generic, MBC_generic, DeltaM_generic, eeMass_generic, decayMode_generic;
  float dsPlusM_continu, dsPlusCharge_continu, DeltaE_continu, MBC_continu, DeltaM_continu, eeMass_continu, decayMode_continu;
  float dsPlusM_physics, dsPlusCharge_physics, DeltaE_physics, MBC_physics, DeltaM_physics, eeMass_physics, decayMode_physics;
  float d0_e_generic, d0_p_generic, z0_e_generic, z0_p_generic, px_e_generic, py_e_generic, pz_e_generic, px_p_generic, py_p_generic, pz_p_generic, E_e_generic, E_p_generic, curv_e_generic, curv_p_generic;
  float d0_e_continu, d0_p_continu, z0_e_continu, z0_p_continu, px_e_continu, py_e_continu, pz_e_continu, px_p_continu, py_p_continu, pz_p_continu, E_e_continu, E_p_continu, curv_e_continu, curv_p_continu;
  float d0_e_physics, d0_p_physics, z0_e_physics, z0_p_physics, px_e_physics, py_e_physics, pz_e_physics, px_p_physics, py_p_physics, pz_p_physics, E_e_physics, E_p_physics, curv_e_physics, curv_p_physics;
  float px_e_generic_MC, py_e_generic_MC, pz_e_generic_MC, E_e_generic_MC, px_p_generic_MC, py_p_generic_MC, pz_p_generic_MC, E_p_generic_MC;
  float px_e_continu_MC, py_e_continu_MC, pz_e_continu_MC, E_e_continu_MC, px_p_continu_MC, py_p_continu_MC, pz_p_continu_MC, E_p_continu_MC;
  
  typedef std::map<int, float> DecayMap;
  DecayMap decayFrequency_generic, decayFrequency_continu, decayFrequency_physics;
  
  TH1D *h_dsPlusM_generic = new TH1D("h_dsPlusM_generic", "m_{D_{S}^{+}} generic Sample; GeV", 100, 1.9, 2.1); h_dsPlusM_generic->SetLineColor(kRed);
  TH1D *h_dsPlusM_continu = new TH1D("h_dsPlusM_continu", "m_{D_{S}^{+}} continuum Background Sample; GeV", 100, 1.9, 2.1); h_dsPlusM_continu->SetLineColor(kBlue);
  TH1D *h_dsPlusM_physics = new TH1D("h_dsPlusM_physics", "m_{D_{S}^{+}} Data; GeV", 100, 1.9, 2.1); h_dsPlusM_physics->SetLineColor(kGreen);
  TH1D *h_DeltaE_generic = new TH1D("h_DeltaE_generic", "#DeltaE generic Sample; GeV", 100, -0.1, 0.2); h_DeltaE_generic->SetLineColor(kRed);
  TH1D *h_DeltaE_continu = new TH1D("h_DeltaE_continu", "#DeltaE continuum Background Sample; GeV", 100, -0.1, 0.2); h_DeltaE_continu->SetLineColor(kBlue);
  TH1D *h_DeltaE_physics = new TH1D("h_DeltaE_physics", "#DeltaE Data; GeV", 100, -0.1, 0.2); h_DeltaE_physics->SetLineColor(kGreen);
  TH1D *h_MBC_generic = new TH1D("h_MBC_generic", "m_{BC} generic Sample; GeV", 100, 2., 2.2); h_MBC_generic->SetLineColor(kRed);
  TH1D *h_MBC_continu = new TH1D("h_MBC_continu", "m_{BC} continuum Background; GeV", 100, 2., 2.2); h_MBC_continu->SetLineColor(kBlue);
  TH1D *h_MBC_physics = new TH1D("h_MBC_physics", "m_{BC} Data; GeV", 100, 2., 2.2); h_MBC_physics->SetLineColor(kGreen);
  TH1D *h_DeltaM_generic = new TH1D("h_DeltaM_generic", "#deltaM generic Sample", 100, 0.0, 0.5); h_DeltaM_generic->SetLineColor(kRed);
  TH1D *h_DeltaM_continu = new TH1D("h_DeltaM_continu", "#deltaM continuum Background Sample; GeV", 100, 0.0, 0.5); h_DeltaM_continu->SetLineColor(kBlue);
  TH1D *h_DeltaM_physics = new TH1D("h_DeltaM_physics", "#deltaM Data; GeV", 100, 0.0, 0.5); h_DeltaM_physics->SetLineColor(kGreen);
  
  TH2D *h_DeltaE_MBC_generic = new TH2D("h_DeltaE_MBC_generic", "h_DeltaE_MBC_generic", 100, -0.1, 0.2, 100, 2., 2.2);
  TH2D *h_DeltaE_MBC_continu = new TH2D("h_DeltaE_MBC_continu", "h_DeltaE_MBC_continu", 100, -0.1, 0.2, 100, 2., 2.2);
  TH2D *h_DeltaE_MBC_physics = new TH2D("h_DeltaE_MBC_physics", "h_DeltaE_MBC_physics", 100, -0.1, 0.2, 100, 2., 2.2);
  TH2D *h_DeltaE_DeltaM_generic = new TH2D("h_DeltaE_DeltaM_generic", "h_DeltaE_DeltaM_generic", 100, -0.1, 0.2, 100, 0.0, 0.2);
  TH2D *h_DeltaE_DeltaM_continu = new TH2D("h_DeltaE_DeltaM_continu", "h_DeltaE_DeltaM_continu", 100, -0.1, 0.2, 100, 0.0, 0.2);
  TH2D *h_DeltaE_DeltaM_physics = new TH2D("h_DeltaE_DeltaM_physics", "h_DeltaE_DeltaM_physics", 100, -0.1, 0.2, 100, 0.0, 0.2);
  TH2D *h_MBC_DeltaM_generic = new TH2D("h_MBC_DeltaM_generic", "h_MBC_DeltaM_generic", 100, 2., 2.2, 100, 0.0, 0.2);
  TH2D *h_MBC_DeltaM_continu = new TH2D("h_MBC_DeltaM_continu", "h_MBC_DeltaM_continu", 100, 2., 2.2, 100, 0.0, 0.2);
  TH2D *h_MBC_DeltaM_physics = new TH2D("h_MBC_DeltaM_physics", "h_MBC_DeltaM_physics", 100, 2., 2.2, 100, 0.0, 0.2);
  
  TH1D *h_d0_e_generic = new TH1D("h_d0_e_generic", "d0_e", 50, -0.01, 0.01); h_d0_e_generic->SetLineColor(kRed);
  TH1D *h_d0_e_continu = new TH1D("h_d0_e_continu", "d0_e", 50, -0.01, 0.01); h_d0_e_continu->SetLineColor(kBlue);
  TH1D *h_d0_e_physics = new TH1D("h_d0_e_physics", "d0_e", 50, -0.01, 0.01); h_d0_e_physics->SetLineColor(kGreen);
  TH1D *h_d0_p_generic = new TH1D("h_d0_p_generic", "d0_p", 50, -0.01, 0.01); h_d0_p_generic->SetLineColor(kRed);
  TH1D *h_d0_p_continu = new TH1D("h_d0_p_continu", "d0_p", 50, -0.01, 0.01); h_d0_p_continu->SetLineColor(kBlue);
  TH1D *h_d0_p_physics = new TH1D("h_d0_p_physics", "d0_p", 50, -0.01, 0.01); h_d0_p_physics->SetLineColor(kGreen);
  TH1D *h_diffD0_generic = new TH1D("h_diffD0_generic", "#Deltad_{0} generic Sample; m", 50, -0.01, 0.01); h_diffD0_generic->SetLineColor(kRed);
  TH1D *h_diffD0_continu = new TH1D("h_diffD0_continu", "#Deltad_{0} continuum Background Sample; m", 50, -0.01, 0.01); h_diffD0_continu->SetLineColor(kBlue);
  TH1D *h_diffD0_physics = new TH1D("h_diffD0_physics", "#Deltad_{0} Data; m", 50, -0.01, 0.01); h_diffD0_physics->SetLineColor(kGreen);
  TH1D *h_dPhi_generic = new TH1D("h_dPhi_generic", "#Delta#Phi generic Sample", 50, -2., 2.); h_dPhi_generic->SetLineColor(kRed);
  TH1D *h_dPhi_continu = new TH1D("h_dPhi_continu", "#Delta#Phi continuum Background Sample", 50, -2., 2.); h_dPhi_continu->SetLineColor(kBlue);
  TH1D *h_dPhi_physics = new TH1D("h_dPhi_physics", "#Delta#Phi Data", 50, -2., 2.); h_dPhi_physics->SetLineColor(kGreen);
  TH2D *h_dPhi_ee_generic = new TH2D("h_dPhi_ee_generic", "dPhi_ee_generic", 50, 0., 0.2, 50, -2., 2.);
  TH2D *h_dPhi_ee_continu = new TH2D("h_dPhi_ee_continu", "dPhi_ee_continu", 50, 0., 0.2, 50, -2., 2.);
  TH2D *h_dPhi_ee_physics = new TH2D("h_dPhi_ee_physics", "dPhi_ee_physics", 50, 0., 0.2, 50, -2., 2.);
  TH2D *h_d0_phi_e_generic = new TH2D("h_d0_phi_e_generic", "h_d0_phi_e_generic", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_p_generic = new TH2D("h_d0_phi_p_generic", "h_d0_phi_p_generic", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_e_continu = new TH2D("h_d0_phi_e_continu", "h_d0_phi_e_continu", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_p_continu = new TH2D("h_d0_phi_p_continu", "h_d0_phi_p_continu", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_e_physics = new TH2D("h_d0_phi_e_physics", "h_d0_phi_e_physics", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_p_physics = new TH2D("h_d0_phi_p_physics", "h_d0_phi_p_physics", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_dPhi_diffD0_generic = new TH2D("h_dPhi_diffD0_generic", "#Delta#Phi vs #Deltad_{0} generic Sample; ; m", 50, -2., 2., 50, -0.01, 0.01); h_dPhi_diffD0_generic->SetLineColor(kRed);
  TH2D *h_dPhi_diffD0_continu = new TH2D("h_dPhi_diffD0_continu", "#Delta#Phi vs #Deltad_{0} continuum Background Sample; ; m", 50, -2., 2., 50, -0.01, 0.01); h_dPhi_diffD0_continu->SetLineColor(kBlue);
  TH2D *h_dPhi_diffD0_physics = new TH2D("h_dPhi_diffD0_physics", "#Delta#Phi vs #Deltad_{0} Data; ; m", 50, -2., 2., 50, -0.01, 0.01); h_dPhi_diffD0_physics->SetLineColor(kGreen);
  TH1D *h_R_generic = new TH1D("h_R_generic", "h_R_generic", 50, -.4, .4); h_R_generic->SetLineColor(kRed);
  TH1D *h_R_continu = new TH1D("h_R_continu", "h_R_continu", 50, -.4, .4); h_R_continu->SetLineColor(kBlue);
  TH1D *h_R_physics = new TH1D("h_R_physics", "h_R_physics", 50, -.4, .4); h_R_physics->SetLineColor(kGreen);
  TH1D *h_ee_generic = new TH1D("h_ee_generic", "h_ee_generic", 20, 0.0, 0.2); h_ee_generic->SetLineColor(kRed);
  TH1D *h_ee_continu = new TH1D("h_ee_continu", "h_ee_continu", 20, 0.0, 0.2); h_ee_continu->SetLineColor(kBlue);
  TH1D *h_ee_physics = new TH1D("h_ee_physics", "h_ee_physics", 20, 0.0, 0.2); h_ee_physics->SetLineColor(kGreen);
  TH2D *h_electronE = new TH2D("h_electronE", "Electron Energy Resolution; Reco Electron Energy (GeV)", 100, 0.0, 0.16, 100, -0.03, 0.03);
  TH2D *h_electronPx = new TH2D("h_electronPx", "Electron Px Resolution; Reco Electron Px (GeV)", 100, 0.0, 0.14, 100, -0.03, 0.03); 
  TH2D *h_electronPy = new TH2D("h_electronPy", "Electron Py Resolution; Reco Electron Py (GeV)", 100, 0.0, 0.14, 100, -0.03, 0.03); 
  TH2D *h_electronPz = new TH2D("h_electronPz", "Electron Pz Resolution; Reco Electron Pz (GeV)", 100, 0.0, 0.14, 100, -0.03, 0.03);
  
  genericTree->SetBranchAddress("Run", &(runNumber_generic));
  genericTree->SetBranchAddress("Event", &(eventNumber_generic));
  genericTree->SetBranchAddress("dsPlusM", &(dsPlusM_generic));
  genericTree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_generic));
  genericTree->SetBranchAddress("DecayMode", &(decayMode_generic));
  genericTree->SetBranchAddress("DeltaE", &(DeltaE_generic));
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
  genericTree->SetBranchAddress("keeMass_reco", &(eeMass_generic));
  
  continuTree->SetBranchAddress("Run", &(runNumber_continu));
  continuTree->SetBranchAddress("Event", &(eventNumber_continu));
  continuTree->SetBranchAddress("dsPlusM", &(dsPlusM_continu));
  continuTree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_continu));
  continuTree->SetBranchAddress("DecayMode", &(decayMode_continu));
  continuTree->SetBranchAddress("DeltaE", &(DeltaE_continu));
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
  continuTree->SetBranchAddress("keeMass_reco", &(eeMass_continu));
  
  physicsTree->SetBranchAddress("Run", &(runNumber_physics));
  physicsTree->SetBranchAddress("Event", &(eventNumber_physics));
  physicsTree->SetBranchAddress("dsPlusM", &(dsPlusM_physics));
  physicsTree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_physics));
  physicsTree->SetBranchAddress("DecayMode", &(decayMode_physics));
  physicsTree->SetBranchAddress("DeltaE", &(DeltaE_physics));
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
  physicsTree->SetBranchAddress("keeMass_reco", &(eeMass_physics));
  
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
    
    if (decayMode_generic==441)
    {
      if (runNumber_generic!=genericNumber.tagCut_run || eventNumber_generic!=genericNumber.tagCut_event)
      {      
        ++genericEvents.tagCut;
        genericNumber.tagCut_run=runNumber_generic;
        genericNumber.tagCut_event=eventNumber_generic;
      }
      h_dsPlusM_generic->Fill(dsPlusM_generic);
      if (dsPlusMCut(dsPlusM_generic))
      {
        if (runNumber_generic!=genericNumber.dsPlusMCut_run || eventNumber_generic!=genericNumber.dsPlusMCut_event)
        {
          ++genericEvents.dsPlusMCut;
          genericNumber.dsPlusMCut_run=runNumber_generic;
          genericNumber.dsPlusMCut_event=eventNumber_generic;
        }
        h_DeltaE_generic->Fill(DeltaE_generic);
 
        // Optimization plots
        h_DeltaE_MBC_generic->Fill(DeltaE_generic, MBC_generic);
        h_DeltaE_DeltaM_generic->Fill(DeltaE_generic, DeltaM_generic);
        h_MBC_DeltaM_generic->Fill(MBC_generic, DeltaM_generic);
        //if (DeltaECut(DeltaE_generic))
        {
          if (runNumber_generic!=genericNumber.deltaECut_run || eventNumber_generic!=genericNumber.deltaECut_event)
          {
            ++genericEvents.deltaECut;
            genericNumber.deltaECut_run=runNumber_generic;
            genericNumber.deltaECut_event=eventNumber_generic;
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
              h_d0_e_generic->Fill(d0_e_generic);
              h_d0_p_generic->Fill(d0_p_generic);
              h_diffD0_generic->Fill(d0_e_generic-d0_p_generic);
              h_dPhi_diffD0_generic->Fill(dPhi, d0_e_generic-d0_p_generic);
              h_d0_phi_e_generic->Fill(phi_e, d0_e_generic);
              h_d0_phi_p_generic->Fill(phi_p, d0_p_generic);
              if (dD0(d0_e_generic-d0_p_generic))
              {
                if (runNumber_generic!=genericNumber.diffD0Cut_run || eventNumber_generic!=genericNumber.diffD0Cut_event) 
                {
                  ++genericEvents.diffD0Cut;
                  genericNumber.diffD0Cut_run=runNumber_generic;
                  genericNumber.diffD0Cut_event=eventNumber_generic;
                }
                h_dPhi_generic->Fill(dPhi);
                h_dPhi_ee_generic->Fill(eeMass_generic, dPhi);
                if (dPhiCut(dPhi))
                {
                  if (runNumber_generic!=genericNumber.dPhiCut_run || eventNumber_generic!=genericNumber.dPhiCut_event) 
                  {
                    ++genericEvents.dPhiCut;
                    genericNumber.dPhiCut_run=runNumber_generic;
                    genericNumber.dPhiCut_event=eventNumber_generic;                   
                    h_electronE->Fill(E_e_generic, E_e_generic-E_e_generic_MC);
                    h_electronE->Fill(E_p_generic, E_p_generic-E_p_generic_MC);
                    h_electronPx->Fill(px_e_generic, px_e_generic-px_e_generic_MC);
                    h_electronPx->Fill(px_p_generic, px_p_generic-px_p_generic_MC);
                    h_electronPy->Fill(py_e_generic, py_e_generic-py_e_generic_MC);
                    h_electronPy->Fill(py_p_generic, py_p_generic-py_p_generic_MC);
                    h_electronPz->Fill(pz_e_generic, pz_e_generic-pz_e_generic_MC);
                    h_electronPz->Fill(pz_p_generic, pz_p_generic-pz_p_generic_MC);
                  }
                  h_ee_generic->Fill(eeMass_generic);
                  double r1=-1/(2*curv_e_generic);
                  double r2=1/(2*curv_p_generic);
                  double R=convert(r1, r2, d0_e_generic, d0_p_generic, phi_e, phi_p);
                  h_R_generic->Fill(R);
                  if (dsPlusCharge_generic<0) oppRecon_generic+=1;
                  if (DeltaM_sideband) h_DeltaM_generic->Fill(DeltaM_generic);
                  if (mBC_sideband) h_MBC_generic->Fill(MBC_generic);
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
    
    if (decayMode_continu==441)
    {
      if (runNumber_continu!=continuNumber.tagCut_run || eventNumber_continu!=continuNumber.tagCut_event)
      {      
        ++continuEvents.tagCut;
        continuNumber.tagCut_run=runNumber_continu;
        continuNumber.tagCut_event=eventNumber_continu;
      }
      h_dsPlusM_continu->Fill(dsPlusM_continu);
      if (dsPlusMCut(dsPlusM_continu))
      {
        if (runNumber_continu!=continuNumber.dsPlusMCut_run || eventNumber_continu!=continuNumber.dsPlusMCut_event)
        {
          ++continuEvents.dsPlusMCut;
          continuNumber.dsPlusMCut_run=runNumber_continu;
          continuNumber.dsPlusMCut_event=eventNumber_continu;
        }
        h_DeltaE_continu->Fill(DeltaE_continu);
 
        // Optimization plots
        h_DeltaE_MBC_continu->Fill(DeltaE_continu, MBC_continu);
        h_DeltaE_DeltaM_continu->Fill(DeltaE_continu, DeltaM_continu);
        h_MBC_DeltaM_continu->Fill(MBC_continu, DeltaM_continu);
        //if (DeltaECut(DeltaE_continu))
        {
          if (runNumber_continu!=continuNumber.deltaECut_run || eventNumber_continu!=continuNumber.deltaECut_event)
          {
            ++continuEvents.deltaECut;
            continuNumber.deltaECut_run=runNumber_continu;
            continuNumber.deltaECut_event=eventNumber_continu;
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
              h_d0_e_continu->Fill(d0_e_continu);
              h_d0_p_continu->Fill(d0_p_continu);
              h_diffD0_continu->Fill(d0_e_continu-d0_p_continu);
              h_dPhi_diffD0_continu->Fill(dPhi, d0_e_continu-d0_p_continu);
              h_d0_phi_e_continu->Fill(phi_e, d0_e_continu);
              h_d0_phi_p_continu->Fill(phi_p, d0_p_continu);
              if (dD0(d0_e_continu-d0_p_continu))
              {
                if (runNumber_continu!=continuNumber.diffD0Cut_run || eventNumber_continu!=continuNumber.diffD0Cut_event) 
                {
                  ++continuEvents.diffD0Cut;
                  continuNumber.diffD0Cut_run=runNumber_continu;
                  continuNumber.diffD0Cut_event=eventNumber_continu;
                }
                h_dPhi_continu->Fill(dPhi);
                h_dPhi_ee_continu->Fill(eeMass_continu, dPhi);
                if (dPhiCut(dPhi))
                {
                  if (runNumber_continu!=continuNumber.dPhiCut_run || eventNumber_continu!=continuNumber.dPhiCut_event) 
                  {
                    ++continuEvents.dPhiCut;
                    continuNumber.dPhiCut_run=runNumber_continu;
                    continuNumber.dPhiCut_event=eventNumber_continu;
                  }
                  h_ee_continu->Fill(eeMass_continu);
                  double r1=-1/(2*curv_e_continu);
                  double r2=1/(2*curv_p_continu);
                  double R=convert(r1, r2, d0_e_continu, d0_p_continu, phi_e, phi_p);
                  h_R_continu->Fill(R);
                  if (dsPlusCharge_continu<0) oppRecon_continu+=1;
                  if (DeltaM_sideband) h_DeltaM_continu->Fill(DeltaM_continu);
                  if (mBC_sideband) h_MBC_continu->Fill(MBC_continu);
                }
              }
            }
          }
        }
      }
    }
    }
  }
  
  
  int nphysicsEvents=physicsTree->GetEntries();
  int oppRecon_physics=0;
  for (int i=0; i<nphysicsEvents; ++i)
  {
    physicsTree->GetEvent(i);
    
    double phi_e=atan2(py_e_physics, px_e_physics);
    double phi_p=atan2(py_p_physics, px_p_physics);
    double dPhi=deltaPhi(phi_e, phi_p);
    
    if (runNumber_physics!=physicsNumber.noCut_run || eventNumber_physics!=physicsNumber.noCut_event)
    {      
      if (decayMode_physics>-1) decayFrequency_physics[int(decayMode_physics)]+=1;
      ++physicsEvents.noCut;
      physicsNumber.noCut_run=runNumber_physics;
      physicsNumber.noCut_event=eventNumber_physics;
    }
    
    if (fabs(d0_e_physics)<0.005 && fabs(d0_p_physics)<0.005 && 
        fabs(z0_e_physics)<0.05 && fabs(z0_p_physics)<0.05)
    {
    
    if (decayMode_physics==441)
    {
      if (runNumber_physics!=physicsNumber.tagCut_run || eventNumber_physics!=physicsNumber.tagCut_event)
      {      
        ++physicsEvents.tagCut;
        physicsNumber.tagCut_run=runNumber_physics;
        physicsNumber.tagCut_event=eventNumber_physics;
      }
      h_dsPlusM_physics->Fill(dsPlusM_physics);
      if (dsPlusMCut(dsPlusM_physics))
      {
        if (runNumber_physics!=physicsNumber.dsPlusMCut_run || eventNumber_physics!=physicsNumber.dsPlusMCut_event)
        {
          ++physicsEvents.dsPlusMCut;
          physicsNumber.dsPlusMCut_run=runNumber_physics;
          physicsNumber.dsPlusMCut_event=eventNumber_physics;
        }
        h_DeltaE_physics->Fill(DeltaE_physics); 
        // Optimization plots
        h_DeltaE_MBC_physics->Fill(DeltaE_physics, MBC_physics);
        h_DeltaE_DeltaM_physics->Fill(DeltaE_physics, DeltaM_physics);
        h_MBC_DeltaM_physics->Fill(MBC_physics, DeltaM_physics);
        //if (DeltaECut(DeltaE_physics))
        {
          if (runNumber_physics!=physicsNumber.deltaECut_run || eventNumber_physics!=physicsNumber.deltaECut_event)
          {
            ++physicsEvents.deltaECut;
            physicsNumber.deltaECut_run=runNumber_physics;
            physicsNumber.deltaECut_event=eventNumber_physics;
          }
          if (!mBC_sideband) h_MBC_physics->Fill(MBC_physics);
          if (mBC_sideband || MBCCut(MBC_physics))
          {
            if (runNumber_physics!=physicsNumber.mbcCut_run || eventNumber_physics!=physicsNumber.mbcCut_event)
            {
              ++physicsEvents.mbcCut;
              physicsNumber.mbcCut_run=runNumber_physics;
              physicsNumber.mbcCut_event=eventNumber_physics;
            }
            if (!DeltaM_sideband) h_DeltaM_physics->Fill(DeltaM_physics);
            if (DeltaM_sideband || DeltaMCut(DeltaM_physics))
            {
              if (runNumber_physics!=physicsNumber.deltaMCut_run || eventNumber_physics!=physicsNumber.deltaMCut_event) 
              {
                ++physicsEvents.deltaMCut;
                physicsNumber.deltaMCut_run=runNumber_physics;
                physicsNumber.deltaMCut_event=eventNumber_physics;
              }
              h_d0_e_physics->Fill(d0_e_physics);
              h_d0_p_physics->Fill(d0_p_physics);
              h_diffD0_physics->Fill(d0_e_physics-d0_p_physics);
              h_dPhi_diffD0_physics->Fill(dPhi, d0_e_physics-d0_p_physics);
              h_d0_phi_e_physics->Fill(phi_e, d0_e_physics);
              h_d0_phi_p_physics->Fill(phi_p, d0_p_physics);
              if (dD0(d0_e_physics-d0_p_physics))
              {
                if (runNumber_physics!=physicsNumber.diffD0Cut_run || eventNumber_physics!=physicsNumber.diffD0Cut_event) 
                {
                  ++physicsEvents.diffD0Cut;
                  physicsNumber.diffD0Cut_run=runNumber_physics;
                  physicsNumber.diffD0Cut_event=eventNumber_physics;
                }
                h_dPhi_physics->Fill(dPhi);
                h_dPhi_ee_physics->Fill(eeMass_physics, dPhi);
                if (dPhiCut(dPhi))
                {
                  if (runNumber_physics!=physicsNumber.dPhiCut_run || eventNumber_physics!=physicsNumber.dPhiCut_event) 
                  {
                    ++physicsEvents.dPhiCut;
                    physicsNumber.dPhiCut_run=runNumber_physics;
                    physicsNumber.dPhiCut_event=eventNumber_physics;
                    std::cout<<runNumber_physics<<" "<<eventNumber_physics<<std::endl;
                  }
                  h_ee_physics->Fill(eeMass_physics);
                  double r1=-1/(2*curv_e_physics);
                  double r2=1/(2*curv_p_physics);
                  double R=convert(r1, r2, d0_e_physics, d0_p_physics, phi_e, phi_p);
                  h_R_physics->Fill(R);
                  if (dsPlusCharge_physics<0) oppRecon_physics+=1;
                  if (DeltaM_sideband) h_DeltaM_physics->Fill(DeltaM_physics);
                  if (mBC_sideband) h_MBC_physics->Fill(MBC_physics);
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
  std::cout<<"Number of generic events after pipi0eta tag = "<<genericEvents.tagCut<<std::endl;
  std::cout<<"Number of generic events after dsPlusMCut = "<<genericEvents.dsPlusMCut<<std::endl;
  std::cout<<"Number of generic events after deltaECut = "<<genericEvents.deltaECut<<std::endl;
  std::cout<<"Number of generic events after mbcCut = "<<genericEvents.mbcCut<<std::endl;
  std::cout<<"Number of generic events after deltaMCut = "<<genericEvents.deltaMCut<<std::endl;
  std::cout<<"Number of generic events after diffD0Cut = "<<genericEvents.diffD0Cut<<std::endl;
  std::cout<<"Number of generic events after dPhiCut = "<<genericEvents.dPhiCut<<std::endl;
  std::cout<<"Number of D_s- reconstructions remaining in generic = "<<oppRecon_generic<<std::endl;
  
  std::cout<<"Number of continuum candidates = "<<ncontinuEvents<<std::endl;
  for (DecayMap::iterator i_decay=decayFrequency_continu.begin(); i_decay!=decayFrequency_continu.end(); ++i_decay)
  {
    std::cout<<"continuum decay mode "<<i_decay->first<<" had "<<i_decay->second<<" decays"<<std::endl;
  }
  std::cout<<"Number of continu events after ntuplizer = "<<continuEvents.noCut<<std::endl;
  std::cout<<"Number of continu events after pipi0eta tag = "<<continuEvents.tagCut<<std::endl;
  std::cout<<"Number of continu events after dsPlusMCut = "<<continuEvents.dsPlusMCut<<std::endl;
  std::cout<<"Number of continu events after deltaECut = "<<continuEvents.deltaECut<<std::endl;
  std::cout<<"Number of continu events after mbcCut = "<<continuEvents.mbcCut<<std::endl;
  std::cout<<"Number of continu events after deltaMCut = "<<continuEvents.deltaMCut<<std::endl;
  std::cout<<"Number of continu events after diffD0Cut = "<<continuEvents.diffD0Cut<<std::endl;
  std::cout<<"Number of continu events after dPhiCut = "<<continuEvents.dPhiCut<<std::endl;
  std::cout<<"Number of D_s- reconstructions remaining in continu = "<<oppRecon_continu<<std::endl;
  
  std::cout<<"Number of data candidates = "<<nphysicsEvents<<std::endl;
  for (DecayMap::iterator i_decay=decayFrequency_physics.begin(); i_decay!=decayFrequency_physics.end(); ++i_decay)
  {
    std::cout<<"physicssion decay mode "<<i_decay->first<<" had "<<i_decay->second<<" decays"<<std::endl;
  }
  std::cout<<"Number of physics events after ntuplizer = "<<physicsEvents.noCut<<std::endl;
  std::cout<<"Number of physics events after pipi0eta tag = "<<physicsEvents.tagCut<<std::endl;
  std::cout<<"Number of physics events after dsPlusMCut = "<<physicsEvents.dsPlusMCut<<std::endl;
  std::cout<<"Number of physics events after deltaECut = "<<physicsEvents.deltaECut<<std::endl;
  std::cout<<"Number of physics events after mbcCut = "<<physicsEvents.mbcCut<<std::endl;
  std::cout<<"Number of physics events after deltaMCut = "<<physicsEvents.deltaMCut<<std::endl;
  std::cout<<"Number of physics events after diffD0Cut = "<<physicsEvents.diffD0Cut<<std::endl;
  std::cout<<"Number of physics events after dPhiCut = "<<physicsEvents.dPhiCut<<std::endl;
  std::cout<<"Number of D_s- reconstructions remaining in physics = "<<oppRecon_physics<<std::endl;
  
  TLine *line;
  float xmin, xmax, ymin, ymax;
  
  gROOT->SetStyle("Plain");
  
  TCanvas *dsPlusM = new TCanvas("dsPlusM", "", 500, 1500);  
  xmin=dsPlusMCut_center-dsPlusMCut_range;
  xmax=dsPlusMCut_center+dsPlusMCut_range;
  dsPlusM->Divide(1,3);  
  dsPlusM->cd(1);
  h_dsPlusM_generic->Draw();
  ymax=(h_dsPlusM_generic->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("dsPlusM_generic_pipi0eta.png"); 
  dsPlusM->cd(2);
  h_dsPlusM_continu->Draw();
  ymax=(h_dsPlusM_continu->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("dsPlusM_continu_pipi0eta.png");
  dsPlusM->cd(3);
  h_dsPlusM_physics->Draw();
  ymax=(h_dsPlusM_physics->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("dsPlusM_physics_pipi0eta.png");
  
  /*
  
  TCanvas *DeltaE_MBC = new TCanvas("DeltaE_MBC", "", 500, 1500);
  xmin=-deltaECut_range; xmax=deltaECut_range;
  ymin=mbcCut_center-mbcCut_range; ymax=mbcCut_center+mbcCut_range;
  DeltaE_MBC->Divide(1,3);
  DeltaE_MBC->cd(1);
  h_DeltaE_MBC_generic->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  DeltaE_MBC->cd(2);
  h_DeltaE_MBC_continu->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  DeltaE_MBC->cd(3);
  h_DeltaE_MBC_physics->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  
  TCanvas *DeltaE_DeltaM = new TCanvas("DeltaE_DeltaM", "", 500, 1500);
  xmin=-deltaECut_range; xmax=deltaECut_range;
  ymin=deltaMCut_center-deltaMCut_range; ymax=deltaMCut_center+deltaMCut_range;
  DeltaE_DeltaM->Divide(1,3);
  DeltaE_DeltaM->cd(1);
  h_DeltaE_DeltaM_generic->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  DeltaE_DeltaM->cd(2);
  h_DeltaE_DeltaM_continu->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  DeltaE_DeltaM->cd(3);
  h_DeltaE_DeltaM_physics->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  
  TCanvas *MBC_DeltaM = new TCanvas("MBC_DeltaM", "", 500, 1500);
  xmin=mbcCut_center-mbcCut_range; xmax=mbcCut_center+mbcCut_range;
  ymin=deltaMCut_center-deltaMCut_range; ymax=deltaMCut_center+deltaMCut_range;
  MBC_DeltaM->Divide(1,3);
  MBC_DeltaM->cd(1);
  h_MBC_DeltaM_generic->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  gPad->Print("mBC_DeltaM_generic_pipi0eta.png");
  MBC_DeltaM->cd(2);
  h_MBC_DeltaM_continu->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  gPad->Print("mBC_DeltaM_continu_pipi0eta.png");
  MBC_DeltaM->cd(3);
  h_MBC_DeltaM_physics->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  gPad->Print("mBC_DeltaM_physics_pipi0eta.png");
  
  */
  
  TCanvas *DeltaE = new TCanvas("DeltaE", "", 500, 1500);
  xmin=deltaECut_center-deltaECut_range; xmax=deltaECut_center+deltaECut_range;
  DeltaE->Divide(1,3);
  DeltaE->cd(1);
  h_DeltaE_generic->Draw();
  ymax=(h_DeltaE_generic->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("DeltaE_generic_pipi0eta.png");
  DeltaE->cd(2);
  h_DeltaE_continu->Draw();
  ymax=(h_DeltaE_continu->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("DeltaE_continu_pipi0eta.png");
  DeltaE->cd(3);
  h_DeltaE_physics->Draw();
  ymax=(h_DeltaE_physics->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("DeltaE_physics_pipi0eta.png");
  
  TCanvas *MBC = new TCanvas("MBC", "", 500, 1500);
  xmin=mbcCut_center-mbcCut_range; xmax=mbcCut_center+mbcCut_range;
  MBC->Divide(1,3);
  MBC->cd(1);
  h_MBC_generic->Draw();
  ymax=(h_MBC_generic->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("mBC_generic_pipi0eta.png");
  MBC->cd(2);
  h_MBC_continu->Draw();
  ymax=(h_MBC_continu->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("mBC_continu_pipi0eta.png");
  MBC->cd(3);
  h_MBC_physics->Draw();
  ymax=(h_MBC_physics->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("mBC_physics_pipi0eta.png");
  
  TCanvas *DeltaM = new TCanvas("DeltaM", "", 500, 1500);
  xmin=deltaMCut_center-deltaMCut_range; xmax=deltaMCut_center+deltaMCut_range;
  DeltaM->Divide(1,3);
  DeltaM->cd(1);
  h_DeltaM_generic->Draw();
  ymax=(h_DeltaM_generic->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("DeltaM_generic_pipi0eta.png");
  DeltaM->cd(2);
  h_DeltaM_continu->Draw();
  ymax=(h_DeltaM_continu->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("DeltaM_continu_pipi0eta.png");
  DeltaM->cd(3);
  h_DeltaM_physics->Draw();
  ymax=(h_DeltaM_physics->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("DeltaM_physics_pipi0eta.png");
  
  TCanvas *d0_phi = new TCanvas("d0_phi");
  d0_phi->Divide(1,4);
  d0_phi->cd(1);
  h_d0_phi_e_generic->Draw();
  d0_phi->cd(2);
  h_d0_phi_p_generic->Draw();
  d0_phi->cd(3);
  h_d0_phi_e_continu->Draw();
  d0_phi->cd(4);
  h_d0_phi_p_continu->Draw();
    
  TCanvas *diffD0 = new TCanvas("diffD0");
  diffD0->Divide(1,3);
  diffD0->cd(1);
  h_diffD0_generic->Draw();
  ymax=(h_diffD0_generic->GetMaximum())*0.75;
  line = new TLine(diffD0Cut,ymax,diffD0Cut,0); line->Draw();
  diffD0->cd(2);
  h_diffD0_continu->Draw();
  ymax=(h_diffD0_continu->GetMaximum())*0.75;
  line = new TLine(diffD0Cut,ymax,diffD0Cut,0); line->Draw();
  diffD0->cd(3);
  h_diffD0_physics->Draw();
  ymax=(h_diffD0_physics->GetMaximum())*0.75;
  line = new TLine(diffD0Cut,ymax,diffD0Cut,0); line->Draw();
  
  TCanvas *dPhi_diffD0 = new TCanvas ("dPhi_diffD0");
  xmin=-2; xmax=dPhiCutLess;
  ymin=diffD0Cut; ymax=0.01;
  dPhi_diffD0->Divide(1,3);
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
  h_dPhi_diffD0_physics->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  
  TCanvas *dPhi = new TCanvas("dPhi", "", 500, 1500);
  dPhi->Divide(1,3);
  dPhi->cd(1);
  h_dPhi_generic->Draw();
  line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
  gPad->Print("dPhi_generic_pipi0eta.png");
  dPhi->cd(2);  
  h_dPhi_continu->Draw();
  ymax=(h_dPhi_continu->GetMaximum())*0.75;
  line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
  gPad->Print("dPhi_continu_pipi0eta.png");
  dPhi->cd(3);
  h_dPhi_physics->Draw();
  ymax=(h_dPhi_physics->GetMaximum())*0.75;
  line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
  gPad->Print("dPhi_physics_pipi0eta.png");
  
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
    std::cout<<"Data: ";
    double physics_beforeSignal=h_MBC_physics->Integral(h_MBC_physics->GetXaxis()->GetFirst(), h_MBC_physics->GetXaxis()->FindBin(mbcCut_center-mbcCut_range));
    double physics_afterSignal=h_MBC_physics->Integral(h_MBC_physics->GetXaxis()->FindBin(mbcCut_center+mbcCut_range), h_MBC_physics->GetXaxis()->GetLast());
    std::cout<<" before signal = "<<physics_beforeSignal<<" MC/Data = "<<(generic_beforeSignal/20+continu_beforeSignal/5)/physics_beforeSignal;
    std::cout<<", after signal = "<<physics_afterSignal<<" MC/Data = "<<(generic_afterSignal/20+continu_afterSignal/5)/physics_afterSignal<<std::endl;
    std::cout<<"===================="<<std::endl;
    std::cout<<" TeX output "<<std::endl;
    std::cout.precision(2);
    std::cout<<" & "<<generic_beforeSignal/20<<" $\\pm$ "<<pow(generic_beforeSignal, 0.5)/20;
    std::cout<<" & "<<continu_beforeSignal/5<<" $\\pm$ "<<pow(continu_beforeSignal, 0.5)/5;
    double mc_beforeSignal=generic_beforeSignal/20+continu_beforeSignal/5;
    double mc_beforeSignal_error=pow(generic_beforeSignal/400+continu_beforeSignal/25, 0.5);
    std::cout<<" & "<<mc_beforeSignal<<" $\\pm$ "<<mc_beforeSignal_error;
    std::cout<<" & "<<physics_beforeSignal<<" $\\pm$ "<<pow(physics_beforeSignal, 0.5);
    double dataMC_beforeSignal=physics_beforeSignal/mc_beforeSignal;
    double dataMC_beforeSignal_error=dataMC_beforeSignal*pow(1/physics_beforeSignal+pow(mc_beforeSignal_error/mc_beforeSignal, 2), 0.5);
    std::cout<<" & "<<dataMC_beforeSignal<<" $\\pm$ "<<dataMC_beforeSignal_error;
    std::cout<<" & "<<generic_afterSignal/20<<" $\\pm$ "<<pow(generic_afterSignal, 0.5)/20;
    std::cout<<" & "<<continu_afterSignal/5<<" $\\pm$ "<<pow(continu_afterSignal, 0.5)/5;
    double mc_afterSignal=generic_afterSignal/20+continu_afterSignal/5;
    double mc_afterSignal_error=pow(generic_afterSignal/400+continu_afterSignal/25, 0.5);
    std::cout<<" & "<<mc_afterSignal<<" $\\pm$ "<<mc_afterSignal_error;
    std::cout<<" & "<<physics_afterSignal<<" $\\pm$ "<<pow(physics_afterSignal, 0.5);
    double dataMC_afterSignal=physics_afterSignal/mc_afterSignal;
    double dataMC_afterSignal_error=dataMC_afterSignal*pow(1/physics_afterSignal+pow(mc_afterSignal_error/mc_afterSignal, 2), 0.5);
    std::cout<<" & "<<dataMC_afterSignal<<" $\\pm$ "<<dataMC_afterSignal_error<<" \\\\  \\hline"<<std::endl;
  }
  
  if (DeltaM_sideband)
  {
    double full=0.5;
    double half=0.3;
    std::cout<<"=== DeltaM sideband ==="<<std::endl;
    std::cout<<"Generic MC: ";
    double generic_full=h_DeltaM_generic->Integral(h_DeltaM_generic->GetXaxis()->FindBin(deltaMCut_center+deltaMCut_range), h_DeltaM_generic->GetXaxis()->FindBin(full));
    double generic_half=h_DeltaM_generic->Integral(h_DeltaM_generic->GetXaxis()->FindBin(deltaMCut_center+deltaMCut_range), h_DeltaM_generic->GetXaxis()->FindBin(half));
    std::cout<<" between 0.0 & "<<full<<" /20 = "<<generic_full/20;
    std::cout<<", between 0.0 & "<<half<<" /20 = "<<generic_half/20<<std::endl;
    std::cout<<"Continuum MC: ";
    double continu_full=h_DeltaM_continu->Integral(h_DeltaM_continu->GetXaxis()->FindBin(deltaMCut_center+deltaMCut_range), h_DeltaM_continu->GetXaxis()->FindBin(full));
    double continu_half=h_DeltaM_continu->Integral(h_DeltaM_continu->GetXaxis()->FindBin(deltaMCut_center+deltaMCut_range), h_DeltaM_continu->GetXaxis()->FindBin(half));
    std::cout<<" between 0.0 & "<<full<<" /5 = "<<continu_full/5;
    std::cout<<", between 0.0 & "<<half<<" /5 = "<<continu_half/5<<std::endl;
    std::cout<<"Sum MC: ";
    std::cout<<" between 0.0 & "<<full<<" = "<<generic_full/20+continu_full/5;
    std::cout<<", between 0.0 & "<<half<<" = "<<generic_half/20+continu_half/5<<std::endl;
    std::cout<<"Data: ";
    double physics_full=h_DeltaM_physics->Integral(h_DeltaM_physics->GetXaxis()->FindBin(deltaMCut_center+deltaMCut_range), h_DeltaM_physics->GetXaxis()->FindBin(full));
    double physics_half=h_DeltaM_physics->Integral(h_DeltaM_physics->GetXaxis()->FindBin(deltaMCut_center+deltaMCut_range), h_DeltaM_physics->GetXaxis()->FindBin(half));
    std::cout<<" between 0.0 & "<<full<<" = "<<physics_full;
    std::cout<<", between 0.0 & "<<half<<" = "<<physics_half<<std::endl;
    std::cout<<" MC/data between 0.0 & "<<full<<" = "<<(generic_full/20+continu_full/5)/physics_full;
    std::cout<<", MC/data beteen 0.0 & "<<half<<" = "<<(generic_half/20+continu_half/5)/physics_half<<std::endl;
    std::cout<<"===================="<<std::endl;
    std::cout<<" TeX output "<<std::endl;
    std::cout.precision(3);
    std::cout<<" & "<<generic_full/20<<" $\\pm$ "<<pow(generic_full, 0.5)/20;
    std::cout<<" & "<<continu_full/5<<" $\\pm$ "<<pow(continu_full, 0.5)/5;
    double mc_full=generic_full/20+continu_full/5;
    double mc_full_error=pow(generic_full/400+continu_full/25, 0.5);
    std::cout<<" & "<<mc_full<<" $\\pm$ "<<mc_full_error;
    std::cout<<" & "<<physics_full<<" $\\pm$ "<<pow(physics_full, 0.5);
    double dataMC_full=physics_full/mc_full;
    double dataMC_full_error=dataMC_full*pow(1/physics_full+pow(mc_full_error/mc_full, 2), 0.5);
    std::cout<<" & "<<dataMC_full<<" $\\pm$ "<<dataMC_full_error;
    std::cout<<" & "<<generic_half/20<<" $\\pm$ "<<pow(generic_half, 0.5)/20;
    std::cout<<" & "<<continu_half/5<<" $\\pm$ "<<pow(continu_half, 0.5)/5;
    double mc_half=generic_half/20+continu_half/5;
    double mc_half_error=pow(generic_half/400+continu_half/25, 0.5);
    std::cout<<" & "<<mc_half<<" $\\pm$ "<<mc_half_error;
    std::cout<<" & "<<physics_half<<" $\\pm$ "<<pow(physics_half, 0.5);
    double dataMC_half=physics_half/mc_half;
    double dataMC_half_error=dataMC_half*pow(1/physics_half+pow(mc_half_error/mc_half, 2), 0.5);
    std::cout<<" & "<<dataMC_half<<" $\\pm$ "<<dataMC_half_error<<" \\\\  \\hline"<<std::endl;
  }
  
  if (!mBC_sideband && !DeltaM_sideband)
  {
    double generic_error=pow(genericEvents.dPhiCut, 0.5)/20;
    double continu_error=pow(continuEvents.dPhiCut, 0.5)/5;
    double physics_error=pow(physicsEvents.dPhiCut, 0.5);
    double mc=genericEvents.dPhiCut/20+continuEvents.dPhiCut/5;
    double mc_error=pow(genericEvents.dPhiCut/400+continuEvents.dPhiCut/25, 0.5);
    
    std::cout<<"=== Signal Region ==="<<std::endl;
    std::cout<<" & "<<genericEvents.dPhiCut/20<<" $\\pm$ "<<generic_error;
    std::cout<<" & "<<continuEvents.dPhiCut/5<<" $\\pm$ "<<continu_error;
    std::cout<<" & "<<mc<<" $\\pm$ "<<mc_error;
    std::cout<<" & "<<physicsEvents.dPhiCut<<" $\\pm$ "<<physics_error<<std::endl;
  }
  
  return 0;
}
