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

double pi=3.14159265358979;
double dsCharge=+1;
double dsPlusMCut_center=1.96849, dsPlusMCut_range=0.011;
double deltaECut_center=0.012, deltaECut_range=0.0; // 0.03
double mbcCut_center=2.112, mbcCut_range=0.005; // 0.01
double deltaMCut_center=0.158, deltaMCut_range=0.009; // 0.013
double diffD0Cut=-0.002; // -0.005
double dPhiCutLess=0.06; // 0.15
/*
double deltaECut_center=0.012, deltaECut_range=0.03; // 0.03
double mbcCut_center=2.112, mbcCut_range=0.01; // 0.01
double deltaMCut_center=0.158, deltaMCut_range=0.013; // 0.013
double diffD0Cut=-0.005; // -0.005
double dPhiCutLess=0.15; // 0.15
*/
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

int DsTaggedAnalysis_pionFit_KKpi()
{
  TFile *signalFile;
  if (dsCharge>0) signalFile = new TFile("/nfs/cor/an2/souvik/MC_vtosll_Dsp_KKpi/DsTaggedDecaysProc_MC_vtosll_Dsp_KKpi.root");
  else signalFile = new TFile("/nfs/cor/an2/souvik/MC_vtosll_Dsm_KKpi/DsTaggedDecaysProc_MC_vtosll_Dsm_KKpi.root");
  signalFile->cd("DsTaggedDecaysProc");
  TTree *signalTree = (TTree*)gDirectory->Get("nt1");
  
  TFile *converFile;
  if (dsCharge>0) converFile = new TFile("/nfs/cor/an2/souvik/MC_gamma_Dsp_KKpi/DsTaggedDecaysProc_MC_gamma_Dsp_KKpi.root");
  else converFile = new TFile("/nfs/cor/an2/souvik/MC_gamma_Dsm_KKpi/DsTaggedDecaysProc_MC_gamma_Dsm_KKpi.root");
  converFile->cd("DsTaggedDecaysProc");
  TTree *converTree = (TTree*)gDirectory->Get("nt1");
  
  TFile *physicFile = new TFile("/nfs/cor/an3/souvik/Dataset39/DsTaggedProc_Data_213586_214863");
  physicFile->cd("DsTaggedDecaysProc");
  TTree *physicTree = (TTree*)gDirectory->Get("nt1");
  
  NEvents signalEvents={0,0,0,0,0,0,0,0};
  NEvents converEvents={0,0,0,0,0,0,0,0};
  NEvents physicEvents={0,0,0,0,0,0,0,0};
  
  EventNumber signalNumber={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  EventNumber converNumber={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  EventNumber physicNumber={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  float runNumber_signal, eventNumber_signal;
  float runNumber_conver, eventNumber_conver;
  float runNumber_physic, eventNumber_physic;
  float dsPlusM_signal, dsPlusCharge_signal, DeltaE_signal, MBC_signal, DeltaM_signal, eeMass_signal, decayMode_signal;
  float dsPlusM_conver, dsPlusCharge_conver, DeltaE_conver, MBC_conver, DeltaM_conver, eeMass_conver, decayMode_conver;
  float dsPlusM_physic, dsPlusCharge_physic, DeltaE_physic, MBC_physic, DeltaM_physic, eeMass_physic, decayMode_physic;
  float d0_e_signal, d0_p_signal, z0_e_signal, z0_p_signal, px_e_signal, py_e_signal, pz_e_signal, px_p_signal, py_p_signal, pz_p_signal, E_e_signal, E_p_signal, curv_e_signal, curv_p_signal;
  float d0_e_conver, d0_p_conver, z0_e_conver, z0_p_conver, px_e_conver, py_e_conver, pz_e_conver, px_p_conver, py_p_conver, pz_p_conver, E_e_conver, E_p_conver, curv_e_conver, curv_p_conver;
  float d0_e_physic, d0_p_physic, z0_e_physic, z0_p_physic, px_e_physic, py_e_physic, pz_e_physic, px_p_physic, py_p_physic, pz_p_physic, E_e_physic, E_p_physic, curv_e_physic, curv_p_physic;
  float px_e_signal_MC, py_e_signal_MC, pz_e_signal_MC, E_e_signal_MC, px_p_signal_MC, py_p_signal_MC, pz_p_signal_MC, E_p_signal_MC;
  float px_e_conver_MC, py_e_conver_MC, pz_e_conver_MC, E_e_conver_MC, px_p_conver_MC, py_p_conver_MC, pz_p_conver_MC, E_p_conver_MC;
  float theta1_signal, theta2_signal;
  float theta1_conver, theta2_conver;
  
  typedef std::map<int, float> DecayMap;
  DecayMap decayFrequency_signal, decayFrequency_conver, decayFrequency_physic;
  
  TH1D *h_dsPlusM_signal = new TH1D("h_dsPlusM_signal", "m_{D_{S}^{+}} Signal Sample; GeV", 100, 1.9, 2.1); h_dsPlusM_signal->SetLineColor(kRed);
  TH1D *h_dsPlusM_conver = new TH1D("h_dsPlusM_conver", "m_{D_{S}^{+}} Conversion Background Sample; GeV", 100, 1.9, 2.1); h_dsPlusM_conver->SetLineColor(kBlue);
  TH1D *h_dsPlusM_physic = new TH1D("h_dsPlusM_physic", "m_{D_{S}^{+}} Data; GeV", 100, 1.9, 2.1); h_dsPlusM_physic->SetLineColor(kGreen);
  TH1D *h_DeltaE_signal = new TH1D("h_DeltaE_signal", "#DeltaE Signal Sample; GeV", 100, -0.1, 0.2); h_DeltaE_signal->SetLineColor(kRed);
  TH1D *h_DeltaE_conver = new TH1D("h_DeltaE_conver", "#DeltaE Conversion Background Sample; GeV", 100, -0.1, 0.2); h_DeltaE_conver->SetLineColor(kBlue);
  TH1D *h_DeltaE_physic = new TH1D("h_DeltaE_physic", "#DeltaE Data; GeV", 100, -0.1, 0.2); h_DeltaE_physic->SetLineColor(kGreen);
  TH1D *h_MBC_signal = new TH1D("h_MBC_signal", "m_{BC} Signal Sample; GeV", 100, 2., 2.2); h_MBC_signal->SetLineColor(kRed);
  TH1D *h_MBC_conver = new TH1D("h_MBC_conver", "m_{BC} Conversion Background; GeV", 100, 2., 2.2); h_MBC_conver->SetLineColor(kBlue);
  TH1D *h_MBC_physic = new TH1D("h_MBC_physic", "m_{BC} Data; GeV", 100, 2., 2.2); h_MBC_physic->SetLineColor(kGreen);
  TH1D *h_DeltaM_signal = new TH1D("h_DeltaM_signal", "#deltaM Signal Sample", 100, 0.0, 0.2); h_DeltaM_signal->SetLineColor(kRed);
  TH1D *h_DeltaM_conver = new TH1D("h_DeltaM_conver", "#deltaM Conversion Background Sample; GeV", 100, 0.0, 0.2); h_DeltaM_conver->SetLineColor(kBlue);
  TH1D *h_DeltaM_physic = new TH1D("h_DeltaM_physic", "#deltaM Data; GeV", 100, 0.0, 0.2); h_DeltaM_physic->SetLineColor(kGreen);
  
  TH2D *h_DeltaE_MBC_signal = new TH2D("h_DeltaE_MBC_signal", "h_DeltaE_MBC_signal", 100, -0.1, 0.2, 100, 2., 2.2);
  TH2D *h_DeltaE_MBC_conver = new TH2D("h_DeltaE_MBC_conver", "h_DeltaE_MBC_conver", 100, -0.1, 0.2, 100, 2., 2.2);
  TH2D *h_DeltaE_MBC_physic = new TH2D("h_DeltaE_MBC_physic", "h_DeltaE_MBC_physic", 100, -0.1, 0.2, 100, 2., 2.2);
  TH2D *h_DeltaE_DeltaM_signal = new TH2D("h_DeltaE_DeltaM_signal", "h_DeltaE_DeltaM_signal", 100, -0.1, 0.2, 100, 0.0, 0.2);
  TH2D *h_DeltaE_DeltaM_conver = new TH2D("h_DeltaE_DeltaM_conver", "h_DeltaE_DeltaM_conver", 100, -0.1, 0.2, 100, 0.0, 0.2);
  TH2D *h_DeltaE_DeltaM_physic = new TH2D("h_DeltaE_DeltaM_physic", "h_DeltaE_DeltaM_physic", 100, -0.1, 0.2, 100, 0.0, 0.2);
  TH2D *h_MBC_DeltaM_signal = new TH2D("h_MBC_DeltaM_signal", "h_MBC_DeltaM_signal", 100, 2., 2.2, 100, 0.0, 0.2);
  TH2D *h_MBC_DeltaM_conver = new TH2D("h_MBC_DeltaM_conver", "h_MBC_DeltaM_conver", 100, 2., 2.2, 100, 0.0, 0.2);
  TH2D *h_MBC_DeltaM_physic = new TH2D("h_MBC_DeltaM_physic", "h_MBC_DeltaM_physic", 100, 2., 2.2, 100, 0.0, 0.2);
  
  TH1D *h_d0_e_signal = new TH1D("h_d0_e_signal", "d0_e", 50, -0.01, 0.01); h_d0_e_signal->SetLineColor(kRed);
  TH1D *h_d0_e_conver = new TH1D("h_d0_e_conver", "d0_e", 50, -0.01, 0.01); h_d0_e_conver->SetLineColor(kBlue);
  TH1D *h_d0_e_physic = new TH1D("h_d0_e_physic", "d0_e", 50, -0.01, 0.01); h_d0_e_physic->SetLineColor(kGreen);
  TH1D *h_d0_p_signal = new TH1D("h_d0_p_signal", "d0_p", 50, -0.01, 0.01); h_d0_p_signal->SetLineColor(kRed);
  TH1D *h_d0_p_conver = new TH1D("h_d0_p_conver", "d0_p", 50, -0.01, 0.01); h_d0_p_conver->SetLineColor(kBlue);
  TH1D *h_d0_p_physic = new TH1D("h_d0_p_physic", "d0_p", 50, -0.01, 0.01); h_d0_p_physic->SetLineColor(kGreen);
  TH1D *h_diffD0_signal = new TH1D("h_diffD0_signal", "#Deltad_{0} Signal Sample; m", 50, -0.01, 0.01); h_diffD0_signal->SetLineColor(kRed);
  TH1D *h_diffD0_conver = new TH1D("h_diffD0_conver", "#Deltad_{0} Conversion Background Sample; m", 50, -0.01, 0.01); h_diffD0_conver->SetLineColor(kBlue);
  TH1D *h_diffD0_physic = new TH1D("h_diffD0_physic", "#Deltad_{0} Data; m", 50, -0.01, 0.01); h_diffD0_physic->SetLineColor(kGreen);
  TH1D *h_dPhi_signal = new TH1D("h_dPhi_signal", "#Delta#Phi Signal Sample", 50, -2., 2.); h_dPhi_signal->SetLineColor(kRed);
  TH1D *h_dPhi_conver = new TH1D("h_dPhi_conver", "#Delta#Phi Conversion Background Sample", 50, -2., 2.); h_dPhi_conver->SetLineColor(kBlue);
  TH1D *h_dPhi_physic = new TH1D("h_dPhi_physic", "#Delta#Phi Data", 50, -2., 2.); h_dPhi_physic->SetLineColor(kGreen);
  TH2D *h_dPhi_ee_signal = new TH2D("h_dPhi_ee_signal", "dPhi_ee_signal", 50, 0., 0.2, 50, -2., 2.);
  TH2D *h_dPhi_ee_conver = new TH2D("h_dPhi_ee_conver", "dPhi_ee_conver", 50, 0., 0.2, 50, -2., 2.);
  TH2D *h_dPhi_ee_physic = new TH2D("h_dPhi_ee_physic", "dPhi_ee_physic", 50, 0., 0.2, 50, -2., 2.);
  TH2D *h_d0_phi_e_signal = new TH2D("h_d0_phi_e_signal", "h_d0_phi_e_signal", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_p_signal = new TH2D("h_d0_phi_p_signal", "h_d0_phi_p_signal", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_e_conver = new TH2D("h_d0_phi_e_conver", "h_d0_phi_e_conver", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_p_conver = new TH2D("h_d0_phi_p_conver", "h_d0_phi_p_conver", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_e_physic = new TH2D("h_d0_phi_e_physic", "h_d0_phi_e_physic", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_p_physic = new TH2D("h_d0_phi_p_physic", "h_d0_phi_p_physic", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_dPhi_diffD0_signal = new TH2D("h_dPhi_diffD0_signal", "#Delta#Phi vs #Deltad_{0} Signal Sample; ; m", 50, -2., 2., 50, -0.01, 0.01); h_dPhi_diffD0_signal->SetLineColor(kRed);
  TH2D *h_dPhi_diffD0_conver = new TH2D("h_dPhi_diffD0_conver", "#Delta#Phi vs #Deltad_{0} Conversion Background Sample; ; m", 50, -2., 2., 50, -0.01, 0.01); h_dPhi_diffD0_conver->SetLineColor(kBlue);
  TH2D *h_dPhi_diffD0_physic = new TH2D("h_dPhi_diffD0_physic", "#Delta#Phi vs #Deltad_{0} Data; ; m", 50, -2., 2., 50, -0.01, 0.01); h_dPhi_diffD0_physic->SetLineColor(kGreen);
  TH1D *h_R_signal = new TH1D("h_R_signal", "h_R_signal", 50, -.4, .4); h_R_signal->SetLineColor(kRed);
  TH1D *h_R_conver = new TH1D("h_R_conver", "h_R_conver", 50, -.4, .4); h_R_conver->SetLineColor(kBlue);
  TH1D *h_R_physic = new TH1D("h_R_physic", "h_R_physic", 50, -.4, .4); h_R_physic->SetLineColor(kGreen);
  TH1D *h_ee_signal = new TH1D("h_ee_signal", "h_ee_signal", 20, 0.0, 0.2); h_ee_signal->SetLineColor(kRed);
  TH1D *h_ee_conver = new TH1D("h_ee_conver", "h_ee_conver", 20, 0.0, 0.2); h_ee_conver->SetLineColor(kBlue);
  TH1D *h_ee_physic = new TH1D("h_ee_physic", "h_ee_physic", 20, 0.0, 0.2); h_ee_physic->SetLineColor(kGreen);
  TH2D *h_electronE = new TH2D("h_electronE", "Electron Energy Resolution; Reco Electron Energy (GeV)", 100, 0.0, 0.16, 100, -0.03, 0.03);
  TH2D *h_electronPx = new TH2D("h_electronPx", "Electron Px Resolution; Reco Electron Px (GeV)", 100, 0.0, 0.14, 100, -0.03, 0.03); 
  TH2D *h_electronPy = new TH2D("h_electronPy", "Electron Py Resolution; Reco Electron Py (GeV)", 100, 0.0, 0.14, 100, -0.03, 0.03); 
  TH2D *h_electronPz = new TH2D("h_electronPz", "Electron Pz Resolution; Reco Electron Pz (GeV)", 100, 0.0, 0.14, 100, -0.03, 0.03); 
  TH1D *h_theta1_signal = new TH1D("h_theta1_signal", "h_theta1_signal", 100, 0., pi/2);
  TH1D *h_theta1_conver = new TH1D("h_theta1_conver", "h_theta1_conver", 100, 0., pi/2);
  TH1D *h_theta2_signal = new TH1D("h_theta2_signal", "h_theta2_signal", 100, 0., pi/2);
  TH1D *h_theta2_conver = new TH1D("h_theta2_conver", "h_theta2_conver", 100, 0., pi/2);
  
  signalTree->SetBranchAddress("Run", &(runNumber_signal));
  signalTree->SetBranchAddress("Event", &(eventNumber_signal));
  signalTree->SetBranchAddress("dsPlusM", &(dsPlusM_signal));
  signalTree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_signal));
  signalTree->SetBranchAddress("DecayMode", &(decayMode_signal));
  signalTree->SetBranchAddress("DeltaE", &(DeltaE_signal));
  signalTree->SetBranchAddress("MBC", &(MBC_signal));
  signalTree->SetBranchAddress("DeltaM", &(DeltaM_signal));
  signalTree->SetBranchAddress("kElectron1E_MC", &(E_e_signal_MC));  
  signalTree->SetBranchAddress("kElectron1Px_MC", &(px_e_signal_MC));
  signalTree->SetBranchAddress("kElectron1Py_MC", &(py_e_signal_MC));
  signalTree->SetBranchAddress("kElectron1Pz_MC", &(pz_e_signal_MC));
  signalTree->SetBranchAddress("kElectron2E_MC", &(E_p_signal_MC));
  signalTree->SetBranchAddress("kElectron2Px_MC", &(px_p_signal_MC));
  signalTree->SetBranchAddress("kElectron2Py_MC", &(py_p_signal_MC));
  signalTree->SetBranchAddress("kElectron2Pz_MC", &(pz_p_signal_MC));
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
  signalTree->SetBranchAddress("keeMass_reco", &(eeMass_signal));
  signalTree->SetBranchAddress("kTheta1", &(theta1_signal));
  signalTree->SetBranchAddress("kTheta2", &(theta2_signal));
  
  converTree->SetBranchAddress("Run", &(runNumber_conver));
  converTree->SetBranchAddress("Event", &(eventNumber_conver));
  converTree->SetBranchAddress("dsPlusM", &(dsPlusM_conver));
  converTree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_conver));
  converTree->SetBranchAddress("DecayMode", &(decayMode_conver));
  converTree->SetBranchAddress("DeltaE", &(DeltaE_conver));
  converTree->SetBranchAddress("MBC", &(MBC_conver));
  converTree->SetBranchAddress("DeltaM", &(DeltaM_conver));
  converTree->SetBranchAddress("kElectron1E_MC", &(E_e_conver_MC));  
  converTree->SetBranchAddress("kElectron1Px_MC", &(px_e_conver_MC));
  converTree->SetBranchAddress("kElectron1Py_MC", &(py_e_conver_MC));
  converTree->SetBranchAddress("kElectron1Pz_MC", &(pz_e_conver_MC));
  converTree->SetBranchAddress("kElectron2E_MC", &(E_p_conver_MC));
  converTree->SetBranchAddress("kElectron2Px_MC", &(px_p_conver_MC));
  converTree->SetBranchAddress("kElectron2Py_MC", &(py_p_conver_MC));
  converTree->SetBranchAddress("kElectron2Pz_MC", &(pz_p_conver_MC));
  converTree->SetBranchAddress("kElectron1D0_reco", &(d0_e_conver));
  converTree->SetBranchAddress("kElectron2D0_reco", &(d0_p_conver));
  converTree->SetBranchAddress("kElectron1Z0_reco", &(z0_e_conver));
  converTree->SetBranchAddress("kElectron2Z0_reco", &(z0_p_conver));
  converTree->SetBranchAddress("kElectron1Px_reco", &(px_e_conver));
  converTree->SetBranchAddress("kElectron1Py_reco", &(py_e_conver));
  converTree->SetBranchAddress("kElectron1Pz_reco", &(pz_e_conver));
  converTree->SetBranchAddress("kElectron2Px_reco", &(px_p_conver));
  converTree->SetBranchAddress("kElectron2Py_reco", &(py_p_conver));
  converTree->SetBranchAddress("kElectron2Pz_reco", &(pz_p_conver));
  converTree->SetBranchAddress("kElectron1E_reco", &(E_e_conver));
  converTree->SetBranchAddress("kElectron2E_reco", &(E_p_conver));
  converTree->SetBranchAddress("kElectron1Curv_reco", &(curv_e_conver));
  converTree->SetBranchAddress("kElectron2Curv_reco", &(curv_p_conver));
  converTree->SetBranchAddress("keeMass_reco", &(eeMass_conver));
  converTree->SetBranchAddress("kTheta1", &(theta1_conver));
  converTree->SetBranchAddress("kTheta2", &(theta2_conver));
  
  physicTree->SetBranchAddress("Run", &(runNumber_physic));
  physicTree->SetBranchAddress("Event", &(eventNumber_physic));
  physicTree->SetBranchAddress("dsPlusM", &(dsPlusM_physic));
  physicTree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_physic));
  physicTree->SetBranchAddress("DecayMode", &(decayMode_physic));
  physicTree->SetBranchAddress("DeltaE", &(DeltaE_physic));
  physicTree->SetBranchAddress("MBC", &(MBC_physic));
  physicTree->SetBranchAddress("DeltaM", &(DeltaM_physic));
  physicTree->SetBranchAddress("kElectron1D0_reco", &(d0_e_physic));
  physicTree->SetBranchAddress("kElectron2D0_reco", &(d0_p_physic));
  physicTree->SetBranchAddress("kElectron1Z0_reco", &(z0_e_physic));
  physicTree->SetBranchAddress("kElectron2Z0_reco", &(z0_p_physic));
  physicTree->SetBranchAddress("kElectron1Px_reco", &(px_e_physic));
  physicTree->SetBranchAddress("kElectron1Py_reco", &(py_e_physic));
  physicTree->SetBranchAddress("kElectron1Pz_reco", &(pz_e_physic));
  physicTree->SetBranchAddress("kElectron2Px_reco", &(px_p_physic));
  physicTree->SetBranchAddress("kElectron2Py_reco", &(py_p_physic));
  physicTree->SetBranchAddress("kElectron2Pz_reco", &(pz_p_physic));
  physicTree->SetBranchAddress("kElectron1E_reco", &(E_e_physic));
  physicTree->SetBranchAddress("kElectron2E_reco", &(E_p_physic));
  physicTree->SetBranchAddress("kElectron1Curv_reco", &(curv_e_physic));
  physicTree->SetBranchAddress("kElectron2Curv_reco", &(curv_p_physic));
  physicTree->SetBranchAddress("keeMass_reco", &(eeMass_physic));
  
  int nSignalEvents=signalTree->GetEntries();
  int oppRecon_signal=0;
  for (int i=0; i<nSignalEvents; ++i)
  {
    signalTree->GetEvent(i);
    
    double phi_e=atan2(py_e_signal, px_e_signal);
    double phi_p=atan2(py_p_signal, px_p_signal);
    double dPhi=deltaPhi(phi_e, phi_p);
    
    if (runNumber_signal!=signalNumber.noCut_run || eventNumber_signal!=signalNumber.noCut_event)
    {      
      if (decayMode_signal>-1) decayFrequency_signal[int(decayMode_signal)]+=1;
      ++signalEvents.noCut;
      signalNumber.noCut_run=runNumber_signal;
      signalNumber.noCut_event=eventNumber_signal;
    }
    
    if (fabs(d0_e_signal)<0.005 && fabs(d0_p_signal)<0.005 &&
        fabs(z0_e_signal)<0.05 && fabs(z0_p_signal)<0.05)
    {
    
    if (decayMode_signal==401)
    {
      if (runNumber_signal!=signalNumber.tagCut_run || eventNumber_signal!=signalNumber.tagCut_event)
      {      
        ++signalEvents.tagCut;
        signalNumber.tagCut_run=runNumber_signal;
        signalNumber.tagCut_event=eventNumber_signal;
      }
      h_dsPlusM_signal->Fill(dsPlusM_signal);
      if (dsPlusMCut(dsPlusM_signal))
      {
        if (runNumber_signal!=signalNumber.dsPlusMCut_run || eventNumber_signal!=signalNumber.dsPlusMCut_event)
        {
          ++signalEvents.dsPlusMCut;
          signalNumber.dsPlusMCut_run=runNumber_signal;
          signalNumber.dsPlusMCut_event=eventNumber_signal;
        }
        h_DeltaE_signal->Fill(DeltaE_signal);
 
        // Optimization plots
        h_DeltaE_MBC_signal->Fill(DeltaE_signal, MBC_signal);
        h_DeltaE_DeltaM_signal->Fill(DeltaE_signal, DeltaM_signal);
        h_MBC_DeltaM_signal->Fill(MBC_signal, DeltaM_signal);
        //if (DeltaECut(DeltaE_signal))
        {
          if (runNumber_signal!=signalNumber.deltaECut_run || eventNumber_signal!=signalNumber.deltaECut_event)
          {
            ++signalEvents.deltaECut;
            signalNumber.deltaECut_run=runNumber_signal;
            signalNumber.deltaECut_event=eventNumber_signal;
          }
          h_MBC_signal->Fill(MBC_signal);
          if (MBCCut(MBC_signal))
          {
            if (runNumber_signal!=signalNumber.mbcCut_run || eventNumber_signal!=signalNumber.mbcCut_event)
            {
              ++signalEvents.mbcCut;
              signalNumber.mbcCut_run=runNumber_signal;
              signalNumber.mbcCut_event=eventNumber_signal;
            }
            h_DeltaM_signal->Fill(DeltaM_signal);
            if (DeltaMCut(DeltaM_signal))
            {
              if (runNumber_signal!=signalNumber.deltaMCut_run || eventNumber_signal!=signalNumber.deltaMCut_event) 
              {
                ++signalEvents.deltaMCut;
                signalNumber.deltaMCut_run=runNumber_signal;
                signalNumber.deltaMCut_event=eventNumber_signal;
              }
              h_d0_e_signal->Fill(d0_e_signal);
              h_d0_p_signal->Fill(d0_p_signal);
              h_diffD0_signal->Fill(d0_e_signal-d0_p_signal);
              h_dPhi_diffD0_signal->Fill(dPhi, d0_e_signal-d0_p_signal);
              h_d0_phi_e_signal->Fill(phi_e, d0_e_signal);
              h_d0_phi_p_signal->Fill(phi_p, d0_p_signal);
              if (dD0(d0_e_signal-d0_p_signal))
              {
                if (runNumber_signal!=signalNumber.diffD0Cut_run || eventNumber_signal!=signalNumber.diffD0Cut_event) 
                {
                  ++signalEvents.diffD0Cut;
                  signalNumber.diffD0Cut_run=runNumber_signal;
                  signalNumber.diffD0Cut_event=eventNumber_signal;
                }
                h_dPhi_signal->Fill(dPhi);
                h_dPhi_ee_signal->Fill(eeMass_signal, dPhi);
                if (dPhiCut(dPhi))
                {
                  if (runNumber_signal!=signalNumber.dPhiCut_run || eventNumber_signal!=signalNumber.dPhiCut_event) 
                  {
                    ++signalEvents.dPhiCut;
                    signalNumber.dPhiCut_run=runNumber_signal;
                    signalNumber.dPhiCut_event=eventNumber_signal;                   
                    h_electronE->Fill(E_e_signal, E_e_signal-E_e_signal_MC);
                    h_electronE->Fill(E_p_signal, E_p_signal-E_p_signal_MC);
                    h_electronPx->Fill(px_e_signal, px_e_signal-px_e_signal_MC);
                    h_electronPx->Fill(px_p_signal, px_p_signal-px_p_signal_MC);
                    h_electronPy->Fill(py_e_signal, py_e_signal-py_e_signal_MC);
                    h_electronPy->Fill(py_p_signal, py_p_signal-py_p_signal_MC);
                    h_electronPz->Fill(pz_e_signal, pz_e_signal-pz_e_signal_MC);
                    h_electronPz->Fill(pz_p_signal, pz_p_signal-pz_p_signal_MC);
                    h_theta1_signal->Fill(theta1_signal);
                    h_theta2_signal->Fill(theta2_signal);
                  }
                  h_ee_signal->Fill(eeMass_signal);
                  double r1=-1/(2*curv_e_signal);
                  double r2=1/(2*curv_p_signal);
                  double R=convert(r1, r2, d0_e_signal, d0_p_signal, phi_e, phi_p);
                  h_R_signal->Fill(R);
                  if ((dsCharge * dsPlusCharge_signal)<0) oppRecon_signal+=1;
                  
                }
              }
            }
          }
        }
      }
    }
    }
  }
  
  
  int nConverEvents=converTree->GetEntries();
  int oppRecon_conver=0;
  for (int i=0; i<nConverEvents; ++i)
  {
    converTree->GetEvent(i);
    
    double phi_e=atan2(py_e_conver, px_e_conver);
    double phi_p=atan2(py_p_conver, px_p_conver);
    double dPhi=deltaPhi(phi_e, phi_p);
    
    if (runNumber_conver!=converNumber.noCut_run || eventNumber_conver!=converNumber.noCut_event)
    {      
      if (decayMode_conver>-1) decayFrequency_conver[int(decayMode_conver)]+=1;
      ++converEvents.noCut;
      converNumber.noCut_run=runNumber_conver;
      converNumber.noCut_event=eventNumber_conver;
    }
    
    if (fabs(d0_e_conver)<0.005 && fabs(d0_p_conver)<0.005 &&
        fabs(z0_e_conver)<0.05 && fabs(z0_p_conver)<0.05)
    {
    
    if (decayMode_conver==401)
    {
      if (runNumber_conver!=converNumber.tagCut_run || eventNumber_conver!=converNumber.tagCut_event)
      {      
        ++converEvents.tagCut;
        converNumber.tagCut_run=runNumber_conver;
        converNumber.tagCut_event=eventNumber_conver;
      }
      h_dsPlusM_conver->Fill(dsPlusM_conver);
      if (dsPlusMCut(dsPlusM_conver))
      {
        if (runNumber_conver!=converNumber.dsPlusMCut_run || eventNumber_conver!=converNumber.dsPlusMCut_event)
        {
          ++converEvents.dsPlusMCut;
          converNumber.dsPlusMCut_run=runNumber_conver;
          converNumber.dsPlusMCut_event=eventNumber_conver;
        }
        h_DeltaE_conver->Fill(DeltaE_conver);
 
        // Optimization plots
        h_DeltaE_MBC_conver->Fill(DeltaE_conver, MBC_conver);
        h_DeltaE_DeltaM_conver->Fill(DeltaE_conver, DeltaM_conver);
        h_MBC_DeltaM_conver->Fill(MBC_conver, DeltaM_conver);
        //if (DeltaECut(DeltaE_conver))
        {
          if (runNumber_conver!=converNumber.deltaECut_run || eventNumber_conver!=converNumber.deltaECut_event)
          {
            ++converEvents.deltaECut;
            converNumber.deltaECut_run=runNumber_conver;
            converNumber.deltaECut_event=eventNumber_conver;
          }
          h_MBC_conver->Fill(MBC_conver);
          if (MBCCut(MBC_conver))
          {
            if (runNumber_conver!=converNumber.mbcCut_run || eventNumber_conver!=converNumber.mbcCut_event)
            {
              ++converEvents.mbcCut;
              converNumber.mbcCut_run=runNumber_conver;
              converNumber.mbcCut_event=eventNumber_conver;
            }
            h_DeltaM_conver->Fill(DeltaM_conver);
            if (DeltaMCut(DeltaM_conver))
            {
              if (runNumber_conver!=converNumber.deltaMCut_run || eventNumber_conver!=converNumber.deltaMCut_event) 
              {
                ++converEvents.deltaMCut;
                converNumber.deltaMCut_run=runNumber_conver;
                converNumber.deltaMCut_event=eventNumber_conver;
              }
              h_d0_e_conver->Fill(d0_e_conver);
              h_d0_p_conver->Fill(d0_p_conver);
              h_diffD0_conver->Fill(d0_e_conver-d0_p_conver);
              h_dPhi_diffD0_conver->Fill(dPhi, d0_e_conver-d0_p_conver);
              h_d0_phi_e_conver->Fill(phi_e, d0_e_conver);
              h_d0_phi_p_conver->Fill(phi_p, d0_p_conver);
              if (dD0(d0_e_conver-d0_p_conver))
              {
                if (runNumber_conver!=converNumber.diffD0Cut_run || eventNumber_conver!=converNumber.diffD0Cut_event) 
                {
                  ++converEvents.diffD0Cut;
                  converNumber.diffD0Cut_run=runNumber_conver;
                  converNumber.diffD0Cut_event=eventNumber_conver;
                }
                h_dPhi_conver->Fill(dPhi);
                h_dPhi_ee_conver->Fill(eeMass_conver, dPhi);
                if (dPhiCut(dPhi))
                {
                  if (runNumber_conver!=converNumber.dPhiCut_run || eventNumber_conver!=converNumber.dPhiCut_event) 
                  {
                    ++converEvents.dPhiCut;
                    converNumber.dPhiCut_run=runNumber_conver;
                    converNumber.dPhiCut_event=eventNumber_conver;
                    h_theta1_conver->Fill(theta1_conver);
                    h_theta2_conver->Fill(theta2_conver);
                  }
                  h_ee_conver->Fill(eeMass_conver);
                  double r1=-1/(2*curv_e_conver);
                  double r2=1/(2*curv_p_conver);
                  double R=convert(r1, r2, d0_e_conver, d0_p_conver, phi_e, phi_p);
                  h_R_conver->Fill(R);
                  if ((dsCharge * dsPlusCharge_conver)<0) oppRecon_conver+=1;
                }
              }
            }
          }
        }
      }
    }
    }
  }
  
  
  int nPhysicEvents=physicTree->GetEntries();
  int oppRecon_physic=0;
  for (int i=0; i<nPhysicEvents; ++i)
  {
    physicTree->GetEvent(i);
    
    double phi_e=atan2(py_e_physic, px_e_physic);
    double phi_p=atan2(py_p_physic, px_p_physic);
    double dPhi=deltaPhi(phi_e, phi_p);
    
    if (runNumber_physic!=physicNumber.noCut_run || eventNumber_physic!=physicNumber.noCut_event)
    {      
      if (decayMode_physic>-1) decayFrequency_physic[int(decayMode_physic)]+=1;
      ++physicEvents.noCut;
      physicNumber.noCut_run=runNumber_physic;
      physicNumber.noCut_event=eventNumber_physic;
    }
    
    if (fabs(d0_e_physic)<0.005 && fabs(d0_p_physic)<0.005 && 
        fabs(z0_e_physic)<0.05 && fabs(z0_p_physic)<0.05)
    {
    
    if (decayMode_physic==401)
    {
      if (runNumber_physic!=physicNumber.tagCut_run || eventNumber_physic!=physicNumber.tagCut_event)
      {      
        ++physicEvents.tagCut;
        physicNumber.tagCut_run=runNumber_physic;
        physicNumber.tagCut_event=eventNumber_physic;
      }
      h_dsPlusM_physic->Fill(dsPlusM_physic);
      if (dsPlusMCut(dsPlusM_physic))
      {
        if (runNumber_physic!=physicNumber.dsPlusMCut_run || eventNumber_physic!=physicNumber.dsPlusMCut_event)
        {
          ++physicEvents.dsPlusMCut;
          physicNumber.dsPlusMCut_run=runNumber_physic;
          physicNumber.dsPlusMCut_event=eventNumber_physic;
        }
        h_DeltaE_physic->Fill(DeltaE_physic); 
        // Optimization plots
        h_DeltaE_MBC_physic->Fill(DeltaE_physic, MBC_physic);
        h_DeltaE_DeltaM_physic->Fill(DeltaE_physic, DeltaM_physic);
        h_MBC_DeltaM_physic->Fill(MBC_physic, DeltaM_physic);
        //if (DeltaECut(DeltaE_physic))
        {
          if (runNumber_physic!=physicNumber.deltaECut_run || eventNumber_physic!=physicNumber.deltaECut_event)
          {
            ++physicEvents.deltaECut;
            physicNumber.deltaECut_run=runNumber_physic;
            physicNumber.deltaECut_event=eventNumber_physic;
          }
          h_MBC_physic->Fill(MBC_physic);
          if (MBCCut(MBC_physic))
          {
            if (runNumber_physic!=physicNumber.mbcCut_run || eventNumber_physic!=physicNumber.mbcCut_event)
            {
              ++physicEvents.mbcCut;
              physicNumber.mbcCut_run=runNumber_physic;
              physicNumber.mbcCut_event=eventNumber_physic;
            }
            h_DeltaM_physic->Fill(DeltaM_physic);
            if (DeltaMCut(DeltaM_physic))
            {
              if (runNumber_physic!=physicNumber.deltaMCut_run || eventNumber_physic!=physicNumber.deltaMCut_event) 
              {
                ++physicEvents.deltaMCut;
                physicNumber.deltaMCut_run=runNumber_physic;
                physicNumber.deltaMCut_event=eventNumber_physic;
              }
              h_d0_e_physic->Fill(d0_e_physic);
              h_d0_p_physic->Fill(d0_p_physic);
              h_diffD0_physic->Fill(d0_e_physic-d0_p_physic);
              h_dPhi_diffD0_physic->Fill(dPhi, d0_e_physic-d0_p_physic);
              h_d0_phi_e_physic->Fill(phi_e, d0_e_physic);
              h_d0_phi_p_physic->Fill(phi_p, d0_p_physic);
              if (dD0(d0_e_physic-d0_p_physic))
              {
                if (runNumber_physic!=physicNumber.diffD0Cut_run || eventNumber_physic!=physicNumber.diffD0Cut_event) 
                {
                  ++physicEvents.diffD0Cut;
                  physicNumber.diffD0Cut_run=runNumber_physic;
                  physicNumber.diffD0Cut_event=eventNumber_physic;
                }
                h_dPhi_physic->Fill(dPhi);
                h_dPhi_ee_physic->Fill(eeMass_physic, dPhi);
                if (dPhiCut(dPhi))
                {
                  if (runNumber_physic!=physicNumber.dPhiCut_run || eventNumber_physic!=physicNumber.dPhiCut_event) 
                  {
                    ++physicEvents.dPhiCut;
                    physicNumber.dPhiCut_run=runNumber_physic;
                    physicNumber.dPhiCut_event=eventNumber_physic;
                    std::cout<<runNumber_physic<<" "<<eventNumber_physic<<std::endl;
                  }
                  h_ee_physic->Fill(eeMass_physic);
                  double r1=-1/(2*curv_e_physic);
                  double r2=1/(2*curv_p_physic);
                  double R=convert(r1, r2, d0_e_physic, d0_p_physic, phi_e, phi_p);
                  h_R_physic->Fill(R);
                  if ((dsCharge*dsPlusCharge_physic)<0) oppRecon_physic+=1;
                  
                  // print out the e+ e- 4 vectors
                  //std::cout<<"Run "<<runNumber_physic<<", event "<<eventNumber_physic<<std::endl;
                  //std::cout<<"Positron (pX, pY, E) = ("<<px_p_physic<<","<<py_p_physic<<","<<pz_p_physic<<", "<<E_p_physic<<")"<<std::endl;
                  //std::cout<<"Electron (pX, pY, E) = ("<<px_e_physic<<","<<py_e_physic<<","<<pz_e_physic<<", "<<E_e_physic<<")"<<std::endl;                  
                }
              }
            }
          }
        }
      }
    }
    }
  }
  
  std::cout<<"Number of signal candidates = "<<nSignalEvents<<std::endl;
  for (DecayMap::iterator i_decay=decayFrequency_signal.begin(); i_decay!=decayFrequency_signal.end(); ++i_decay)
  {
    std::cout<<"Signal decay mode "<<i_decay->first<<" had "<<i_decay->second<<" decays"<<std::endl;
  }
  std::cout<<"Number of signal events after ntuplizer = "<<signalEvents.noCut<<std::endl;
  std::cout<<"Number of signal events after KKpi tag = "<<signalEvents.tagCut<<std::endl;
  std::cout<<"Number of signal events after dsPlusMCut = "<<signalEvents.dsPlusMCut<<std::endl;
  std::cout<<"Number of signal events after deltaECut = "<<signalEvents.deltaECut<<std::endl;
  std::cout<<"Number of signal events after mbcCut = "<<signalEvents.mbcCut<<std::endl;
  std::cout<<"Number of signal events after deltaMCut = "<<signalEvents.deltaMCut<<std::endl;
  std::cout<<"Number of signal events after diffD0Cut = "<<signalEvents.diffD0Cut<<std::endl;
  std::cout<<"Number of signal events after dPhiCut = "<<signalEvents.dPhiCut<<std::endl;
  std::cout<<"Number of opposite sign D_s reconstructions remaining in signal = "<<oppRecon_signal<<std::endl;
  
  std::cout<<"Number of conversion candidates = "<<nConverEvents<<std::endl;
  for (DecayMap::iterator i_decay=decayFrequency_conver.begin(); i_decay!=decayFrequency_conver.end(); ++i_decay)
  {
    std::cout<<"Conversion decay mode "<<i_decay->first<<" had "<<i_decay->second<<" decays"<<std::endl;
  }
  std::cout<<"Number of conver events after ntuplizer = "<<converEvents.noCut<<std::endl;
  std::cout<<"Number of conver events after KKpi tag = "<<converEvents.tagCut<<std::endl;
  std::cout<<"Number of conver events after dsPlusMCut = "<<converEvents.dsPlusMCut<<std::endl;
  std::cout<<"Number of conver events after deltaECut = "<<converEvents.deltaECut<<std::endl;
  std::cout<<"Number of conver events after mbcCut = "<<converEvents.mbcCut<<std::endl;
  std::cout<<"Number of conver events after deltaMCut = "<<converEvents.deltaMCut<<std::endl;
  std::cout<<"Number of conver events after diffD0Cut = "<<converEvents.diffD0Cut<<std::endl;
  std::cout<<"Number of conver events after dPhiCut = "<<converEvents.dPhiCut<<std::endl;
  std::cout<<"Number of opposite sign D_s reconstructions remaining in conver = "<<oppRecon_conver<<std::endl;
  
  std::cout<<"Number of data candidates = "<<nPhysicEvents<<std::endl;
  for (DecayMap::iterator i_decay=decayFrequency_physic.begin(); i_decay!=decayFrequency_physic.end(); ++i_decay)
  {
    std::cout<<"physicsion decay mode "<<i_decay->first<<" had "<<i_decay->second<<" decays"<<std::endl;
  }
  std::cout<<"Number of physic events after ntuplizer = "<<physicEvents.noCut<<std::endl;
  std::cout<<"Number of physic events after KKpi tag = "<<physicEvents.tagCut<<std::endl;
  std::cout<<"Number of physic events after dsPlusMCut = "<<physicEvents.dsPlusMCut<<std::endl;
  std::cout<<"Number of physic events after deltaECut = "<<physicEvents.deltaECut<<std::endl;
  std::cout<<"Number of physic events after mbcCut = "<<physicEvents.mbcCut<<std::endl;
  std::cout<<"Number of physic events after deltaMCut = "<<physicEvents.deltaMCut<<std::endl;
  std::cout<<"Number of physic events after diffD0Cut = "<<physicEvents.diffD0Cut<<std::endl;
  std::cout<<"Number of physic events after dPhiCut = "<<physicEvents.dPhiCut<<std::endl;
  std::cout<<"Number of opposite sign D_s reconstructions remaining in physic = "<<oppRecon_physic<<std::endl;
  
  TLine *line;
  float xmin, xmax, ymin, ymax;
  
  gROOT->SetStyle("Plain");
  
  TCanvas *dsPlusM = new TCanvas("dsPlusM", "", 500, 1500);  
  xmin=dsPlusMCut_center-dsPlusMCut_range;
  xmax=dsPlusMCut_center+dsPlusMCut_range;
  dsPlusM->Divide(1,3);  
  dsPlusM->cd(1);
  h_dsPlusM_signal->Draw("SAME");
  ymax=(h_dsPlusM_signal->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("dsPlusM_signal_KKpi.png"); 
  dsPlusM->cd(2);
  h_dsPlusM_conver->Draw("SAME");
  ymax=(h_dsPlusM_conver->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("dsPlusM_conver_KKpi.png");
  dsPlusM->cd(3);
  h_dsPlusM_physic->Draw("SAME");
  ymax=(h_dsPlusM_physic->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("dsPlusM_physic_KKpi.png");
  
  TCanvas *DeltaE_MBC = new TCanvas("DeltaE_MBC", "", 500, 1500);
  xmin=-deltaECut_range; xmax=deltaECut_range;
  ymin=mbcCut_center-mbcCut_range; ymax=mbcCut_center+mbcCut_range;
  DeltaE_MBC->Divide(1,3);
  DeltaE_MBC->cd(1);
  h_DeltaE_MBC_signal->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  DeltaE_MBC->cd(2);
  h_DeltaE_MBC_conver->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  DeltaE_MBC->cd(3);
  h_DeltaE_MBC_physic->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  
  TCanvas *DeltaE_DeltaM = new TCanvas("DeltaE_DeltaM", "", 500, 1500);
  xmin=-deltaECut_range; xmax=deltaECut_range;
  ymin=deltaMCut_center-deltaMCut_range; ymax=deltaMCut_center+deltaMCut_range;
  DeltaE_DeltaM->Divide(1,3);
  DeltaE_DeltaM->cd(1);
  h_DeltaE_DeltaM_signal->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  DeltaE_DeltaM->cd(2);
  h_DeltaE_DeltaM_conver->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  DeltaE_DeltaM->cd(3);
  h_DeltaE_DeltaM_physic->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  
  TCanvas *MBC_DeltaM = new TCanvas("MBC_DeltaM", "", 500, 1500);
  xmin=mbcCut_center-mbcCut_range; xmax=mbcCut_center+mbcCut_range;
  ymin=deltaMCut_center-deltaMCut_range; ymax=deltaMCut_center+deltaMCut_range;
  MBC_DeltaM->Divide(1,3);
  MBC_DeltaM->cd(1);
  h_MBC_DeltaM_signal->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  gPad->Print("mBC_DeltaM_signal_KKpi.png");
  MBC_DeltaM->cd(2);
  h_MBC_DeltaM_conver->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  gPad->Print("mBC_DeltaM_conver_KKpi.png");
  MBC_DeltaM->cd(3);
  h_MBC_DeltaM_physic->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  gPad->Print("mBC_DeltaM_physic_KKpi.png");
  
  TCanvas *DeltaE = new TCanvas("DeltaE", "", 500, 1500);
  xmin=deltaECut_center-deltaECut_range; xmax=deltaECut_center+deltaECut_range;
  DeltaE->Divide(1,3);
  DeltaE->cd(1);
  h_DeltaE_signal->Draw("SAME");
  ymax=(h_DeltaE_signal->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("DeltaE_signal_KKpi.png");
  DeltaE->cd(2);
  h_DeltaE_conver->Draw("SAME");
  ymax=(h_DeltaE_conver->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("DeltaE_conver_KKpi.png");
  DeltaE->cd(3);
  h_DeltaE_physic->Draw("SAME");
  ymax=(h_DeltaE_physic->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("DeltaE_physic_KKpi.png");
  
  TCanvas *MBC = new TCanvas("MBC", "", 500, 1500);
  xmin=mbcCut_center-mbcCut_range; xmax=mbcCut_center+mbcCut_range;
  MBC->Divide(1,3);
  MBC->cd(1);
  h_MBC_signal->Draw("SAME");
  ymax=(h_MBC_signal->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("mBC_signal_KKpi.png");
  MBC->cd(2);
  h_MBC_conver->Draw("SAME");
  ymax=(h_MBC_conver->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("mBC_conver_KKpi.png");
  MBC->cd(3);
  h_MBC_physic->Draw("SAME");
  ymax=(h_MBC_physic->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("mBC_physic_KKpi.png");
  
  TCanvas *DeltaM = new TCanvas("DeltaM", "", 500, 1500);
  xmin=deltaMCut_center-deltaMCut_range; xmax=deltaMCut_center+deltaMCut_range;
  DeltaM->Divide(1,3);
  DeltaM->cd(1);
  h_DeltaM_signal->Draw("SAME");
  ymax=(h_DeltaM_signal->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("DeltaM_signal_KKpi.png");
  DeltaM->cd(2);
  h_DeltaM_conver->Draw("SAME");
  ymax=(h_DeltaM_conver->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("DeltaM_conver_KKpi.png");
  DeltaM->cd(3);
  h_DeltaM_physic->Draw("SAME");
  ymax=(h_DeltaM_physic->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  gPad->Print("DeltaM_physic_KKpi.png");
  
  TCanvas *d0 = new TCanvas("d0");
  d0->Divide(1,4);
  d0->cd(1);
  //h_d0_e_signal->Scale(1/h_d0_e_signal->GetEntries());
  //h_d0_e_conver->Scale(1/h_d0_e_conver->GetEntries());
  h_d0_e_signal->Draw("SAME");
  d0->cd(2);
  h_d0_e_conver->Draw("SAME");
  d0->cd(3);
  //h_d0_p_signal->Scale(1/h_d0_p_signal->GetEntries());
  //h_d0_p_conver->Scale(1/h_d0_p_conver->GetEntries());
  h_d0_p_signal->Draw("SAME");
  d0->cd(4);
  h_d0_p_conver->Draw("SAME");
  
  TCanvas *d0_phi = new TCanvas("d0_phi");
  d0_phi->Divide(1,4);
  d0_phi->cd(1);
  h_d0_phi_e_signal->Draw();
  d0_phi->cd(2);
  h_d0_phi_p_signal->Draw();
  d0_phi->cd(3);
  h_d0_phi_e_conver->Draw();
  d0_phi->cd(4);
  h_d0_phi_p_conver->Draw();
    
  TCanvas *diffD0 = new TCanvas("diffD0");
  diffD0->Divide(1,3);
  diffD0->cd(1);
  h_diffD0_signal->Draw("SAME");
  ymax=(h_diffD0_signal->GetMaximum())*0.75;
  line = new TLine(diffD0Cut,ymax,diffD0Cut,0); line->Draw();
  diffD0->cd(2);
  h_diffD0_conver->Draw("SAME");
  ymax=(h_diffD0_conver->GetMaximum())*0.75;
  line = new TLine(diffD0Cut,ymax,diffD0Cut,0); line->Draw();
  diffD0->cd(3);
  h_diffD0_physic->Draw("SAME");
  ymax=(h_diffD0_physic->GetMaximum())*0.75;
  line = new TLine(diffD0Cut,ymax,diffD0Cut,0); line->Draw();
  
  TCanvas *dPhi_diffD0 = new TCanvas ("dPhi_diffD0");
  xmin=-2; xmax=dPhiCutLess;
  ymin=diffD0Cut; ymax=0.01;
  dPhi_diffD0->Divide(1,3);
  dPhi_diffD0->cd(1);
  h_dPhi_diffD0_signal->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  dPhi_diffD0->cd(2);
  h_dPhi_diffD0_conver->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  dPhi_diffD0->cd(3);
  h_dPhi_diffD0_physic->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  
  TCanvas *dPhi = new TCanvas("dPhi", "", 500, 1500);
  dPhi->Divide(1,3);
  dPhi->cd(1);
  h_dPhi_signal->Draw("SAME");
  TF1 *f_dPhi_signal=new TF1("f_dPhi_signal", "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)+[6]*exp(-0.5*((x-[7])/[8])**2)", -2., 2.);
  f_dPhi_signal->SetParLimits(2, 0.1, 0.5);
  f_dPhi_signal->SetParLimits(5, 0.5, 5);
  f_dPhi_signal->SetParLimits(8, 0.05, 0.2);
  h_dPhi_signal->Fit("f_dPhi_signal", "RL");
  ymax=(h_dPhi_signal->GetMaximum())*0.75;
  line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
  double signal_integral=f_dPhi_signal->Integral(-2., 2.);
  std::cout<<"Signal integral = "<<signal_integral<<" => "<<signal_integral*12.5<<std::endl;
  gPad->Print("dPhi_signal_KKpi.png");
  dPhi->cd(2);  
  h_dPhi_conver->Draw("SAME");
  TF1 *f_dPhi_conver=new TF1("f_dPhi_conver", "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)", -2, 2);
  f_dPhi_conver->SetParLimits(1, 0.14, 0.26); 
  f_dPhi_conver->SetParLimits(2, 0.1, 0.5);
  f_dPhi_conver->SetParLimits(4, 0.14, 0.26);
  f_dPhi_conver->SetParLimits(5, 0.05, 0.1);
  h_dPhi_conver->Fit("f_dPhi_conver", "RL");
  ymax=(h_dPhi_conver->GetMaximum())*0.75;
  line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
  double conver_integral=f_dPhi_conver->Integral(-2., 2.);
  std::cout<<"Conversion integral = "<<conver_integral<<" => "<<conver_integral*12.5<<std::endl;
  gPad->Print("dPhi_conver_KKpi.png");
  dPhi->cd(3);
  TF1 *f_dPhi_physic=new TF1("f_dPhi_physic", "[0]*([1]*f_dPhi_signal+f_dPhi_conver)", -2, 2);
  
  //f_dPhi_physic->SetParLimits(0, 1e-6, 1e-3);
  //f_dPhi_physic->SetParLimits(1, 1e-6, 1);
  f_dPhi_physic->FixParameter(2, f_dPhi_signal->GetParameter(0)/signal_integral);
  f_dPhi_physic->FixParameter(3, f_dPhi_signal->GetParameter(1));
  f_dPhi_physic->FixParameter(4, f_dPhi_signal->GetParameter(2));
  f_dPhi_physic->FixParameter(5, f_dPhi_signal->GetParameter(3)/signal_integral);
  f_dPhi_physic->FixParameter(6, f_dPhi_signal->GetParameter(4));
  f_dPhi_physic->FixParameter(7, f_dPhi_signal->GetParameter(5));
  f_dPhi_physic->FixParameter(8, f_dPhi_signal->GetParameter(6)/signal_integral);
  f_dPhi_physic->FixParameter(9, f_dPhi_signal->GetParameter(7));
  f_dPhi_physic->FixParameter(10, f_dPhi_signal->GetParameter(8));
  
  f_dPhi_physic->FixParameter(11, f_dPhi_conver->GetParameter(0)/conver_integral);
  f_dPhi_physic->FixParameter(12, f_dPhi_conver->GetParameter(1));
  f_dPhi_physic->FixParameter(13, f_dPhi_conver->GetParameter(2));
  f_dPhi_physic->FixParameter(14, f_dPhi_conver->GetParameter(3)/conver_integral);
  f_dPhi_physic->FixParameter(15, f_dPhi_conver->GetParameter(4));
  f_dPhi_physic->FixParameter(16, f_dPhi_conver->GetParameter(5));
  
  h_dPhi_physic->Draw();
  h_dPhi_physic->Fit("f_dPhi_physic", "L");
  ymax=(h_dPhi_physic->GetMaximum())*0.75;
  line = new TLine(dPhiCutLess,ymax,dPhiCutLess,0); line->Draw();
  gPad->Print("dPhi_physic_KKpi.png");
  
  TCanvas *dPhi_ee = new TCanvas("dPhi_ee");
  dPhi_ee->Divide(1,2);
  dPhi_ee->cd(1);
  h_dPhi_ee_signal->Draw();
  dPhi_ee->cd(2);
  h_dPhi_ee_conver->Draw();
  
  TCanvas *R = new TCanvas("R");
  R->Divide(1,3);
  R->cd(1);
  h_R_signal->Draw();
  R->cd(2);
  h_R_conver->Draw();
  R->cd(3);
  h_R_physic->Draw();
  
  TCanvas *ee_Mass = new TCanvas("ee_Mass");
  ee_Mass->Divide(1,2);
  ee_Mass->cd(1);
  h_ee_signal->Draw();
  ee_Mass->cd(2);
  h_ee_conver->Draw();
  
  TCanvas *electronEResolution = new TCanvas("electronEResolution");
  electronEResolution->Divide(1,2);
  electronEResolution->cd(1);
  h_electronE->Draw("box");
  electronEResolution->cd(2);
  TProfile *p_electronE=h_electronE->ProfileX();
  TF1 *f1=new TF1("f1", "[0]+[1]/(x*x)", 0.05, 0.13);
  p_electronE->Fit("f1", "R");
  
  TCanvas *electronPxResolution = new TCanvas("electronPxResolution");
  electronPxResolution->Divide(1,2);
  electronPxResolution->cd(1);
  h_electronPx->Draw("box");
  electronPxResolution->cd(2);
  TProfile *p_electronPx=h_electronPx->ProfileX();
  p_electronPx->Draw();
  
  TCanvas *electronPyResolution = new TCanvas("electronPyResolution");
  electronPyResolution->Divide(1,2);
  electronPyResolution->cd(1);
  h_electronPy->Draw("box");
  electronPyResolution->cd(2);
  TProfile *p_electronPy=h_electronPy->ProfileX();
  p_electronPy->Draw();
  
  TCanvas *electronPzResolution = new TCanvas("electronPzResolution");
  electronPzResolution->Divide(1,2);
  electronPzResolution->cd(1);
  h_electronPz->Draw("box");
  electronPzResolution->cd(2);
  TProfile *p_electronPz=h_electronPz->ProfileX();
  p_electronPz->Draw();
  
  TCanvas *theta1 = new TCanvas("theta1");
  theta1->Divide(1,2);
  theta1->cd(1);
  h_theta1_signal->Draw();
  theta1->cd(2);
  h_theta1_conver->Draw();
  
  TCanvas *theta2 = new TCanvas("theta2");
  theta2->Divide(1,2);
  theta2->cd(1);
  h_theta2_signal->Draw();
  theta2->cd(2);
  h_theta2_conver->Draw();
  
  return 0;
}
