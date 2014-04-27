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
double dsPlusMCut_center=1.96849, dsPlusMCut_range=0.011; //0.011
double deltaECut_range=0.05;
double mbcCut_center=2.112, mbcCut_range=0.02;
double deltaMCut_center=0.1455, deltaMCut_range=0.0085;
double diffD0Cut=-0.004;
double dPhiCutLess=0.1;

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
  return (fabs(DeltaE)<deltaECut_range);
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

double invariantMassSq(double px1, double py1, double pz1, double E1, double px2, double py2, double pz2, double E2)
{
  double inv=pow(E1+E2,2)-pow(px1+px2,2)-pow(py1+py2,2)-pow(pz1+pz2,2);
  return pow(inv, 0.5);
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

int DsTaggedAnalysis_all_KKpi()
{

  TFile *signalFile = new TFile("DssteeBIGFILES/MyDChainFile_MC_vtosll_Dsp_KKpi.root");
  signalFile->cd("DsTaggedDecaysProc");
  TTree *signalTree1 = (TTree*)gDirectory->Get("nt1"); // pion
  TTree *signalTree2 = (TTree*)gDirectory->Get("nt2"); // paramterized pion
  //TFile *signalFile1 = new TFile("DssteeBIGFILES/MyDChainFile_MC_vtosll_electronFit_KKpi.root");
  //signalFile1->cd("DsTaggedDecaysProc");
  TTree *signalTree3 = (TTree*)gDirectory->Get("nt3"); // electron
  
  TFile *converFile = new TFile("DssteeBIGFILES/MyDChainFile_MCgamma_KKpi.root");
  converFile->cd("DsTaggedDecaysProc");
  TTree *converTree = (TTree*)gDirectory->Get("nt5");
  
  TFile *physicFile = new TFile("DssteeBIGFILES/MyDChainFile_Data_dtag_10.root");
  //TFile *physicFile = new TFile("DssteeBIGFILES/MyDChainFile_GenericMC_60.root");
  physicFile->cd("DsTaggedDecaysProc");
  TTree *physicTree = (TTree*)gDirectory->Get("nt5");
  
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
  
  TH1D *h_dsPlusM_signal = new TH1D("h_dsPlusM_signal", "m_{D_{S}^{+}} Signal Sample; GeV", 100, 1.9, 2.1); h_dsPlusM_signal->SetLineColor(kRed);
  TH1D *h_dsPlusM_conver = new TH1D("h_dsPlusM_conver", "m_{D_{S}^{+}} Conversion Background Sample; GeV", 100, 1.9, 2.1); h_dsPlusM_conver->SetLineColor(kBlue);
  TH1D *h_dsPlusM_physic = new TH1D("h_dsPlusM_physic", "m_{D_{S}^{+}} Data; GeV", 100, 1.9, 2.1); h_dsPlusM_physic->SetLineColor(kGreen);
  TH1D *h_DeltaE_signal1 = new TH1D("h_DeltaE_signal1", "#DeltaE Signal Sample; GeV", 100, -0.1, 0.2); h_DeltaE_signal1->SetLineColor(kRed);
  TH1D *h_DeltaE_signal2 = new TH1D("h_DeltaE_signal2", "#DeltaE Signal Sample; GeV", 100, -0.1, 0.2); h_DeltaE_signal2->SetLineColor(kGreen);
  TH1D *h_DeltaE_signal3 = new TH1D("h_DeltaE_signal3", "#DeltaE Signal Sample; GeV", 100, -0.1, 0.2); h_DeltaE_signal3->SetLineColor(kBlue);
  TH1D *h_DeltaE_conver = new TH1D("h_DeltaE_conver", "#DeltaE Conversion Background Sample; GeV", 100, -0.1, 0.2); h_DeltaE_conver->SetLineColor(kBlue);
  TH1D *h_DeltaE_physic = new TH1D("h_DeltaE_physic", "#DeltaE Data; GeV", 100, -0.1, 0.2); h_DeltaE_physic->SetLineColor(kGreen);
  TH1D *h_MBC_signal1 = new TH1D("h_MBC_signal1", "m_{BC} Signal Sample; GeV", 100, 2., 2.2); h_MBC_signal1->SetLineColor(kRed);
  TH1D *h_MBC_signal3 = new TH1D("h_MBC_signal3", "m_{BC} Signal Sample; GeV", 100, 2., 2.2); h_MBC_signal3->SetLineColor(kBlue);
  TH1D *h_MBC_conver = new TH1D("h_MBC_conver", "m_{BC} Conversion Background; GeV", 100, 2., 2.2); h_MBC_conver->SetLineColor(kBlue);
  TH1D *h_MBC_physic = new TH1D("h_MBC_physic", "m_{BC} Data; GeV", 100, 2., 2.2); h_MBC_physic->SetLineColor(kGreen);
  TH1D *h_DeltaM_signal1 = new TH1D("h_DeltaM_signal1", "#deltaM Signal Sample", 100, 0.0, 0.2); h_DeltaM_signal1->SetLineColor(kRed);
  TH1D *h_DeltaM_signal2 = new TH1D("h_DeltaM_signal2", "#deltaM Signal Sample", 100, 0.0, 0.2); h_DeltaM_signal2->SetLineColor(kGreen);
  TH1D *h_DeltaM_signal3 = new TH1D("h_DeltaM_signal3", "#deltaM Signal Sample", 100, 0.0, 0.2); h_DeltaM_signal3->SetLineColor(kBlue);
  TH1D *h_DeltaM_conver = new TH1D("h_DeltaM_conver", "#deltaM Conversion Background Sample; GeV", 100, 0.0, 0.2); h_DeltaM_conver->SetLineColor(kBlue);
  TH1D *h_DeltaM_physic = new TH1D("h_DeltaM_physic", "#deltaM Data; GeV", 100, 0.0, 0.2); h_DeltaM_physic->SetLineColor(kGreen);
  
  TH2D *h_electronE_signal1 = new TH2D("h_electronE_signal1", "Electron Energy Resolution; Reco Electron Energy (GeV)", 100, 0.0, 0.16, 100, -0.03, 0.03);
  TH2D *h_electronE_signal2 = new TH2D("h_electronE_signal2", "Electron Energy Resolution; Reco Electron Energy (GeV)", 100, 0.0, 0.16, 100, -0.03, 0.03);
  TH2D *h_electronE_signal3 = new TH2D("h_electronE_signal3", "Electron Energy Resolution; Reco Electron Energy (GeV)", 100, 0.0, 0.16, 100, -0.03, 0.03);

  // Pion Fit
  signalTree1->SetBranchAddress("dsPlusM", &(dsPlusM_signal));
  signalTree1->SetBranchAddress("DecayMode", &(decayMode_signal));
  signalTree1->SetBranchAddress("DeltaE", &(DeltaE_signal));
  signalTree1->SetBranchAddress("MBC", &(MBC_signal));
  signalTree1->SetBranchAddress("DeltaM", &(DeltaM_signal));
  signalTree1->SetBranchAddress("kElectron1E_MC", &(E_e_signal_MC));  
  signalTree1->SetBranchAddress("kElectron1Px_MC", &(px_e_signal_MC));
  signalTree1->SetBranchAddress("kElectron1Py_MC", &(py_e_signal_MC));
  signalTree1->SetBranchAddress("kElectron1Pz_MC", &(pz_e_signal_MC));
  signalTree1->SetBranchAddress("kElectron2E_MC", &(E_p_signal_MC));
  signalTree1->SetBranchAddress("kElectron2Px_MC", &(px_p_signal_MC));
  signalTree1->SetBranchAddress("kElectron2Py_MC", &(py_p_signal_MC));
  signalTree1->SetBranchAddress("kElectron2Pz_MC", &(pz_p_signal_MC));
  signalTree1->SetBranchAddress("kElectron1D0_reco", &(d0_e_signal));
  signalTree1->SetBranchAddress("kElectron2D0_reco", &(d0_p_signal));
  signalTree1->SetBranchAddress("kElectron1Z0_reco", &(z0_e_signal));
  signalTree1->SetBranchAddress("kElectron2Z0_reco", &(z0_p_signal));
  signalTree1->SetBranchAddress("kElectron1Px_reco", &(px_e_signal));
  signalTree1->SetBranchAddress("kElectron1Py_reco", &(py_e_signal));
  signalTree1->SetBranchAddress("kElectron1Pz_reco", &(pz_e_signal));
  signalTree1->SetBranchAddress("kElectron2Px_reco", &(px_p_signal));
  signalTree1->SetBranchAddress("kElectron2Py_reco", &(py_p_signal));
  signalTree1->SetBranchAddress("kElectron2Pz_reco", &(pz_p_signal));
  signalTree1->SetBranchAddress("kElectron1E_reco", &(E_e_signal));
  signalTree1->SetBranchAddress("kElectron2E_reco", &(E_p_signal));
  signalTree1->SetBranchAddress("keeMass_reco", &(eeMass_signal));
  
  for (int i=0; i<signalTree1->GetEntries(); ++i)
  {
    signalTree1->GetEvent(i);   
    
    //if (fabs(d0_e_signal)<0.005 && fabs(d0_p_signal)<0.005 &&
    //    fabs(z0_e_signal)<0.05 && fabs(z0_p_signal)<0.05)
    {
      if (decayMode_signal==401)
      { 
        h_dsPlusM_signal->Fill(dsPlusM_signal);
        if (dsPlusMCut(dsPlusM_signal) && DeltaECut(DeltaE_signal) && MBCCut(MBC_signal))
        {
          h_DeltaE_signal1->Fill(DeltaE_signal);
          h_MBC_signal1->Fill(MBC_signal); 
          h_DeltaM_signal1->Fill(DeltaM_signal);
          h_electronE_signal1->Fill(E_e_signal, E_e_signal-E_e_signal_MC);
          h_electronE_signal1->Fill(E_p_signal, E_p_signal-E_p_signal_MC);
        }
      }
    }
  }
  std::cout<<"Number of pion fitted signal events = "<<h_DeltaE_signal1->GetEntries()<<std::endl;
  
  
  // Parameterized Pion Fit
  signalTree2->SetBranchAddress("dsPlusM", &(dsPlusM_signal));
  signalTree2->SetBranchAddress("DecayMode", &(decayMode_signal));
  signalTree2->SetBranchAddress("DeltaE", &(DeltaE_signal));
  signalTree2->SetBranchAddress("MBC", &(MBC_signal));
  signalTree2->SetBranchAddress("DeltaM", &(DeltaM_signal));
  signalTree2->SetBranchAddress("kElectron1E_MC", &(E_e_signal_MC));  
  signalTree2->SetBranchAddress("kElectron1Px_MC", &(px_e_signal_MC));
  signalTree2->SetBranchAddress("kElectron1Py_MC", &(py_e_signal_MC));
  signalTree2->SetBranchAddress("kElectron1Pz_MC", &(pz_e_signal_MC));
  signalTree2->SetBranchAddress("kElectron2E_MC", &(E_p_signal_MC));
  signalTree2->SetBranchAddress("kElectron2Px_MC", &(px_p_signal_MC));
  signalTree2->SetBranchAddress("kElectron2Py_MC", &(py_p_signal_MC));
  signalTree2->SetBranchAddress("kElectron2Pz_MC", &(pz_p_signal_MC));
  signalTree2->SetBranchAddress("kElectron1D0_reco", &(d0_e_signal));
  signalTree2->SetBranchAddress("kElectron2D0_reco", &(d0_p_signal));
  signalTree2->SetBranchAddress("kElectron1Z0_reco", &(z0_e_signal));
  signalTree2->SetBranchAddress("kElectron2Z0_reco", &(z0_p_signal));
  signalTree2->SetBranchAddress("kElectron1Px_reco", &(px_e_signal));
  signalTree2->SetBranchAddress("kElectron1Py_reco", &(py_e_signal));
  signalTree2->SetBranchAddress("kElectron1Pz_reco", &(pz_e_signal));
  signalTree2->SetBranchAddress("kElectron2Px_reco", &(px_p_signal));
  signalTree2->SetBranchAddress("kElectron2Py_reco", &(py_p_signal));
  signalTree2->SetBranchAddress("kElectron2Pz_reco", &(pz_p_signal));
  signalTree2->SetBranchAddress("kElectron1E_reco", &(E_e_signal));
  signalTree2->SetBranchAddress("kElectron2E_reco", &(E_p_signal));
  signalTree2->SetBranchAddress("keeMass_reco", &(eeMass_signal));
  
  for (int i=0; i<signalTree2->GetEntries(); ++i)
  {
    signalTree2->GetEvent(i);   
    
    //if (fabs(d0_e_signal)<0.005 && fabs(d0_p_signal)<0.005 &&
    //    fabs(z0_e_signal)<0.05 && fabs(z0_p_signal)<0.05)
    {
      if (decayMode_signal==401)
      { 
        h_dsPlusM_signal->Fill(dsPlusM_signal);
        if (dsPlusMCut(dsPlusM_signal) && DeltaECut(DeltaE_signal) && MBCCut(MBC_signal))
        {
          h_DeltaE_signal2->Fill(DeltaE_signal);
          h_DeltaM_signal2->Fill(DeltaM_signal);
          h_electronE_signal2->Fill(E_e_signal, E_e_signal-E_e_signal_MC);
          h_electronE_signal2->Fill(E_p_signal, E_p_signal-E_p_signal_MC);
        }
      }
    }
  }
  std::cout<<"Number of parameterized-pion fitted signal events = "<<h_DeltaE_signal2->GetEntries()<<std::endl;
  
  // Electron Fit
  signalTree3->SetBranchAddress("dsPlusM", &(dsPlusM_signal));
  signalTree3->SetBranchAddress("DecayMode", &(decayMode_signal));
  signalTree3->SetBranchAddress("DeltaE", &(DeltaE_signal));
  signalTree3->SetBranchAddress("MBC", &(MBC_signal));
  signalTree3->SetBranchAddress("DeltaM", &(DeltaM_signal));
  signalTree3->SetBranchAddress("kElectron1E_MC", &(E_e_signal_MC));  
  signalTree3->SetBranchAddress("kElectron1Px_MC", &(px_e_signal_MC));
  signalTree3->SetBranchAddress("kElectron1Py_MC", &(py_e_signal_MC));
  signalTree3->SetBranchAddress("kElectron1Pz_MC", &(pz_e_signal_MC));
  signalTree3->SetBranchAddress("kElectron2E_MC", &(E_p_signal_MC));
  signalTree3->SetBranchAddress("kElectron2Px_MC", &(px_p_signal_MC));
  signalTree3->SetBranchAddress("kElectron2Py_MC", &(py_p_signal_MC));
  signalTree3->SetBranchAddress("kElectron2Pz_MC", &(pz_p_signal_MC));
  signalTree3->SetBranchAddress("kElectron1D0_reco", &(d0_e_signal));
  signalTree3->SetBranchAddress("kElectron2D0_reco", &(d0_p_signal));
  signalTree3->SetBranchAddress("kElectron1Z0_reco", &(z0_e_signal));
  signalTree3->SetBranchAddress("kElectron2Z0_reco", &(z0_p_signal));
  signalTree3->SetBranchAddress("kElectron1Px_reco", &(px_e_signal));
  signalTree3->SetBranchAddress("kElectron1Py_reco", &(py_e_signal));
  signalTree3->SetBranchAddress("kElectron1Pz_reco", &(pz_e_signal));
  signalTree3->SetBranchAddress("kElectron2Px_reco", &(px_p_signal));
  signalTree3->SetBranchAddress("kElectron2Py_reco", &(py_p_signal));
  signalTree3->SetBranchAddress("kElectron2Pz_reco", &(pz_p_signal));
  signalTree3->SetBranchAddress("kElectron1E_reco", &(E_e_signal));
  signalTree3->SetBranchAddress("kElectron2E_reco", &(E_p_signal));
  signalTree3->SetBranchAddress("keeMass_reco", &(eeMass_signal));
  
  for (int i=0; i<signalTree3->GetEntries(); ++i)
  {
    signalTree3->GetEvent(i);
    
    if (fabs(d0_e_signal)<0.005 && fabs(d0_p_signal)<0.005 &&
        fabs(z0_e_signal)<0.05 && fabs(z0_p_signal)<0.05)
    {
      if (decayMode_signal==401)
      { 
        h_dsPlusM_signal->Fill(dsPlusM_signal);
        if (dsPlusMCut(dsPlusM_signal) && DeltaECut(DeltaE_signal) && MBCCut(MBC_signal))
        {
          h_DeltaE_signal3->Fill(DeltaE_signal);
          h_MBC_signal3->Fill(MBC_signal);
          h_DeltaM_signal3->Fill(DeltaM_signal);
          h_electronE_signal3->Fill(E_e_signal, E_e_signal-E_e_signal_MC);
          h_electronE_signal3->Fill(E_p_signal, E_p_signal-E_p_signal_MC);
        }
      }
    }
  }
  std::cout<<"Number of electron fitted signal events = "<<h_DeltaE_signal3->GetEntries()<<std::endl;
  
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
  
  int nConverEvents=converTree->GetEntries();
  for (int i=0; i<nConverEvents; ++i)
  {
    converTree->GetEvent(i);
    
    if (fabs(d0_e_conver)<0.005 && fabs(d0_p_conver)<0.005 &&
        fabs(z0_e_conver)<0.05 && fabs(z0_p_conver)<0.05)
    {    
      if (decayMode_conver==401)
      { 
        h_dsPlusM_conver->Fill(dsPlusM_conver); 
        h_DeltaE_conver->Fill(DeltaE_conver);
        h_MBC_conver->Fill(MBC_conver); 
        h_DeltaM_conver->Fill(DeltaM_conver);
      }
    }
  }
  
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
  physicTree->SetBranchAddress("kElectron2Px_reco", &(px_p_physic));
  physicTree->SetBranchAddress("kElectron2Py_reco", &(py_p_physic));
  physicTree->SetBranchAddress("kElectron1E_reco", &(E_e_physic));
  physicTree->SetBranchAddress("kElectron2E_reco", &(E_p_physic));
  physicTree->SetBranchAddress("kElectron1Curv_reco", &(curv_e_physic));
  physicTree->SetBranchAddress("kElectron2Curv_reco", &(curv_p_physic));
  physicTree->SetBranchAddress("keeMass_reco", &(eeMass_physic));
  
  int nPhysicEvents=physicTree->GetEntries();
  int oppRecon_physic=0;
  for (int i=0; i<nPhysicEvents; ++i)
  {
    physicTree->GetEvent(i);
    
    if (fabs(d0_e_physic)<0.005 && fabs(d0_p_physic)<0.005 && 
        fabs(z0_e_physic)<0.05 && fabs(z0_p_physic)<0.05)
    {    
      if (decayMode_physic==401)
      {
        h_dsPlusM_physic->Fill(dsPlusM_physic); 
        h_DeltaE_physic->Fill(DeltaE_physic); 
        h_MBC_physic->Fill(MBC_physic); 
        h_DeltaM_physic->Fill(DeltaM_physic);
      }
    }
  }
  
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
  
  TCanvas *DeltaE = new TCanvas("DeltaE", "", 500, 1500);
  xmin=-deltaECut_range; xmax=deltaECut_range;
  DeltaE->Divide(1,3);
  DeltaE->cd(1);
  h_DeltaE_signal3->Draw("SAME");
  h_DeltaE_signal1->Draw("SAME");
  h_DeltaE_signal2->Draw("SAME");  
  ymax=(h_DeltaE_signal1->GetMaximum())*0.75;
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
  h_MBC_signal3->Draw("SAME");
  h_MBC_signal1->Draw("SAME");
  ymax=(h_MBC_signal1->GetMaximum())*0.75;
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
  h_DeltaM_signal3->Draw("SAME");
  h_DeltaM_signal1->Draw("SAME");
  h_DeltaM_signal2->Draw("SAME");  
  ymax=(h_DeltaM_signal1->GetMaximum())*0.75;
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
  
  TCanvas *electronEResolution = new TCanvas("electronEResolution", "", 500, 1500);
  electronEResolution->Divide(1,3);
  electronEResolution->cd(1);
  h_electronE_signal1->Draw("box");
  electronEResolution->cd(2);
  h_electronE_signal2->Draw("box");
  gPad->Print("electronEResolution_param.png");
  electronEResolution->cd(3);
  h_electronE_signal3->Draw("box");
  gPad->Print("electronEResolution_eFit.png");
  
  return 0;
}
