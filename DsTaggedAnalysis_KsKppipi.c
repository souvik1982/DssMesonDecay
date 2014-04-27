#include <iostream>
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TFile.h"
#include "TPaveStats.h"
#include <set>
#include <map>

double pi=3.14159265358979;
double dsPlusMCut_center=1.96849, dsPlusMCut_range=0.02;
double deltaECut_range=0.019;
double mbcCut_center=2.112, mbcCut_range=0.008;
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

int DsTaggedAnalysis_KsKppipi()
{

  TFile *signalFile = new TFile("MyDChainFile_MC_vtosll_KsKppipi.root");
  //TFile *signalFile = new TFile("MyDChainFile_Data_dtag_1.root");
  signalFile->cd("DsTaggedDecaysProc");
  TTree *signalTree = (TTree*)gDirectory->Get("nt5");
  
  TFile *converFile = new TFile("MyDChainFile_MCgamma_KsKppipi.root");
  converFile->cd("DsTaggedDecaysProc");
  TTree *converTree = (TTree*)gDirectory->Get("nt5");
  
  TFile *physicFile = new TFile("MyDChainFile_Data_dtag_10.root");
  physicFile->cd("DsTaggedDecaysProc");
  TTree *physicTree = (TTree*)gDirectory->Get("nt5");
  
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
  float d0_e_signal, d0_p_signal, px_e_signal, py_e_signal, px_p_signal, py_p_signal, E_e_signal, E_p_signal, curv_e_signal, curv_p_signal;
  float d0_e_conver, d0_p_conver, px_e_conver, py_e_conver, px_p_conver, py_p_conver, E_e_conver, E_p_conver, curv_e_conver, curv_p_conver;
  float d0_e_physic, d0_p_physic, px_e_physic, py_e_physic, px_p_physic, py_p_physic, E_e_physic, E_p_physic, curv_e_physic, curv_p_physic;
  
  typedef std::map<int, float> DecayMap;
  DecayMap decayFrequency_signal, decayFrequency_conver, decayFrequency_physic;
  
  TH1D *h_phiM_signal = new TH1D("h_phiM_signal", "h_phiM_signal", 100, .95, 1.15); h_phiM_signal->SetLineColor(kRed);
  TH1D *h_phiM_conver = new TH1D("h_phiM_conver", "h_phiM_conver", 100, .95, 1.15); h_phiM_conver->SetLineColor(kBlue);
  TH1D *h_phiM_physic = new TH1D("h_phiM_physic", "h_phiM_physic", 100, .95, 1.15); h_phiM_physic->SetLineColor(kGreen);
  TH1D *h_dsPlusM_signal = new TH1D("h_dsPlusM_signal", "h_dsPlusM_signal", 100, 1.9, 2.1); h_dsPlusM_signal->SetLineColor(kRed);
  TH1D *h_dsPlusM_conver = new TH1D("h_dsPlusM_conver", "h_dsPlusM_conver", 100, 1.9, 2.1); h_dsPlusM_conver->SetLineColor(kBlue);
  TH1D *h_dsPlusM_physic = new TH1D("h_dsPlusM_physic", "h_dsPlusM_physic", 100, 1.9, 2.1); h_dsPlusM_physic->SetLineColor(kGreen);
  TH1D *h_DeltaE_signal = new TH1D("h_DeltaE_signal", "DeltaE", 100, -0.1, 0.2); h_DeltaE_signal->SetLineColor(kRed);
  TH1D *h_DeltaE_conver = new TH1D("h_DeltaE_conver", "DeltaE", 100, -0.1, 0.2); h_DeltaE_conver->SetLineColor(kBlue);
  TH1D *h_DeltaE_physic = new TH1D("h_DeltaE_physic", "DeltaE", 100, -0.1, 0.2); h_DeltaE_physic->SetLineColor(kGreen);
  TH1D *h_MBC_signal = new TH1D("h_MBC_signal", "MBC", 100, 2., 2.2); h_MBC_signal->SetLineColor(kRed);
  TH1D *h_MBC_conver = new TH1D("h_MBC_conver", "MBC", 100, 2., 2.2); h_MBC_conver->SetLineColor(kBlue);
  TH1D *h_MBC_physic = new TH1D("h_MBC_physic", "MBC", 100, 2., 2.2); h_MBC_physic->SetLineColor(kGreen);
  TH1D *h_DeltaM_signal = new TH1D("h_DeltaM_signal", "DeltaM", 100, 0.0, 0.2); h_DeltaM_signal->SetLineColor(kRed);
  TH1D *h_DeltaM_conver = new TH1D("h_DeltaM_conver", "DeltaM", 100, 0.0, 0.2); h_DeltaM_conver->SetLineColor(kBlue);
  TH1D *h_DeltaM_physic = new TH1D("h_DeltaM_physic", "DeltaM", 100, 0.0, 0.2); h_DeltaM_physic->SetLineColor(kGreen);
  
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
  TH1D *h_diffD0_signal = new TH1D("h_diffD0_signal", "h_diffD0_signal", 50, -0.01, 0.01); h_diffD0_signal->SetLineColor(kRed);
  TH1D *h_diffD0_conver = new TH1D("h_diffD0_conver", "h_diffD0_conver", 50, -0.01, 0.01); h_diffD0_conver->SetLineColor(kBlue);
  TH1D *h_diffD0_physic = new TH1D("h_diffD0_physic", "h_diffD0_physic", 50, -0.01, 0.01); h_diffD0_physic->SetLineColor(kGreen);
  TH1D *h_dPhi_signal = new TH1D("h_dPhi_signal", "dPhi", 50, -2., 2.); h_dPhi_signal->SetLineColor(kRed);
  TH1D *h_dPhi_conver = new TH1D("h_dPhi_conver", "dPhi", 50, -2., 2.); h_dPhi_conver->SetLineColor(kBlue);
  TH1D *h_dPhi_physic = new TH1D("h_dPhi_physic", "dPhi", 50, -2., 2.); h_dPhi_physic->SetLineColor(kGreen);
  TH2D *h_dPhi_ee_signal = new TH2D("h_dPhi_ee_signal", "dPhi_ee_signal", 50, 0., 0.2, 50, -2., 2.);
  TH2D *h_dPhi_ee_conver = new TH2D("h_dPhi_ee_conver", "dPhi_ee_conver", 50, 0., 0.2, 50, -2., 2.);
  TH2D *h_dPhi_ee_physic = new TH2D("h_dPhi_ee_physic", "dPhi_ee_physic", 50, 0., 0.2, 50, -2., 2.);
  TH2D *h_d0_phi_e_signal = new TH2D("h_d0_phi_e_signal", "h_d0_phi_e_signal", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_p_signal = new TH2D("h_d0_phi_p_signal", "h_d0_phi_p_signal", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_e_conver = new TH2D("h_d0_phi_e_conver", "h_d0_phi_e_conver", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_p_conver = new TH2D("h_d0_phi_p_conver", "h_d0_phi_p_conver", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_e_physic = new TH2D("h_d0_phi_e_physic", "h_d0_phi_e_physic", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_p_physic = new TH2D("h_d0_phi_p_physic", "h_d0_phi_p_physic", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_dPhi_sumPhi_signal = new TH2D("h_dPhi_sumPhi_signal", "dPhi_sumPhi_signal", 50, -2*pi, 2*pi, 50, -2., 2.);
  TH2D *h_dPhi_sumPhi_conver = new TH2D("h_dPhi_sumPhi_conver", "dPhi_sumPhi_conver", 50, -2*pi, 2*pi, 50, -2., 2.);
  TH2D *h_dPhi_sumPhi_physic = new TH2D("h_dPhi_sumPhi_physic", "dPhi_sumPhi_physic", 50, -2*pi, 2*pi, 50, -2., 2.);
  TH2D *h_dPhi_diffD0_signal = new TH2D("h_dPhi_diffD0_signal", "h_dPhi_diffD0_signal", 50, -2., 2., 50, -0.01, 0.01);
  TH2D *h_dPhi_diffD0_conver = new TH2D("h_dPhi_diffD0_conver", "h_dPhi_diffD0_conver", 50, -2., 2., 50, -0.01, 0.01);
  TH2D *h_dPhi_diffD0_physic = new TH2D("h_dPhi_diffD0_physic", "h_dPhi_diffD0_physic", 50, -2., 2., 50, -0.01, 0.01);
  TH1D *h_R_signal = new TH1D("h_R_signal", "h_R_signal", 50, -.4, .4); h_R_signal->SetLineColor(kRed);
  TH1D *h_R_conver = new TH1D("h_R_conver", "h_R_conver", 50, -.4, .4); h_R_conver->SetLineColor(kBlue);
  TH1D *h_R_physic = new TH1D("h_R_physic", "h_R_physic", 50, -.4, .4); h_R_physic->SetLineColor(kGreen);
  TH1D *h_ee_signal = new TH1D("h_ee_signal", "h_ee_signal", 20, 0.0, 0.2); h_ee_signal->SetLineColor(kRed);
  TH1D *h_ee_conver = new TH1D("h_ee_conver", "h_ee_conver", 20, 0.0, 0.2); h_ee_conver->SetLineColor(kBlue);
  TH1D *h_ee_physic = new TH1D("h_ee_physic", "h_ee_physic", 20, 0.0, 0.2); h_ee_physic->SetLineColor(kGreen);
  
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
  signalTree->SetBranchAddress("kElectron1Px_reco", &(px_e_signal));
  signalTree->SetBranchAddress("kElectron1Py_reco", &(py_e_signal));
  signalTree->SetBranchAddress("kElectron2Px_reco", &(px_p_signal));
  signalTree->SetBranchAddress("kElectron2Py_reco", &(py_p_signal));
  signalTree->SetBranchAddress("kElectron1E_reco", &(E_e_signal));
  signalTree->SetBranchAddress("kElectron2E_reco", &(E_p_signal));
  signalTree->SetBranchAddress("kElectron1Curv_reco", &(curv_e_signal));
  signalTree->SetBranchAddress("kElectron2Curv_reco", &(curv_p_signal));
  signalTree->SetBranchAddress("keeMass_reco", &(eeMass_signal));
  
  converTree->SetBranchAddress("Run", &(runNumber_conver));
  converTree->SetBranchAddress("Event", &(eventNumber_conver));
  converTree->SetBranchAddress("dsPlusM", &(dsPlusM_conver));
  converTree->SetBranchAddress("dsPlusCharge", &(dsPlusCharge_conver));
  converTree->SetBranchAddress("DecayMode", &(decayMode_conver));
  converTree->SetBranchAddress("DeltaE", &(DeltaE_conver));
  converTree->SetBranchAddress("MBC", &(MBC_conver));
  converTree->SetBranchAddress("DeltaM", &(DeltaM_conver));
  converTree->SetBranchAddress("kElectron1D0_reco", &(d0_e_conver));
  converTree->SetBranchAddress("kElectron2D0_reco", &(d0_p_conver));
  converTree->SetBranchAddress("kElectron1Px_reco", &(px_e_conver));
  converTree->SetBranchAddress("kElectron1Py_reco", &(py_e_conver));
  converTree->SetBranchAddress("kElectron2Px_reco", &(px_p_conver));
  converTree->SetBranchAddress("kElectron2Py_reco", &(py_p_conver));
  converTree->SetBranchAddress("kElectron1E_reco", &(E_e_conver));
  converTree->SetBranchAddress("kElectron2E_reco", &(E_p_conver));
  converTree->SetBranchAddress("kElectron1Curv_reco", &(curv_e_conver));
  converTree->SetBranchAddress("kElectron2Curv_reco", &(curv_p_conver));
  converTree->SetBranchAddress("keeMass_reco", &(eeMass_conver));
  
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
  physicTree->SetBranchAddress("kElectron1Px_reco", &(px_e_physic));
  physicTree->SetBranchAddress("kElectron1Py_reco", &(py_e_physic));
  physicTree->SetBranchAddress("kElectron2Px_reco", &(px_p_physic));
  physicTree->SetBranchAddress("kElectron2Py_reco", &(py_p_physic));
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
    double sumPhi=phi_e+phi_p;
    
    if (runNumber_signal!=signalNumber.noCut_run || eventNumber_signal!=signalNumber.noCut_event)
    {      
      if (decayMode_signal>-1) decayFrequency_signal[int(decayMode_signal)]+=1;
      ++signalEvents.noCut;
      signalNumber.noCut_run=runNumber_signal;
      signalNumber.noCut_event=eventNumber_signal;
    }
    
    if (decayMode_signal==405)
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
        if (DeltaECut(DeltaE_signal))
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
                h_dPhi_sumPhi_signal->Fill(sumPhi, dPhi);
                h_dPhi_ee_signal->Fill(eeMass_signal, dPhi);
                if (dPhiCut(dPhi))
                {
                  if (runNumber_signal!=signalNumber.dPhiCut_run || eventNumber_signal!=signalNumber.dPhiCut_event) 
                  {
                    ++signalEvents.dPhiCut;
                    signalNumber.dPhiCut_run=runNumber_signal;
                    signalNumber.dPhiCut_event=eventNumber_signal;
                  }
                  h_ee_signal->Fill(eeMass_signal);
                  double r1=-1/(2*curv_e_signal);
                  double r2=1/(2*curv_p_signal);
                  double R=convert(r1, r2, d0_e_signal, d0_p_signal, phi_e, phi_p);
                  h_R_signal->Fill(R);
                  if (dsPlusCharge_signal<0) oppRecon_signal+=1;
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
    double sumPhi=phi_e+phi_p;
    
    if (runNumber_conver!=converNumber.noCut_run || eventNumber_conver!=converNumber.noCut_event)
    {      
      if (decayMode_conver>-1) decayFrequency_conver[int(decayMode_conver)]+=1;
      ++converEvents.noCut;
      converNumber.noCut_run=runNumber_conver;
      converNumber.noCut_event=eventNumber_conver;
    }
    
    if (decayMode_conver==405)
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
        if (DeltaECut(DeltaE_conver))
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
                h_dPhi_sumPhi_conver->Fill(sumPhi, dPhi);
                h_dPhi_ee_conver->Fill(eeMass_conver, dPhi);
                if (dPhiCut(dPhi))
                {
                  if (runNumber_conver!=converNumber.dPhiCut_run || eventNumber_conver!=converNumber.dPhiCut_event) 
                  {
                    ++converEvents.dPhiCut;
                    converNumber.dPhiCut_run=runNumber_conver;
                    converNumber.dPhiCut_event=eventNumber_conver;
                  }
                  h_ee_conver->Fill(eeMass_conver);
                  double r1=-1/(2*curv_e_conver);
                  double r2=1/(2*curv_p_conver);
                  double R=convert(r1, r2, d0_e_conver, d0_p_conver, phi_e, phi_p);
                  h_R_conver->Fill(R);
                  if (dsPlusCharge_conver<0) oppRecon_conver+=1;
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
    double sumPhi=phi_e+phi_p;
    
    if (runNumber_physic!=physicNumber.noCut_run || eventNumber_physic!=physicNumber.noCut_event)
    {      
      if (decayMode_physic>-1) decayFrequency_physic[int(decayMode_physic)]+=1;
      ++physicEvents.noCut;
      physicNumber.noCut_run=runNumber_physic;
      physicNumber.noCut_event=eventNumber_physic;
    }
    
    if (decayMode_physic==405)
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
        if (DeltaECut(DeltaE_physic))
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
                h_dPhi_sumPhi_physic->Fill(sumPhi, dPhi);
                h_dPhi_ee_physic->Fill(eeMass_physic, dPhi);
                if (dPhiCut(dPhi))
                {
                  if (runNumber_physic!=physicNumber.dPhiCut_run || eventNumber_physic!=physicNumber.dPhiCut_event) 
                  {
                    ++physicEvents.dPhiCut;
                    physicNumber.dPhiCut_run=runNumber_physic;
                    physicNumber.dPhiCut_event=eventNumber_physic;
                  }
                  h_ee_physic->Fill(eeMass_physic);
                  double r1=-1/(2*curv_e_physic);
                  double r2=1/(2*curv_p_physic);
                  double R=convert(r1, r2, d0_e_physic, d0_p_physic, phi_e, phi_p);
                  h_R_physic->Fill(R);
                  if (dsPlusCharge_physic<0) oppRecon_physic+=1;
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
  std::cout<<"Number of signal events after KsKppipi tag = "<<signalEvents.tagCut<<std::endl;
  std::cout<<"Number of signal events after dsPlusMCut = "<<signalEvents.dsPlusMCut<<std::endl;
  std::cout<<"Number of signal events after deltaECut = "<<signalEvents.deltaECut<<std::endl;
  std::cout<<"Number of signal events after mbcCut = "<<signalEvents.mbcCut<<std::endl;
  std::cout<<"Number of signal events after deltaMCut = "<<signalEvents.deltaMCut<<std::endl;
  std::cout<<"Number of signal events after diffD0Cut = "<<signalEvents.diffD0Cut<<std::endl;
  std::cout<<"Number of signal events after dPhiCut = "<<signalEvents.dPhiCut<<std::endl;
  std::cout<<"Number of D_s- reconstructions remaining in signal = "<<oppRecon_signal<<std::endl;
  
  std::cout<<"Number of conversion candidates = "<<nConverEvents<<std::endl;
  for (DecayMap::iterator i_decay=decayFrequency_conver.begin(); i_decay!=decayFrequency_conver.end(); ++i_decay)
  {
    std::cout<<"Conversion decay mode "<<i_decay->first<<" had "<<i_decay->second<<" decays"<<std::endl;
  }
  std::cout<<"Number of conver events after ntuplizer = "<<converEvents.noCut<<std::endl;
  std::cout<<"Number of conver events after KsKppipi tag = "<<converEvents.tagCut<<std::endl;
  std::cout<<"Number of conver events after dsPlusMCut = "<<converEvents.dsPlusMCut<<std::endl;
  std::cout<<"Number of conver events after deltaECut = "<<converEvents.deltaECut<<std::endl;
  std::cout<<"Number of conver events after mbcCut = "<<converEvents.mbcCut<<std::endl;
  std::cout<<"Number of conver events after deltaMCut = "<<converEvents.deltaMCut<<std::endl;
  std::cout<<"Number of conver events after diffD0Cut = "<<converEvents.diffD0Cut<<std::endl;
  std::cout<<"Number of conver events after dPhiCut = "<<converEvents.dPhiCut<<std::endl;
  std::cout<<"Number of D_s- reconstructions remaining in conver = "<<oppRecon_conver<<std::endl;
  
  std::cout<<"Number of data candidates = "<<nPhysicEvents<<std::endl;
  for (DecayMap::iterator i_decay=decayFrequency_physic.begin(); i_decay!=decayFrequency_physic.end(); ++i_decay)
  {
    std::cout<<"physicsion decay mode "<<i_decay->first<<" had "<<i_decay->second<<" decays"<<std::endl;
  }
  std::cout<<"Number of physic events after ntuplizer = "<<physicEvents.noCut<<std::endl;
  std::cout<<"Number of physic events after KsKppipi tag = "<<physicEvents.tagCut<<std::endl;
  std::cout<<"Number of physic events after dsPlusMCut = "<<physicEvents.dsPlusMCut<<std::endl;
  std::cout<<"Number of physic events after deltaECut = "<<physicEvents.deltaECut<<std::endl;
  std::cout<<"Number of physic events after mbcCut = "<<physicEvents.mbcCut<<std::endl;
  std::cout<<"Number of physic events after deltaMCut = "<<physicEvents.deltaMCut<<std::endl;
  std::cout<<"Number of physic events after diffD0Cut = "<<physicEvents.diffD0Cut<<std::endl;
  std::cout<<"Number of physic events after dPhiCut = "<<physicEvents.dPhiCut<<std::endl;
  std::cout<<"Number of D_s- reconstructions remaining in physic = "<<oppRecon_physic<<std::endl;
  
  TLine *line;
  float xmin, xmax, ymin, ymax;
  
  TCanvas *dsPlusM = new TCanvas("dsPlusM");
  xmin=dsPlusMCut_center-dsPlusMCut_range;
  xmax=dsPlusMCut_center+dsPlusMCut_range;
  dsPlusM->Divide(1,3);
  dsPlusM->cd(1);
  h_dsPlusM_signal->Draw("SAME");
  line = new TLine(xmin,500,xmin,0); line->Draw();
  line = new TLine(xmax,500,xmax,0); line->Draw();  
  dsPlusM->cd(2);
  h_dsPlusM_conver->Draw("SAME");
  line = new TLine(xmin,500,xmin,0); line->Draw();
  line = new TLine(xmax,500,xmax,0); line->Draw(); 
  dsPlusM->cd(3);
  h_dsPlusM_physic->Draw("SAME");
  line = new TLine(xmin,70,xmin,0); line->Draw();
  line = new TLine(xmax,70,xmax,0); line->Draw(); 
  
  TCanvas *DeltaE_MBC = new TCanvas("DeltaE_MBC");
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
  
  TCanvas *DeltaE_DeltaM = new TCanvas("DeltaE_DeltaM");
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
  
  TCanvas *MBC_DeltaM = new TCanvas("MBC_DeltaM");
  xmin=mbcCut_center-mbcCut_range; xmax=mbcCut_center+mbcCut_range;
  ymin=deltaMCut_center-deltaMCut_range; ymax=deltaMCut_center+deltaMCut_range;
  MBC_DeltaM->Divide(1,3);
  MBC_DeltaM->cd(1);
  h_MBC_DeltaM_signal->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  MBC_DeltaM->cd(2);
  h_MBC_DeltaM_conver->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  MBC_DeltaM->cd(3);
  h_MBC_DeltaM_physic->Draw("box");
  line = new TLine(xmin, ymax, xmax, ymax); line->Draw();
  line = new TLine(xmax, ymax, xmax, ymin); line->Draw();
  line = new TLine(xmax, ymin, xmin, ymin); line->Draw();
  line = new TLine(xmin, ymin, xmin, ymax); line->Draw();
  
  TCanvas *DeltaE = new TCanvas("DeltaE");
  xmin=-deltaECut_range; xmax=deltaECut_range;
  DeltaE->Divide(1,3);
  DeltaE->cd(1);
  h_DeltaE_signal->Draw("SAME");
  line = new TLine(xmin,150,xmin,0); line->Draw();
  line = new TLine(xmax,150,xmax,0); line->Draw();
  DeltaE->cd(2);
  h_DeltaE_conver->Draw("SAME");
  line = new TLine(xmin,25,xmin,0); line->Draw();
  line = new TLine(xmax,25,xmax,0); line->Draw();
  DeltaE->cd(3);
  h_DeltaE_physic->Draw("SAME");
  line = new TLine(xmin,4,xmin,0); line->Draw();
  line = new TLine(xmax,4,xmax,0); line->Draw();
  
  TCanvas *MBC = new TCanvas("MBC");
  xmin=mbcCut_center-mbcCut_range; xmax=mbcCut_center+mbcCut_range;
  MBC->Divide(1,3);
  MBC->cd(1);
  h_MBC_signal->Draw("SAME");
  line = new TLine(xmin,250,xmin,0); line->Draw();
  line = new TLine(xmax,250,xmax,0); line->Draw();
  MBC->cd(2);
  h_MBC_conver->Draw("SAME");
  line = new TLine(xmin,15,xmin,0); line->Draw();
  line = new TLine(xmax,15,xmax,0); line->Draw();
  MBC->cd(3);
  h_MBC_physic->Draw("SAME");
  line = new TLine(xmin,3,xmin,0); line->Draw();
  line = new TLine(xmax,3,xmax,0); line->Draw();
  
  TCanvas *DeltaM = new TCanvas("DeltaM");
  xmin=deltaMCut_center-deltaMCut_range; xmax=deltaMCut_center+deltaMCut_range;
  DeltaM->Divide(1,3);
  DeltaM->cd(1);
  h_DeltaM_signal->Draw("SAME");
  line = new TLine(xmin,200,xmin,0); line->Draw();
  line = new TLine(xmax,200,xmax,0); line->Draw();
  DeltaM->cd(2);
  h_DeltaM_conver->Draw("SAME");
  line = new TLine(xmin,10,xmin,0); line->Draw();
  line = new TLine(xmax,10,xmax,0); line->Draw();
  DeltaM->cd(3);
  h_DeltaM_physic->Draw("SAME");
  line = new TLine(xmin,1,xmin,0); line->Draw();
  line = new TLine(xmax,1,xmax,0); line->Draw();
  
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
  line = new TLine(diffD0Cut,60,diffD0Cut,0); line->Draw();
  diffD0->cd(2);
  h_diffD0_conver->Draw("SAME");
  line = new TLine(diffD0Cut,5,diffD0Cut,0); line->Draw();
  diffD0->cd(3);
  h_diffD0_physic->Draw("SAME");
  line = new TLine(diffD0Cut,2,diffD0Cut,0); line->Draw();
  
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
  
  TCanvas *dPhi = new TCanvas("dPhi");
  dPhi->Divide(1,3);
  dPhi->cd(1);
  h_dPhi_signal->Draw("SAME");
  line = new TLine(dPhiCutLess,150,dPhiCutLess,0); line->Draw();
  dPhi->cd(2);
  h_dPhi_conver->Draw("SAME");
  line = new TLine(dPhiCutLess,10,dPhiCutLess,0); line->Draw();
  dPhi->cd(3);
  h_dPhi_physic->Draw("SAME");
  line = new TLine(dPhiCutLess,2,dPhiCutLess,0); line->Draw();
  
  TCanvas *dPhi_ee = new TCanvas("dPhi_ee");
  dPhi_ee->Divide(1,2);
  dPhi_ee->cd(1);
  h_dPhi_ee_signal->Draw();
  dPhi_ee->cd(2);
  h_dPhi_ee_conver->Draw();
  
  TCanvas *dPhi_sumPhi = new TCanvas("dPhi_sumPhi");
  dPhi_sumPhi->Divide(1,2);
  dPhi_sumPhi->cd(1);
  h_dPhi_sumPhi_signal->Draw();
  dPhi_sumPhi->cd(2);
  h_dPhi_sumPhi_conver->Draw();
  
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
  
  
  
  return 0;
}
