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

double pi=3.14159265358979;

double deltaPhi(double phi1, double phi2)
{
  double dphi=phi1-phi2;
  if (dphi<-pi) dphi=2*pi+dphi;
  if (dphi>pi) dphi=2*pi-dphi;
  return dphi;
}

bool phiMCut(float phiM)
{
  return (fabs(phiM-1.0195)<0.015);
}

bool dsPlusMCut(float dsPlusM)
{
  return (fabs(dsPlusM-1.96849)<0.02);
}

bool DeltaECut(float DeltaE)
{
  return (fabs(DeltaE)<0.016);
}

bool MBCCut(float MBC)
{
  return (fabs(MBC-2.112)<0.005);
}

bool DeltaMCut(float DeltaM)
{
  return (fabs(DeltaM-0.144)<0.005);
}

bool dPhiCut(float dPhi)
{
  return (dPhi<0.1);
}

double convert(double r1, double r2, double d1, double d2, double phi1, double phi2)
{
  double a1=r1+d1;
  double a2=r2+d2;
  double theta=phi2-phi1;
  double R=a1*a2*sin(theta)/pow(a1*a1+a2*a2+2*a1*a2*cos(theta), 0.5);
  return R;
}

int DsAnalysisDiffD0Cut()
{

  TFile *signalFile = new TFile("MyDChainFile_MC_vtosll.root");
  signalFile->cd("MyDChainProc");
  TTree *signalTree = (TTree*)gDirectory->Get("nt5");
  
  TFile *converFile = new TFile("MyDChainFile_MCgamma_eRefit.root");
  converFile->cd("MyDChainProc");
  TTree *converTree = (TTree*)gDirectory->Get("nt5");  
  
  float phiM_signal, dsPlusM_signal, DeltaE_signal, MBC_signal, DeltaM_signal, eeMass_signal;
  float phiM_conver, dsPlusM_conver, DeltaE_conver, MBC_conver, DeltaM_conver, eeMass_conver;
  float d0_e_signal, d0_p_signal, px_e_signal, py_e_signal, px_p_signal, py_p_signal, E_e_signal, E_p_signal, curv_e_signal, curv_p_signal;
  float d0_e_conver, d0_p_conver, px_e_conver, py_e_conver, px_p_conver, py_p_conver, E_e_conver, E_p_conver, curv_e_conver, curv_p_conver;
  
  TH1D *h_phiM_signal = new TH1D("h_phiM_signal", "h_phiM_signal", 100, .95, 1.15); h_phiM_signal->SetLineColor(kRed);
  TH1D *h_phiM_conver = new TH1D("h_phiM_conver", "h_phiM_conver", 100, .95, 1.15); h_phiM_conver->SetLineColor(kBlue);
  TH1D *h_dsPlusM_signal = new TH1D("h_dsPlusM_signal", "h_dsPlusM_signal", 100, 1.9, 2.1); h_dsPlusM_signal->SetLineColor(kRed);
  TH1D *h_dsPlusM_conver = new TH1D("h_dsPlusM_conver", "h_dsPlusM_conver", 100, 1.9, 2.1); h_dsPlusM_conver->SetLineColor(kBlue);
  TH1D *h_DeltaE_signal = new TH1D("h_DeltaE_signal", "DeltaE", 100, -0.1, 0.2); h_DeltaE_signal->SetLineColor(kRed);
  TH1D *h_DeltaE_conver = new TH1D("h_DeltaE_conver", "DeltaE", 100, -0.1, 0.2); h_DeltaE_conver->SetLineColor(kBlue);
  TH1D *h_MBC_signal = new TH1D("h_MBC_signal", "MBC", 100, 2., 2.2); h_MBC_signal->SetLineColor(kRed);
  TH1D *h_MBC_conver = new TH1D("h_MBC_conver", "MBC", 100, 2., 2.2); h_MBC_conver->SetLineColor(kBlue);
  TH1D *h_DeltaM_signal = new TH1D("h_DeltaM_signal", "DeltaM", 100, 0.0, 0.2); h_DeltaM_signal->SetLineColor(kRed);
  TH1D *h_DeltaM_conver = new TH1D("h_DeltaM_conver", "DeltaM", 100, 0.0, 0.2); h_DeltaM_conver->SetLineColor(kBlue);
  TH1D *h_d0_e_signal = new TH1D("h_d0_e_signal", "d0_e", 50, -0.01, 0.01); h_d0_e_signal->SetLineColor(kRed);
  TH1D *h_d0_e_conver = new TH1D("h_d0_e_conver", "d0_e", 50, -0.01, 0.01); h_d0_e_conver->SetLineColor(kBlue);
  TH1D *h_d0_p_signal = new TH1D("h_d0_p_signal", "d0_p", 50, -0.01, 0.01); h_d0_p_signal->SetLineColor(kRed);
  TH1D *h_d0_p_conver = new TH1D("h_d0_p_conver", "d0_p", 50, -0.01, 0.01); h_d0_p_conver->SetLineColor(kBlue);
  TH1D *h_sumD0_signal = new TH1D("h_sumD0_signal", "h_sumD0_signal", 50, -0.01, 0.01); h_sumD0_signal->SetLineColor(kRed);
  TH1D *h_sumD0_conver = new TH1D("h_sumD0_conver", "h_sumD0_conver", 50, -0.01, 0.01); h_sumD0_conver->SetLineColor(kBlue);
  TH1D *h_diffD0_signal = new TH1D("h_diffD0_signal", "h_diffD0_signal", 50, -0.01, 0.01); h_diffD0_signal->SetLineColor(kRed);
  TH1D *h_diffD0_conver = new TH1D("h_diffD0_conver", "h_diffD0_conver", 50, -0.01, 0.01); h_diffD0_conver->SetLineColor(kRed);
  TH1D *h_dPhi_signal = new TH1D("h_dPhi_signal", "dPhi", 50, -2., 2.); h_dPhi_signal->SetLineColor(kRed);
  TH1D *h_dPhi_conver = new TH1D("h_dPhi_conver", "dPhi", 50, -2., 2.); h_dPhi_conver->SetLineColor(kBlue);
  TH2D *h_dPhi_ee_signal = new TH2D("h_dPhi_ee_signal", "dPhi_ee_signal", 50, 0., 0.2, 50, -2., 2.);
  TH2D *h_dPhi_ee_conver = new TH2D("h_dPhi_ee_conver", "dPhi_ee_conver", 50, 0., 0.2, 50, -2., 2.);
  TH2D *h_d0_phi_e_signal = new TH2D("h_d0_phi_e_signal", "h_d0_phi_e_signal", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_p_signal = new TH2D("h_d0_phi_p_signal", "h_d0_phi_p_signal", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_e_conver = new TH2D("h_d0_phi_e_conver", "h_d0_phi_e_conver", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_d0_phi_p_conver = new TH2D("h_d0_phi_p_signal", "h_d0_phi_p_signal", 50, -pi, pi, 50, -0.01, 0.01);
  TH2D *h_dPhi_sumPhi_signal = new TH2D("h_dPhi_sumPhi_signal", "dPhi_sumPhi_signal", 50, -2*pi, 2*pi, 50, -2., 2.);
  TH2D *h_dPhi_sumPhi_conver = new TH2D("h_dPhi_sumPhi_conver", "dPhi_sumPhi_conver", 50, -2*pi, 2*pi, 50, -2., 2.);
  TH2D *h_dPhi_diffD0_signal = new TH2D("h_dPhi_diffD0_signal", "h_dPhi_diffD0_signal", 50, -2., 2., 50, -0.01, 0.01);
  TH2D *h_dPhi_diffD0_conver = new TH2D("h_dPhi_diffD0_conver", "h_dPhi_diffD0_conver", 50, -2., 2., 50, -0.01, 0.01);  
  TH1D *h_EnergyAsymmetry_signal = new TH1D("h_EnergyAsymmetry_signal", "h_EnergyAsymmetry_signal", 50, 0., .2); h_EnergyAsymmetry_signal->SetLineColor(kRed);
  TH1D *h_EnergyAsymmetry_conver = new TH1D("h_EnergyAsymmetry_conver", "h_EnergyAsymmetry_conver", 50, 0., .2); h_EnergyAsymmetry_conver->SetLineColor(kBlue);
  TH1D *h_R_signal = new TH1D("h_R_signal", "h_R_signal", 50, -.4, .4); h_R_signal->SetLineColor(kRed);
  TH1D *h_R_conver = new TH1D("h_R_conver", "h_R_conver", 50, -.4, .4); h_R_conver->SetLineColor(kBlue);
  TH1D *h_ee_signal = new TH1D("h_ee_signal", "h_ee_signal", 20, 0.0, 0.2); h_ee_signal->SetLineColor(kRed);
  TH1D *h_ee_conver = new TH1D("h_ee_conver", "h_ee_conver", 20, 0.0, 0.2); h_ee_conver->SetLineColor(kBlue); 
  
  signalTree->SetBranchAddress("phiM", &(phiM_signal));
  signalTree->SetBranchAddress("dsPlusM", &(dsPlusM_signal));
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
  
  converTree->SetBranchAddress("phiM", &(phiM_conver));
  converTree->SetBranchAddress("dsPlusM", &(dsPlusM_conver));
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
  
  int nSignalEvents=signalTree->GetEntries();
  int nConverEvents=converTree->GetEntries();
  
  std::cout<<"Number of signal events = "<<nSignalEvents<<std::endl;
  std::cout<<"Number of conversion events = "<<nConverEvents<<std::endl;
  
  for (int i=0; i<nSignalEvents; ++i)
  {
    signalTree->GetEvent(i);
    
    double phi_e=atan2(py_e_signal, px_e_signal);
    double phi_p=atan2(py_p_signal, px_p_signal);
    double dPhi=deltaPhi(phi_e, phi_p);
    double sumPhi=phi_e+phi_p;
    
    h_phiM_signal->Fill(phiM_signal);
    
    if (phiMCut(phiM_signal))
    {
      h_dsPlusM_signal->Fill(dsPlusM_signal);
      if (dsPlusMCut(dsPlusM_signal))
      {
        h_DeltaE_signal->Fill(DeltaE_signal);
        if (DeltaECut(DeltaE_signal))
        {
          h_MBC_signal->Fill(MBC_signal);
          if (MBCCut(MBC_signal))
          {
            h_DeltaM_signal->Fill(DeltaM_signal);
            if (DeltaMCut(DeltaM_signal))
            {
              h_d0_e_signal->Fill(d0_e_signal);
              h_d0_p_signal->Fill(d0_p_signal);
              h_dPhi_diffD0_signal->Fill(dPhi, d0_e_signal-d0_p_signal);
              h_d0_phi_e_signal->Fill(phi_e, d0_e_signal);
              h_d0_phi_p_signal->Fill(phi_p, d0_p_signal);
              h_dPhi_signal->Fill(dPhi);
              h_dPhi_sumPhi_signal->Fill(sumPhi, dPhi);
              h_dPhi_ee_signal->Fill(eeMass_signal, dPhi);
              if (dPhiCut(dPhi))
              {
                h_diffD0_signal->Fill(d0_e_signal-d0_p_signal);
                h_ee_signal->Fill(eeMass_signal);
                h_EnergyAsymmetry_signal->Fill(fabs(E_e_signal-E_p_signal));
                double r1=-1/(2*curv_e_signal);
                double r2=1/(2*curv_p_signal);
                double R=convert(r1, r2, d0_e_signal, d0_p_signal, phi_e, phi_p);
                h_R_signal->Fill(R);
                if (R<0.06)
                {
                  h_sumD0_signal->Fill(d0_e_signal+d0_p_signal);
                }
              }
            }
          }
        }
      }
    }    
  }
  
  for (int i=0; i<nConverEvents; ++i)
  {
    converTree->GetEvent(i);
    
    double phi_e=atan2(py_e_conver, px_e_conver);
    double phi_p=atan2(py_p_conver, px_p_conver);
    double dPhi=deltaPhi(phi_e, phi_p);
    double sumPhi=phi_e+phi_p;
    
    h_phiM_conver->Fill(phiM_conver);
    
    if (phiMCut(phiM_conver))
    {
      h_dsPlusM_conver->Fill(dsPlusM_conver);
      if (dsPlusMCut(dsPlusM_conver))
      {
        h_DeltaE_conver->Fill(DeltaE_conver);
        if (DeltaECut(DeltaE_conver))
        {
          h_MBC_conver->Fill(MBC_conver);
          if (MBCCut(MBC_conver))
          {
            h_DeltaM_conver->Fill(DeltaM_conver);
            if (DeltaMCut(DeltaM_conver))
            {
              h_d0_e_conver->Fill(d0_e_conver);
              h_d0_p_conver->Fill(d0_p_conver);
              h_dPhi_diffD0_conver->Fill(dPhi, d0_e_conver-d0_p_conver);
              h_d0_phi_e_conver->Fill(phi_e, d0_e_conver);
              h_d0_phi_p_conver->Fill(phi_p, d0_p_conver);
              h_dPhi_conver->Fill(dPhi);
              h_dPhi_sumPhi_conver->Fill(sumPhi, dPhi);
              h_dPhi_ee_conver->Fill(eeMass_conver, dPhi);
              if (dPhiCut(dPhi))
              {
                h_diffD0_conver->Fill(d0_e_conver-d0_p_conver);
                h_ee_conver->Fill(eeMass_conver);
                h_EnergyAsymmetry_conver->Fill(fabs(E_e_conver-E_p_conver));
                double r1=-1/(2*curv_e_conver);
                double r2=1/(2*curv_p_conver);
                double R=convert(r1, r2, d0_e_conver, d0_p_conver, phi_e, phi_p);
                h_R_conver->Fill(R);
                if (R<0.06)
                {
                  h_sumD0_conver->Fill(d0_e_signal+d0_p_conver);
                }
              }
            }
          }
        }
      }
    } 
  }
  
  TCanvas *phiM = new TCanvas("phiM");
  //h_phiM_signal->Scale(1/h_phiM_signal->GetEntries());
  //h_phiM_conver->Scale(1/h_phiM_conver->GetEntries());
  phiM->Divide(1,2);
  phiM->cd(1);
  h_phiM_signal->Draw("SAME");
  phiM->cd(2);
  h_phiM_conver->Draw("SAME");
  //TLine a(0.2,0.2,0.8,0.8);
  //a.Draw();
  //TPaveStats *s=(TPaveStats*)h_phiM_conver->GetListOfFunctions()->FindObject("stats");
  //s->SetX1NDC(0.5); s->SetX2NDC(0.7);
  
  
  TCanvas *dsPlusM = new TCanvas("dsPlusM");
  //h_dsPlusM_signal->Scale(1/h_dsPlusM_signal->GetEntries());
  //h_dsPlusM_conver->Scale(1/h_dsPlusM_conver->GetEntries());
  dsPlusM->Divide(1,2);
  dsPlusM->cd(1);
  h_dsPlusM_signal->Draw("SAME");
  dsPlusM->cd(2);
  h_dsPlusM_conver->Draw("SAME");
  
  TCanvas *DeltaE = new TCanvas("DeltaE");
  //h_DeltaE_signal->Scale(1/h_DeltaE_signal->GetEntries());
  //h_DeltaE_conver->Scale(1/h_DeltaE_conver->GetEntries());
  DeltaE->Divide(1,2);
  DeltaE->cd(1);
  h_DeltaE_signal->Draw("SAME");
  DeltaE->cd(2);
  h_DeltaE_conver->Draw("SAME");
  
  TCanvas *MBC = new TCanvas("MBC");
  //h_MBC_signal->Scale(1/h_MBC_signal->GetEntries());
  //h_MBC_conver->Scale(1/h_MBC_conver->GetEntries());
  MBC->Divide(1,2);
  MBC->cd(1);
  h_MBC_signal->Draw("SAME");
  MBC->cd(2);
  h_MBC_conver->Draw("SAME");
  
  TCanvas *DeltaM = new TCanvas("DeltaM");
  //h_DeltaM_signal->Scale(1/h_DeltaM_signal->GetEntries());
  //h_DeltaM_conver->Scale(1/h_DeltaM_conver->GetEntries());
  DeltaM->Divide(1,2);
  DeltaM->cd(1);
  h_DeltaM_signal->Draw("SAME");
  DeltaM->cd(2);
  h_DeltaM_conver->Draw("SAME");
  
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
  
  TCanvas *sumD0 = new TCanvas("sumD0");
  //h_sumD0_signal->Scale(1/h_sumD0_signal->GetEntries());
  //h_sumD0_conver->Scale(1/h_sumD0_conver->GetEntries());
  sumD0->Divide(1,2);
  sumD0->cd(1);
  h_sumD0_signal->Draw("SAME");
  sumD0->cd(2);
  h_sumD0_conver->Draw("SAME");
  
  TCanvas *diffD0 = new TCanvas("diffD0");
  diffD0->Divide(1,2);
  diffD0->cd(1);
  h_diffD0_signal->Draw("SAME");
  diffD0->cd(2);
  h_diffD0_conver->Draw("SAME");
  
  TCanvas *dPhi_diffD0 = new TCanvas ("dPhi_diffD0");
  dPhi_diffD0->Divide(1,2);
  dPhi_diffD0->cd(1);
  h_dPhi_diffD0_signal->Draw();
  dPhi_diffD0->cd(2);
  h_dPhi_diffD0_conver->Draw();
  
  TCanvas *dPhi = new TCanvas("dPhi");
  //h_dPhi_signal->Scale(1/h_dPhi_signal->GetEntries());
  //h_dPhi_conver->Scale(1/h_dPhi_conver->GetEntries());
  dPhi->Divide(1,2);
  dPhi->cd(1);
  h_dPhi_signal->Draw("SAME");
  dPhi->cd(2);
  h_dPhi_conver->Draw("SAME");
  
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
  
  TCanvas *EnergyAsymmetry = new TCanvas("EnergyAsymmetry");
  //h_EnergyAsymmetry_signal->Scale(1/h_EnergyAsymmetry_signal->GetEntries());
  //h_EnergyAsymmetry_conver->Scale(1/h_EnergyAsymmetry_conver->GetEntries());
  EnergyAsymmetry->Divide(1,2);
  EnergyAsymmetry->cd(1);
  h_EnergyAsymmetry_signal->Draw("SAME");
  EnergyAsymmetry->cd(2);
  h_EnergyAsymmetry_conver->Draw("SAME");
  
  TCanvas *R = new TCanvas("R");
  R->Divide(1,2);
  R->cd(1);
  h_R_signal->Draw();
  R->cd(2);
  h_R_conver->Draw();
  
  TCanvas *ee_Mass = new TCanvas("ee_Mass");
  ee_Mass->Divide(1,2);
  ee_Mass->cd(1);
  h_ee_signal->Draw();
  ee_Mass->cd(2);
  h_ee_conver->Draw();
  
  
  return 0;
}
