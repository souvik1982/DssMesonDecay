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

int CompareElectronPionFits()
{

  //TFile *file = new TFile("/nfs/cor/an2/souvik/Dataset48/DsTaggedDecaysProc_ReTaggedData_232264_234607.root");
  TFile *file = new TFile("/nfs/cor/an2/souvik/Dataset47/DsTaggedDecaysProc_ReTaggedData_230474_232255.root");
  file->cd("DsTaggedDecaysProc");
  TTree *tree = (TTree*)gDirectory->Get("nt5");
  
  float run, event, trackNumber;
  float e_electron, px_electron, py_electron, pz_electron;
  float e_pion, px_pion, py_pion, pz_pion;
  
  tree->SetBranchAddress("Run", &(run));
  tree->SetBranchAddress("Event", &(event));
  tree->SetBranchAddress("trackNumber", &(trackNumber));
  tree->SetBranchAddress("trackE_electron", &(e_electron));
  tree->SetBranchAddress("trackPx_electron", &(px_electron));
  tree->SetBranchAddress("trackPy_electron", &(py_electron));
  tree->SetBranchAddress("trackPz_electron", &(pz_electron));
  tree->SetBranchAddress("trackE_pion", &(e_pion));
  tree->SetBranchAddress("trackPx_pion", &(px_pion));
  tree->SetBranchAddress("trackPy_pion", &(py_pion));
  tree->SetBranchAddress("trackPz_pion", &(pz_pion));
  
  TH2F *h_dE = new TH2F("h_dE", "Electron Energy Deviation between Fit Hypotheses; E(e-fit track) GeV; E(pi-fit track)-E(e-fit track) GeV", 100, 0., .16, 100, -0.02, 0.02);
  
  std::cout<<tree->GetEntries()<<std::endl;
  for (int i=0; i<tree->GetEntries(); ++i)
  {
    tree->GetEvent(i);
    if (trackNumber!=-1)
    {
      //std::cout<<"e_pion = "<<e_pion<<", e_electron = "<<e_electron<<std::endl;
      h_dE->Fill(e_electron, e_pion-e_electron);
    }
  }
  
  gROOT->SetStyle("Plain");
  
  TCanvas *dE = new TCanvas("dE");
  dE->Divide(1,2);
  dE->cd(1);
  h_dE->Draw("box");
  dE->cd(2);
  TProfile *h_dE_profile=h_dE->ProfileX();
  h_dE_profile->Draw();
  
  return 0;
}
