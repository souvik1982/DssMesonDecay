#include "TCanvas.h"
#include "TH1D.h"
#include "TROOT.h"
#include "THStack.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLine.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

int eEnergy_Compare()
{

  TFile *file_MC=new TFile("eEnergy_MC.root");
  TH1D *h_eEnergy_MC = (TH1D*)gDirectory->Get("h_e_Reco_e"); h_eEnergy_MC->Sumw2();
  
  TFile *file_Reco=new TFile("eEnergy_data.root");
  TH1D *h_eEnergy_data = (TH1D*)gDirectory->Get("h_e_Reco_e"); h_eEnergy_data->Sumw2();
  
  TH1D *h_eEnergy_Ratio = new TH1D("h_eEnergy_Ratio", "h_eEnergy_Ratio", 100, 0., 0.3);
  
  h_eEnergy_Ratio->Divide(h_eEnergy_data, h_eEnergy_MC);
  
  TCanvas *c_eEnergy = new TCanvas("c_eEnergy", "c_eEnergy");
  c_eEnergy->Divide(1,3);
  c_eEnergy->cd(1);
  h_eEnergy_MC->Draw();
  c_eEnergy->cd(2);
  h_eEnergy_data->Draw();
  c_eEnergy->cd(3);
  h_eEnergy_Ratio->Draw();
  TF1 *f_line = new TF1("f_line", "[0]+[1]*x", 0.02, 0.14);
  h_eEnergy_Ratio->Fit(f_line, "R");
  
}

