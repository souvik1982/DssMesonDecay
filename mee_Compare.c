#include "TCanvas.h"
#include "TH1D.h"
#include "TROOT.h"
#include "THStack.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

int mee_Compare()
{
  const double m_e=0.000511;
  const double m_Ds=1.969;
  const double m_Dsstr=2.112;
  const double qmin=2*m_e;
  const double qmax=(m_Dsstr-m_Ds);
  
  double xmin=log(qmin*1000)/log(10);
  double xmax=log(qmax*1000)/log(10);

  TFile *file=new TFile("sum_mee.root");
  TH1D *h_mee_signal = (TH1D*)gDirectory->Get("h_mee_signal");
  TH1D *h_mee_conver_Dsp = (TH1D*)gDirectory->Get("h_mee_conver_Dsp");
  TH1D *h_mee_conver_Dsm = (TH1D*)gDirectory->Get("h_mee_conver_Dsm");
  TH1D *h_mee_generic = (TH1D*)gDirectory->Get("h_mee_generic");
  TH1D *h_mee_continu = (TH1D*)gDirectory->Get("h_mee_continu");
  TH1D *h_mee_physics = (TH1D*)gDirectory->Get("h_mee_physics"); h_mee_physics->SetLineColor(kBlack); h_mee_physics->Sumw2();
  
  TH1D *h_mee_conver = (TH1D*)h_mee_conver_Dsp->Clone();
  h_mee_conver->Add(h_mee_conver_Dsm);
  
  h_mee_conver_Dsp->SetFillColor(kRed);  h_mee_conver_Dsp->SetFillStyle(3244); h_mee_conver_Dsp->SetStats(false);
  h_mee_conver_Dsm->SetFillColor(kRed); h_mee_conver_Dsm->SetFillStyle(3244); h_mee_conver_Dsm->SetStats(false);
  h_mee_conver->SetFillColor(kRed); h_mee_conver->SetFillStyle(3244); h_mee_conver->SetStats(false);
  h_mee_generic->SetFillColor(kGreen); h_mee_generic->SetFillStyle(3245); h_mee_generic->SetStats(false);
  h_mee_continu->SetFillColor(kBlue); h_mee_continu->SetFillStyle(3003); h_mee_continu->SetStats(false);
  h_mee_signal->SetFillColor(kCyan); h_mee_signal->SetFillStyle(3016); h_mee_signal->SetStats(false);
  
  THStack *s_mee = new THStack("s_mee", "; log_{10}(m_{ee} / 1 MeV); # Events");
  s_mee->Add(h_mee_continu, "hist");
  s_mee->Add(h_mee_generic, "hist");
  s_mee->Add(h_mee_conver, "hist");
  s_mee->Add(h_mee_signal, "hist");
  
  gROOT->SetStyle("Plain");
  gStyle->SetErrorX(0);
  
  TArrow *line;
  
  TCanvas *c_mee = new TCanvas("c_mee", "c_mee");
  s_mee->SetMaximum(15);
  s_mee->Draw();
  h_mee_physics->SetMarkerStyle(20);
  h_mee_physics->SetMarkerSize(0.9);
  h_mee_physics->Draw("SAMEE1");
  //line = new TArrow(xmin, 12, xmin, 0.5, 0.02, ">"); line->Draw();
  //line = new TArrow(xmax, 12, xmax, 0.5, 0.02, ">"); line->Draw();
  TLegend *legend=new TLegend(.6, .9, .9, .7);
  legend->AddEntry(h_mee_signal, "Signal Simulation", "f");
  legend->AddEntry(h_mee_continu, "non-c#bar{c} Simulation", "f");
  legend->AddEntry(h_mee_generic, "c#bar{c} w/o Conversions Simulation", "f");
  legend->AddEntry(h_mee_conver_Dsp, "Conversions Simulation", "f");
  legend->AddEntry(h_mee_physics, "Data");
  legend->SetFillColor(kWhite);
  legend->Draw();
  
  // Kolmogorov-Smirnov Test
  TH1D *h_mee_sum = new TH1D(*h_mee_signal);
  h_mee_sum->Add(h_mee_continu);
  h_mee_sum->Add(h_mee_generic);
  h_mee_sum->Add(h_mee_conver_Dsm);
  h_mee_sum->Add(h_mee_conver_Dsp);
  h_mee_sum->Sumw2();
  double ksProb=h_mee_sum->KolmogorovTest(h_mee_physics, "D");
  std::cout<<"ksProb = "<<ksProb<<std::endl;
  
  TCanvas *c_check = new TCanvas("c_check", "c_check");
  h_mee_sum->SetMaximum(18);
  h_mee_sum->Draw();
  h_mee_physics->Draw("SAME");
  
}
