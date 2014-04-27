#include "TCanvas.h"
#include "TH1D.h"
#include "TROOT.h"
#include "THStack.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

double deltaMCut_center=0.1438;
double deltaMCut_range=0.007;
double sideband_dist=0.5;

double plotScale=200;

double T_1(double x)
{
  return x;
}

double T_2(double x)
{
  return 2*x*x-1;
}

double T_3(double x)
{
  return 4*x*x*x-3*x;
}

Double_t sum_background(Double_t *x, Double_t *par)
{
  return par[0]+par[1]*T_1(x[0])+par[2]*T_2(x[0])+par[3]*T_3(x[0]);
}

Double_t physics(Double_t *x, Double_t *par)
{
  if (x[0]>0.1298 && x[0]<0.1578)
  //if (x[0]>0.1368 && x[0]<0.1508)
  {
    TF1::RejectPoint();
    //return 0;
  }
  return par[0]+par[1]*T_1(x[0])+par[2]*T_2(x[0])+par[3]*T_3(x[0]);
}

Double_t physics_scaledMC(Double_t *x, Double_t *par)
{
  if (x[0]>0.1298 && x[0]<0.1578) TF1::RejectPoint();
  return par[4]*(par[0]+par[1]*T_1(x[0])+par[2]*T_2(x[0])+par[3]*T_3(x[0]));
}

int SidebandFitter_DeltaM()
{
  
  // TFile *file=new TFile("Sum_Sideband_DeltaM_005_mBCDoubled.root");
  TFile *file=new TFile("Sum_Unblinded_DeltaM_005.root");
  
  TH1D *h_DeltaM_signal = (TH1D*)gDirectory->Get("h_DeltaM_signal");
  TH1D *h_DeltaM_conver_Dsp = (TH1D*)gDirectory->Get("h_DeltaM_conver_Dsp");
  TH1D *h_DeltaM_conver_Dsm = (TH1D*)gDirectory->Get("h_DeltaM_conver_Dsm");
  TH1D *h_DeltaM_generic = (TH1D*)gDirectory->Get("h_DeltaM_generic");
  TH1D *h_DeltaM_continu = (TH1D*)gDirectory->Get("h_DeltaM_continu");
  TH1D *h_DeltaM_physics = (TH1D*)gDirectory->Get("h_DeltaM_physics"); h_DeltaM_physics->SetStats(false); h_DeltaM_physics->SetLineColor(kBlack);
  
  TH1D *h_DeltaM_conver = (TH1D*)h_DeltaM_conver_Dsp->Clone();
  h_DeltaM_conver->Add(h_DeltaM_conver_Dsm);
  
  TH1D *h_Sum_background = (TH1D*)h_DeltaM_conver_Dsp->Clone(); h_Sum_background->SetStats(false);
  h_Sum_background->SetLineColor(kMagenta);
  h_Sum_background->Add(h_DeltaM_conver_Dsm);
  h_Sum_background->Add(h_DeltaM_generic);
  h_Sum_background->Add(h_DeltaM_continu);
  
  h_DeltaM_conver_Dsp->SetFillColor(kRed);  h_DeltaM_conver_Dsp->SetFillStyle(3244); h_DeltaM_conver_Dsp->SetStats(false);
  h_DeltaM_conver_Dsm->SetFillColor(kRed); h_DeltaM_conver_Dsm->SetFillStyle(3244); h_DeltaM_conver_Dsm->SetStats(false);
  h_DeltaM_conver->SetFillColor(kRed); h_DeltaM_conver->SetFillStyle(3244); h_DeltaM_conver->SetStats(false);
  h_DeltaM_generic->SetFillColor(kGreen); h_DeltaM_generic->SetFillStyle(3245); h_DeltaM_generic->SetStats(false);
  h_DeltaM_continu->SetFillColor(kBlue); h_DeltaM_continu->SetFillStyle(3003); h_DeltaM_continu->SetStats(false);
  h_DeltaM_signal->SetFillColor(kCyan); h_DeltaM_signal->SetFillStyle(3016); h_DeltaM_signal->SetStats(false);
  
  double lowerSide_xmin=0.1;
  double lowerSide_xmax=deltaMCut_center-deltaMCut_range-sideband_dist*2*deltaMCut_range;
  double xmin=deltaMCut_center-deltaMCut_range;
  double xmax=deltaMCut_center+deltaMCut_range;
  double upperSide_xmin=deltaMCut_center+deltaMCut_range+sideband_dist*2*deltaMCut_range;
  double upperSide_xmax=0.25;  
  double ymax=(h_DeltaM_signal->GetMaximum())*0.95;
  TLine *line;
  
  gROOT->SetStyle("Plain");
  
  TCanvas *c_FitBackground = new TCanvas("c_FitBackground");
  h_Sum_background->Draw();
  TF1 *f_Sum_background=new TF1("f_Sum_background", sum_background, lowerSide_xmin, upperSide_xmax, 4);
  f_Sum_background->SetLineColor(kBlack);
  h_Sum_background->Fit("f_Sum_background", "RLLEFM");
  line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
  line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
  line=new TLine(xmin, ymax, xmin, 0); line->Draw();
  line=new TLine(xmax, ymax, xmax, 0); line->Draw();
  line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
  line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
  double total_bg=f_Sum_background->Integral(lowerSide_xmin, upperSide_xmax)*plotScale;
  double total_bg_count=h_Sum_background->Integral(h_Sum_background->FindBin(lowerSide_xmin), h_Sum_background->FindBin(upperSide_xmax));
  double s1_bg=f_Sum_background->Integral(lowerSide_xmin, lowerSide_xmax)*plotScale;
  double s2_bg=f_Sum_background->Integral(xmin, xmax)*plotScale;
  double s3_bg=f_Sum_background->Integral(upperSide_xmin, upperSide_xmax)*plotScale;
  double s1_bg_error=pow(s1_bg*(1-s1_bg/total_bg), 0.5);
  double s2_bg_error=pow(s2_bg*(1-s2_bg/total_bg), 0.5);
  double s3_bg_error=pow(s3_bg*(1-s3_bg/total_bg), 0.5);
  double R_bg=s2_bg/(s1_bg+s3_bg);
  double R_bg_error=R_bg*pow(pow(s2_bg_error/s2_bg, 2)+(pow(s1_bg_error, 2)+pow(s3_bg_error, 2))/pow(s1_bg+s3_bg, 2), 0.5);
  std::cout<<"total_bg_count = "<<total_bg_count<<std::endl;
  std::cout<<"total_bg = "<<total_bg<<std::endl;
  std::cout<<"s1_bg = "<<s1_bg<<" \pm "<<s1_bg_error<<std::endl;
  std::cout<<"s2_bg = "<<s2_bg<<" \pm "<<s2_bg_error<<std::endl;
  std::cout<<"s3_bg = "<<s3_bg<<" \pm "<<s3_bg_error<<std::endl;
  std::cout<<"R_bg = "<<R_bg<<" \pm "<<R_bg_error<<std::endl;
  
  TCanvas *c_FitData = new TCanvas("c_FitData");
  
  // TF1 *f_physics=new TF1("f_physics", physics, lowerSide_xmin, upperSide_xmax, 4);
  TF1 *f_physics=new TF1("f_physics", physics_scaledMC, lowerSide_xmin, upperSide_xmax, 5);
  f_physics->FixParameter(0, f_Sum_background->GetParameter(0));
  f_physics->FixParameter(1, f_Sum_background->GetParameter(1));
  f_physics->FixParameter(2, f_Sum_background->GetParameter(2));
  f_physics->FixParameter(3, f_Sum_background->GetParameter(3));
  
  f_physics->SetLineColor(kBlack);
  f_physics->SetLineWidth(2);
  h_DeltaM_physics->Fit(f_physics, "RLLEFM");
  line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
  line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
  line=new TLine(xmin, ymax, xmin, 0); line->Draw();
  line=new TLine(xmax, ymax, xmax, 0); line->Draw();
  line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
  line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
  double total_physics=f_physics->Integral(lowerSide_xmin, upperSide_xmax)*plotScale;
  double s1_physics=f_physics->Integral(lowerSide_xmin, lowerSide_xmax)*plotScale;
  double s2_physics=f_physics->Integral(xmin, xmax)*plotScale;
  double s3_physics=f_physics->Integral(upperSide_xmin, upperSide_xmax)*plotScale;
  double s1_physics_error=pow(s1_physics*(1-s1_physics/total_physics), 0.5);
  double s2_physics_error=pow(s2_physics*(1-s2_physics/total_physics), 0.5);
  double s3_physics_error=pow(s3_physics*(1-s3_physics/total_physics), 0.5);
  double R_physics=s2_physics/(s1_physics+s3_physics);
  double R_physics_error=R_physics*pow(pow(s2_physics_error/s2_physics, 2)+(pow(s1_physics_error, 2)+pow(s3_physics_error, 2))/pow(s1_physics+s3_physics, 2), 0.5);
  std::cout<<"total_physics = "<<total_physics<<std::endl;
  std::cout<<"s1_physics = "<<s1_physics<<" \pm "<<s1_physics_error<<std::endl;
  std::cout<<"s2_physics = "<<s2_physics<<" \pm "<<s2_physics_error<<std::endl;
  std::cout<<"s3_physics = "<<s3_physics<<" \pm "<<s3_physics_error<<std::endl;
  std::cout<<"R_physics = "<<R_physics<<" \pm "<<R_physics_error<<std::endl;
  
  gStyle->SetErrorX(0);
  std::string title="#deltam Distribution for Sum of Backgrounds in D_{s}^{*#pm} #rightarrow D_{s}^{#pm} e^{+} e^{-}";
  TCanvas *c_Stacked = new TCanvas("c_Stacked");
  THStack *s_DeltaM_sideband=new THStack("s_DeltaM_sideband", "#deltam Sidebands");
  s_DeltaM_sideband->Add(h_DeltaM_continu, "hist");
  s_DeltaM_sideband->Add(h_DeltaM_generic, "hist");
  s_DeltaM_sideband->Add(h_DeltaM_conver, "hist");
  s_DeltaM_sideband->Add(h_DeltaM_signal, "hist");
  s_DeltaM_sideband->SetMaximum(30.0/1.1);
  s_DeltaM_sideband->Draw();
  // s_DeltaM_sideband->SetTitle(title.c_str());
  s_DeltaM_sideband->SetTitle();
  s_DeltaM_sideband->GetXaxis()->SetTitle("#deltaM (GeV)");
  s_DeltaM_sideband->GetXaxis()->CenterTitle();
  s_DeltaM_sideband->GetXaxis()->SetRangeUser(lowerSide_xmin, upperSide_xmax);
  s_DeltaM_sideband->GetYaxis()->SetTitle("Number of Events / 5 MeV");
  s_DeltaM_sideband->GetYaxis()->CenterTitle();
  h_DeltaM_physics->SetMarkerStyle(20);
  h_DeltaM_physics->SetMarkerSize(0.9);
  h_DeltaM_physics->Draw("SAMEE1");
  f_physics->Draw("SAME");
  //line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
  line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
  //line=new TLine(xmin, ymax, xmin, 0); line->Draw();
  //line=new TLine(xmax, ymax, xmax, 0); line->Draw();
  line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
  //line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
  TLegend *legendDeltaM=new TLegend(.6, .9, .9, .6);
  legendDeltaM->AddEntry(h_DeltaM_signal, "Signal Simulation", "f");
  legendDeltaM->AddEntry(h_DeltaM_continu, "non-c#bar{c} Simulation", "f");
  legendDeltaM->AddEntry(h_DeltaM_generic, "c#bar{c} w/o Conversions Simulation", "f");
  legendDeltaM->AddEntry(h_DeltaM_conver_Dsp, "Conversions Simulation", "f");
  legendDeltaM->AddEntry(h_DeltaM_physics, "Data");
  legendDeltaM->SetFillColor(kWhite);
  legendDeltaM->Draw();
}
