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

double plotScale=500;

double mbcCut_center=2.112;
double mbcCut_range=0.006;
double sideband_dist=0.5;

double lowerSide_xmin=2.06;
double lowerSide_xmax=mbcCut_center-(1+2*sideband_dist)*mbcCut_range;
double xmin=mbcCut_center-mbcCut_range;
double xmax=mbcCut_center+mbcCut_range;
double upperSide_xmin=mbcCut_center+(1+2*sideband_dist)*mbcCut_range;
double upperSide_xmax=2.155;

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
  //return (par[0]+par[1]*T_1(x[0])+par[2]*T_2(x[0]))*pow(2.155-x[0], 0.5);
  return (par[0]+par[1]*x[0])*pow(2.155-x[0], 0.5);
}

Double_t physics(Double_t *x, Double_t *par)
{
  if (x[0]>lowerSide_xmax && x[0]<upperSide_xmin) TF1::RejectPoint();
  //return (par[0]+par[1]*T_1(x[0])+par[2]*T_2(x[0]))*pow(2.155-x[0], 0.5);
  return (par[0]+par[1]*x[0])*pow(2.155-x[0], 0.5);
}

Double_t physics_scaledMC(Double_t *x, Double_t *par)
{
  if (x[0]>lowerSide_xmax && x[0]<upperSide_xmin) TF1::RejectPoint();
  return (par[0]+par[1]*x[0])*pow(2.155-x[0], 0.5)*par[2];
}

int SidebandFitter_MBC()
{
  
  //TFile *file=new TFile("Sum_Sideband_MBC_002.root");
  TFile *file=new TFile("Sum_Unblinded_MBC_002.root");
  
  TH1D *h_MBC_signal = (TH1D*)gDirectory->Get("h_MBC_signal");         
  TH1D *h_MBC_conver_Dsp = (TH1D*)gDirectory->Get("h_MBC_conver_Dsp"); 
  TH1D *h_MBC_conver_Dsm = (TH1D*)gDirectory->Get("h_MBC_conver_Dsm"); 
  TH1D *h_MBC_generic = (TH1D*)gDirectory->Get("h_MBC_generic");       
  TH1D *h_MBC_continu = (TH1D*)gDirectory->Get("h_MBC_continu");       
  TH1D *h_MBC_physics = (TH1D*)gDirectory->Get("h_MBC_physics"); h_MBC_physics->SetStats(false); h_MBC_physics->SetLineColor(kBlack);
  
  TH1D *h_MBC_conver = (TH1D*)h_MBC_conver_Dsp->Clone();
  h_MBC_conver->Add(h_MBC_conver_Dsm);
  
  TH1D *h_Sum_background = (TH1D*)h_MBC_conver_Dsp->Clone();
  h_Sum_background->SetLineColor(kMagenta);
  h_Sum_background->Add(h_MBC_conver_Dsm);
  h_Sum_background->Add(h_MBC_generic);
  h_Sum_background->Add(h_MBC_continu);
  
  h_MBC_conver_Dsp->SetFillColor(kRed);  h_MBC_conver_Dsp->SetFillStyle(3244); h_MBC_conver_Dsp->SetStats(false);
  h_MBC_conver_Dsm->SetFillColor(kRed); h_MBC_conver_Dsm->SetFillStyle(3244); h_MBC_conver_Dsm->SetStats(false);
  h_MBC_conver->SetFillColor(kRed); h_MBC_conver->SetFillStyle(3244); h_MBC_conver->SetStats(false);
  h_MBC_generic->SetFillColor(kGreen); h_MBC_generic->SetFillStyle(3245); h_MBC_generic->SetStats(false);
  h_MBC_continu->SetFillColor(kBlue); h_MBC_continu->SetFillStyle(3003); h_MBC_continu->SetStats(false);
  h_MBC_signal->SetFillColor(kCyan); h_MBC_signal->SetFillStyle(3016); h_MBC_signal->SetStats(false);
  
  /*
  h_MBC_continu->SetFillColor(7);
  h_MBC_generic->SetFillColor(8);
  h_MBC_conver->SetFillColor(5);
  h_MBC_signal->SetFillColor(0);
  */
  double ymax=(h_MBC_signal->GetMaximum())*0.95;
  TLine *line;
  
  gROOT->SetStyle("Plain");
  
  TCanvas *c_FitBackground = new TCanvas("c_FitBackground");
  h_Sum_background->Draw();
  TF1 *f_Sum_background=new TF1("f_Sum_background", sum_background, lowerSide_xmin, upperSide_xmax, 2);
  f_Sum_background->SetLineColor(kBlack);
  f_Sum_background->SetLineWidth(2);
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
  
  // TF1 *f_physics=new TF1("f_physics", physics, lowerSide_xmin, upperSide_xmax, 2);
  // Scaled version of MC fit follows
  TF1 *f_physics=new TF1("f_physics", physics_scaledMC, lowerSide_xmin, upperSide_xmax, 3);
  f_physics->FixParameter(0, f_Sum_background->GetParameter(0));
  f_physics->FixParameter(1, f_Sum_background->GetParameter(1));
  
  f_physics->SetLineWidth(2);
  f_physics->SetLineColor(kBlack);
  h_MBC_physics->Fit(f_physics, "RLLEFM");
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
  std::string title="m_{BC} Distribution for Sum of Backgrounds in D_{s}^{*#pm} #rightarrow D_{s}^{#pm} e^{+} e^{-}";
  TCanvas *c_Stacked = new TCanvas("c_Stacked");
  THStack *s_MBC_sideband=new THStack("s_MBC_sideband", "\m_{BC} Sidebands");
  s_MBC_sideband->Add(h_MBC_continu, "hist");
  s_MBC_sideband->Add(h_MBC_generic, "hist");
  s_MBC_sideband->Add(h_MBC_conver, "hist");
  s_MBC_sideband->Add(h_MBC_signal, "hist");
  s_MBC_sideband->Draw();
  h_MBC_physics->SetMarkerStyle(20);
  h_MBC_physics->SetMarkerSize(0.9);
  h_MBC_physics->Draw("SAMEE1");
  f_physics->Draw("SAME");
  // s_MBC_sideband->SetTitle(title.c_str());
  s_MBC_sideband->SetTitle();
  s_MBC_sideband->GetXaxis()->SetTitle("M_{BC} (GeV)");
  s_MBC_sideband->GetXaxis()->CenterTitle();
  s_MBC_sideband->GetXaxis()->SetRangeUser(lowerSide_xmin, upperSide_xmax);
  s_MBC_sideband->GetYaxis()->SetTitle("Number of Events / 2 MeV");
  s_MBC_sideband->GetYaxis()->CenterTitle();
  //line=new TLine(lowerSide_xmin, ymax, lowerSide_xmin, 0); line->Draw();
  line=new TLine(lowerSide_xmax, ymax, lowerSide_xmax, 0); line->Draw();
  // line=new TLine(xmin, ymax, xmin, 0); line->Draw();
  // line=new TLine(xmax, ymax, xmax, 0); line->Draw();
  line=new TLine(upperSide_xmin, ymax, upperSide_xmin, 0); line->Draw();
  //line=new TLine(upperSide_xmax, ymax, upperSide_xmax, 0); line->Draw();
  TLegend *legendMBC=new TLegend(.6, .9, .9, .6);
  legendMBC->AddEntry(h_MBC_signal, "Signal Simulation",  "f");
  legendMBC->AddEntry(h_MBC_continu, "non-c#bar{c} Simulation", "f");
  legendMBC->AddEntry(h_MBC_generic, "c#bar{c} w/o Conversions Simulation", "f");
  legendMBC->AddEntry(h_MBC_conver_Dsp, "Conversions Simulation", "f");
  legendMBC->AddEntry(h_MBC_physics, "Data");
  legendMBC->SetFillColor(kWhite);
  legendMBC->Draw();
  s_MBC_sideband->SetMaximum(25.0/1.1);
}
