#include "TCanvas.h"
#include "TH1D.h"
#include "TROOT.h"
#include "THStack.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"x
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include <iostream>

double constant_;

std::string decay="KKpi";
// options: KKpi, KsK, pieta, pietaprime, KKpipi0, pipipi, KsKmpipi, pipi0eta, pietaprimerho

std::string texDsstgamma="D_{s}^{*} #rightarrow D_{s} #gamma";
double mbcCut_center, mbcCut_range;
double mbc_x1, mbc_x2, mbc_y1, mbc_y2;
double xmin, xmax;
double par_min[32], par_max[32];

double luminosity=586; // /pb
double luminosity_error=6;
double prodCrossSection_DsDss=948;
double prodCrossSection_DsDss_error=36;
double branchingFr_mode;
double branchingFr_mode_error=0;
double branchingFr_mode_generic;

double nConversionSample_Dsp=99880;
double nConversionSample_Dsm=99880;

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

double T_4(double x)
{
  return 8*x*x*x*x-8*x*x+1;
}

double T_5(double x)
{
  return 16*x*x*x*x*x-20*x*x*x+5*x;
}

Double_t argus(Double_t *x, Double_t *par) // 6 Parameters
{
  Double_t result=0;
  // if (par[0]>x[0]) result=(par[2]+par[3]*T_1(x[0])+par[4]*T_2(x[0])+par[5]*T_3(x[0]))*pow(par[0]-x[0], par[1]);
  if (par[0]>x[0]) result=(par[2]+par[3]*x[0]+par[4]*pow(x[0],2)+par[5]*pow(x[0], 3))*pow(par[0]-x[0], par[1]);
  return result;
}

Double_t crystalBall(Double_t *x, Double_t *par) // with a soft Gaussian - 8 Parameters
{
  // Double shouldered 
  Double_t std=(x[0]-par[0])/par[1];
  Double_t A=pow(par[3]/par[2], par[3])*exp(-0.5*pow(par[2], 2));
  Double_t B=par[3]/par[2]-par[2];
  Double_t result=0;
  if (std<=par[2]) // Gaussian Region
  {
    result=exp(-0.5*pow(std, 2));
  }
  else // Power Law Region
  {
    result=A/pow(B+std, par[3]);
  }
  
  // Soft Gaussian on the right
  result+=par[7]*exp(-0.5*pow((x[0]-par[5])/par[6], 2));
  
  // Scale factor
  result=result*par[4];
  
  return result;
}

Double_t converFit(Double_t *x, Double_t *par) // 8+6=14 Parameters
{
  Double_t result;
  
  result=crystalBall(x, par);
  result=result+argus(x, par+8);
  
  return result;
}

Double_t wrongDsMatchedPhoton(Double_t *x, Double_t *par) // 11 Parameters
{
  Double_t result;
  
  result=crystalBall(x, par); // 8 parameters (Right side CB with soft Gaussian on tail)
  result+=par[8]*exp(-0.5*pow((x[0]-par[9])/par[10], 2)); // +3
  
  return result;
}

Double_t wrongDsUnmatchedPhoton(Double_t *x, Double_t *par) // 6 Parameters
{
  Double_t result=0;
  
  Double_t std=(x[0]-par[0])/par[1];
  Double_t A=pow(par[3]/par[2], par[3])*exp(-0.5*pow(par[2], 2));
  Double_t B=par[3]/par[2]-par[2];
  Double_t a=(1+pow(par[4]/par[1], 2))*exp(-0.5*pow(par[4], 2));
  Double_t b=-(par[4]/pow(par[1], 2))*exp(-0.5*pow(par[4], 2));
  if (std<=par[2] && std>=par[4]) // Gaussian Region
  {
    result=exp(-0.5*pow(std, 2));
  }
  else if (std>par[2]) // Power Law Region
  {
    result=A/pow(B+std, par[3]);
  }
  else if (std<par[4])
  {
    result=a+(b*std);
  }
  result=result*par[5];
  
  // std::cout<<"std = "<<std<<", par[4] = "<<par[4]<<", a = "<<a<<", b = "<<b<<std::endl;
  
  return result;
}

Double_t crystalBall_double(Double_t *x, Double_t *par) // 7 Parameters
{
  // Double shouldered 
  Double_t std=(x[0]-par[0])/par[1];
  Double_t A=pow(par[3]/par[2], par[3])*exp(-0.5*pow(par[2], 2));
  Double_t B=par[3]/par[2]-par[2];
  Double_t C=pow(par[5]/fabs(par[4]), par[5])*exp(-0.5*pow(par[4], 2));
  Double_t D=par[5]/fabs(par[4])-fabs(par[4]);
  Double_t result=0;
  if (std>=par[4] && std<=par[2]) // Gaussian Region
  {
    result=exp(-0.5*pow(std, 2));
  }
  else if (std>par[2]) // Power Law Region
  {
    result=A/pow(B+std, par[3]);
  }
  else if (std<par[4]) // Power Law Region
  {
    result=C/pow(D-std, par[5]);
  }
  result=result*par[6];
  
  return result;
}

Double_t wrongConverFit(Double_t *x, Double_t *par) // 14 Parameters
{
  Double_t result;
  
  result=crystalBall_double(x, par); // 7 parameters
  
  // Now add the second Gaussian
  result+=par[7]*exp(-0.5*pow((x[0]-par[8])/par[9], 2));
  
  // Now add the background function
  Double_t width=(x[0]-par[10])/par[11];
  if (x[0]>par[10]) result=result+exp(-pow(width, 2))*pow(width, par[12])*par[13];
  
  return result;
}

/*
Double_t dataFit_wrongConver(Double_t *x, Double_t *par) // 6+14+1=21 Parameters
{
  Double_t result;
  result=argus(x, par); // 6 Parameters
  result+=wrongConverFit(x, par+6)*par[20];
  
  return result;
}
*/

Double_t dataFit_wrongDsUnmatchedPhoton(Double_t *x, Double_t *par) // 12 parameters
{
  Double_t result;
  result=argus(x, par); // 6 parameters
  result+=wrongDsUnmatchedPhoton(x, par+6); // 6 parameters
  return result;
}

Double_t dataFit_wrongDsMatchedPhoton(Double_t *x, Double_t *par) // 24 parameters
{
  Double_t result;
  result=dataFit_wrongDsUnmatchedPhoton(x, par); // 12 parameters
  result+=wrongDsMatchedPhoton(x, par+12)*par[23]; // 12 parameters
  return result;
}

Double_t dataFit(Double_t *x, Double_t *par) // 32 parameters
{
  Double_t result;
  result=dataFit_wrongDsMatchedPhoton(x, par); // 24 parameters
  result+=crystalBall(x, par+24);              // 8 parameters
  return result;
}

Double_t dataFit1(Double_t *x, Double_t *par) // 32 parameters
{
  Double_t result;
  result=dataFit_wrongDsMatchedPhoton(x, par); // 24 parameters
  
  Double_t result1=0;
  Double_t std=(x[0]-par[24])/par[25];
  Double_t A=pow(par[27]/par[26], par[27])*exp(-0.5*pow(par[26], 2));
  Double_t B=par[27]/par[26]-par[26];
  if (std<=par[26]) // Gaussian Region
  {
    result1=exp(-0.5*pow(std, 2));
  }
  else // Power Law Region
  {
    result1=A/pow(B+std, par[27]);
  }
  
  // Soft Gaussian on the right
  result1+=par[31]*exp(-0.5*pow((x[0]-par[29])/par[30], 2));
  
  // Scale factor
  result1=result1*(par[23]/constant_);
  //std::cout<<" constant_ = "<<constant_;
  
  return result+result1;
}

void setValues()
{
  if (decay=="KKpi")
  {
    mbcCut_center=2.112; mbcCut_range=0.004;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
    xmin=2.08; xmax=2.150;
    
    branchingFr_mode=0.055;
    branchingFr_mode_error=0.0028;
    branchingFr_mode_generic=0.0537;
    
    nConversionSample_Dsp=99880;
    nConversionSample_Dsm=99880;
    
    par_xmin={2.140, 0.1, 0., 0., 0., 0., 2.125, 0.001, 0.1, 1.0, -0.1, 0., 2.111, 0.001, 0.1, 1.0, 0., 2.125, 0.005, 0., 0., 2.110, 0.001, 0., 2.111, 0.001, 1.1, 2.0, 0., 2.11, 0.005, 0.};
    par_xmax={2.155, 1.0, 0., 0., 0., 0., 2.135, 0.01,  2.0, 5.0, 0.0,  1., 2.115, 0.01,  2.0, 5.0, 0., 2.135, 0.01,  0., 0., 2.112, 0.1,   0., 2.113, 0.1,   1.6, 4.0, 0., 2.13, 0.02,  0.};
  }
  else if (decay=="KsK")
  {
    mbcCut_center=2.112; mbcCut_range=0.007;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.35; mbc_y2=0.6;
    
    xmin=2.08; xmax=2.150;
    
    branchingFr_mode=0.0149;
    branchingFr_mode_error=0.0009;
    branchingFr_mode_generic=0.0293*0.5;
  }
  else if (decay=="pieta")
  {
    mbcCut_center=2.112; mbcCut_range=0.008;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
    xmin=2.08; xmax=2.150;
    
    branchingFr_mode=0.0158*0.3931;
    branchingFr_mode_error=branchingFr_mode*pow(pow(0.0021/0.0158, 2)+pow(0.0020/0.3931, 2), 0.5);
    branchingFr_mode_generic=0.0154*0.39466;
  }
  else if (decay=="pietaprime")
  {
    mbcCut_center=2.112; mbcCut_range=0.011;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
    xmin=2.08; xmax=2.150;
    
    branchingFr_mode=0.038*0.446*0.3931;
    branchingFr_mode_error=branchingFr_mode*pow(pow(0.004/0.038, 2)+pow(0.0014/0.446, 2)+pow(0.0020/0.3931, 2), 0.5);
    branchingFr_mode_generic=0.0367*0.4370*0.39466;
  }
  else if (decay=="KKpipi0")
  {
    mbcCut_center=2.112; mbcCut_range=0.004;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
    xmin=2.08; xmax=2.150;
    
    branchingFr_mode=0.056;
    branchingFr_mode_error=0.005;
    branchingFr_mode_generic=0;
  }
  else if (decay=="pipipi")
  {
    mbcCut_center=2.112; mbcCut_range=0.004;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
    xmin=2.08; xmax=2.150;
    
    branchingFr_mode=0.0111;
    branchingFr_mode_error=0.0008;
  }
  else if (decay=="KsKmpipi")
  {
    mbcCut_center=2.112; mbcCut_range=0.005;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
    xmin=2.08; xmax=2.150;
    
    branchingFr_mode=0.0164;
    branchingFr_mode_error=0.0012;
    branchingFr_mode_generic=0.0741*0.6657*0.660*0.5;
  }
  else if (decay=="pipi0eta")
  {
    mbcCut_center=2.112; mbcCut_range=0.004;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
    xmin=2.08; xmax=2.150;
    
    branchingFr_mode=0.130*1.*0.3931;
    branchingFr_mode_error=branchingFr_mode*pow(pow(0.022/0.130, 2)+pow(0.0020/0.3931, 2), 0.5);
    branchingFr_mode_generic=0.0758*1.*0.3931;
  }
  else if (decay=="pietaprimerho")
  {
    mbcCut_center=2.112; mbcCut_range=0.004;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
   xmin=2.08; xmax=2.150;
    
    branchingFr_mode=0.038*0.294;
    branchingFr_mode_error=branchingFr_mode*pow(pow(0.004/0.038, 2)+pow(0.009/0.294, 2), 0.5);
    branchingFr_mode_generic=0.0367*0.302;
  }
  
}

void DsGammaAnalysis_MBC_Fitter()
{

  setValues();
  
  std::string filename;
  filename=decay;
  filename+="_DsGamma_MBC.root";
  TFile *file=new TFile(filename.c_str());
  
  TH1D *h_MBC_conver = (TH1D*)gDirectory->Get("h_MBC_conver");
  TH1D *h_MBC_generic_veto = (TH1D*)gDirectory->Get("h_MBC_generic_veto");
  TH1D *h_MBC_generic_right = (TH1D*)gDirectory->Get("h_MBC_generic");
  TH1D *h_MBC_generic_wrong = (TH1D*)gDirectory->Get("h_MBC_generic_wrong");
  TH1D *h_MBC_generic = new TH1D(*h_MBC_generic_right);
  TH1D *h_MBC_continu = (TH1D*)gDirectory->Get("h_MBC_continu");
  h_MBC_generic->Add(h_MBC_generic_wrong);
  h_MBC_generic->Add(h_MBC_generic_veto);
  h_MBC_generic->Add(h_MBC_continu);
  TH1D *h_MBC_physics = (TH1D*)gDirectory->Get("h_MBC_physics");
  TH1D *h_MBC_wrongConver = (TH1D*)gDirectory->Get("h_MBC_wrongConver");
  //filename=decay;
  //filename+="_DsGamma_MBC_MatchedOnly_1.root";
  //TFile *file1=new TFile(filename.c_str());
  TH1D *h_MBC_conver_matched = (TH1D*)gDirectory->Get("h_MBC_conver_matched");
  TH1D *h_MBC_conver_unmatched = (TH1D*)gDirectory->Get("h_MBC_conver_unmatched");
  TH1D *h_MBC_wrongConver_matched = (TH1D*)gDirectory->Get("h_MBC_wrongConver_matched");
  TH1D *h_MBC_wrongConver_unmatched = (TH1D*)gDirectory->Get("h_MBC_wrongConver_unmatched");
  TH1D *h_MBC_wrongConver_strictUnmatched = (TH1D*)gDirectory->Get("h_MBC_wrongConver_strictUnmatched");
  
  /*
  h_MBC_generic_veto->SetFillColor(kGreen);
  h_MBC_generic_right->SetFillColor(kGreen);
  h_MBC_generic_wrong->SetLineColor(kRed);
  h_MBC_continu->SetFillColor(kBlue);
  */
  
  double ymax;
  TLine *line;
  
  gROOT->SetStyle("Plain");
  
  // First fit the MC matched Ds gamma signal sample
  // with a single shouldered Crystal Ball function
  TCanvas *c_MBC_conver_matched=new TCanvas("c_MBC_conver_matched", "c_MBC_conver_matched");
  h_MBC_conver_matched->Draw();
  TF1 *f_converFit_matched=new TF1("f_converFit_matched", crystalBall, xmin, xmax, 8);
  f_converFit_matched->SetLineWidth(0);
  f_converFit_matched->SetParLimits(0, 2.111, 2.113);
  f_converFit_matched->SetParLimits(1, 0.001, 0.1);
  f_converFit_matched->SetParLimits(2, 1.1, 1.6);
  f_converFit_matched->SetParLimits(3, 2.0, 4.0);
  f_converFit_matched->SetParLimits(5, 2.11, 2.13); // Gaussian center
  f_converFit_matched->SetParLimits(6, 0.005, 0.02);
  h_MBC_conver_matched->Fit(f_converFit_matched, "REFM");
  
  // Now fit the right-sign Ds with photon unmatched
  TCanvas *c_MBC_conver_unmatched=new TCanvas("c_MBC_conver_unmatched");
  h_MBC_conver_unmatched->Draw();
  TF1 *f_converFit_unmatched=new TF1("f_converFit_unmatched", argus, xmin, xmax, 6);
  f_converFit_unmatched->SetLineWidth(0);
  f_converFit_unmatched->SetParLimits(0, 2.140, 2.155);
  f_converFit_unmatched->SetParLimits(1, 0.1, 1.0);
  h_MBC_conver_unmatched->Fit(f_converFit_unmatched, "REFM");
  
  // Fit the conver sample
  std::string title;
  title="m_{BC} Distribution in Signal Sample of ";
  title+=texDsstgamma;
  title+=", D_{s} #rightarrow ";
  title+=decay;
  TCanvas *c_MBC_conver=new TCanvas("c_MBC_conver", "c_MBC_conver");
  h_MBC_conver->SetTitle(title.c_str());
  h_MBC_conver->GetYaxis()->SetTitle("Efficiency / MeV");
  h_MBC_conver->GetYaxis()->CenterTitle();
  h_MBC_conver->GetYaxis()->SetTitleOffset(1.3);
  h_MBC_conver->Draw();
  TF1 *f_converFit=new TF1("f_converFit", converFit, xmin, xmax, 14);
  // define another CB function.
  f_converFit->FixParameter(0, f_converFit_matched->GetParameter(0));
  f_converFit->FixParameter(1, f_converFit_matched->GetParameter(1));
  f_converFit->FixParameter(2, f_converFit_matched->GetParameter(2));
  f_converFit->FixParameter(3, f_converFit_matched->GetParameter(3));
  f_converFit->FixParameter(5, f_converFit_matched->GetParameter(5));
  f_converFit->FixParameter(6, f_converFit_matched->GetParameter(6));
  f_converFit->FixParameter(7, f_converFit_matched->GetParameter(7));
  f_converFit->SetParLimits(8, 2.140, 2.155);
  f_converFit->SetParLimits(9, 0.1, 1.0);
  f_converFit->SetLineWidth(0);
  f_converFit->SetNpx(1000);
  h_MBC_conver->Fit(f_converFit, "REFM");
  TF1 *f_converFit_bg=new TF1("f_converFit_bg", argus, xmin, xmax, 6);
  Double_t converPar[14];
  f_converFit->GetParameters(converPar);
  f_converFit_bg->SetParameters(converPar+8);
  f_converFit_bg->SetLineWidth(0);
  f_converFit_bg->Draw("SAME");
  ymax=(h_MBC_conver->GetMaximum())*0.95;
  line=new TLine(xmin, 0, xmin, ymax); line->Draw();
  line=new TLine(xmax, 0, xmax, ymax); line->Draw();
  std::string filename_eps=decay;
  std::string filename_png=decay;
  filename_eps+="_DsGammaEff_MBC.eps";
  filename_png+="_DsGammaEff_MBC.png";
  c_MBC_conver->Print(filename_eps.c_str());
  c_MBC_conver->Print(filename_png.c_str());
  double nRegion_MC=(f_converFit->Integral(xmin, xmax))*(120.0/0.12);
  double nRegion_MC_errorFr=(f_converFit->GetParError(4))/(f_converFit->GetParameter(4));
  double nCombo_MC=(f_converFit_bg->Integral(xmin, xmax))*(120.0/0.12);
  double nSignal_MC=nRegion_MC-nCombo_MC;
  double nSignal_MC_error=nSignal_MC*nRegion_MC_errorFr;
  std::cout<<"- MC -"<<std::endl;
  std::cout<<"Number of signal MC events in Region = "<<nRegion_MC<<std::endl;
  std::cout<<"Number of signal MC events as Combinatorics in the Region = "<<nCombo_MC<<std::endl;;
  std::cout<<"Number of signal MC events identified as signal = "<<nSignal_MC<<"+-"<<nSignal_MC_error<<std::endl;
  std::cout<<"---"<<std::endl;
  
  // Wrong sign legwork follows
  
  // Matched wrong sign Ds + matched photon
  TCanvas *c_wrongConver_matched = new TCanvas("c_wrongConver_matched");
  h_MBC_wrongConver_matched->Draw();
  TF1 *f_MBC_wrongConver_matched=new TF1("f_MBC_wrongConver_matched", wrongDsMatchedPhoton, xmin, xmax, 13);
  f_MBC_wrongConver_matched->SetLineWidth(0);
  f_MBC_wrongConver_matched->SetParLimits(0, 2.111, 2.115);
  f_MBC_wrongConver_matched->SetParLimits(1, 0.001, 0.01);
  f_MBC_wrongConver_matched->SetParLimits(2, 0.1, 2.0);
  f_MBC_wrongConver_matched->SetParLimits(3, 1.0, 5.0);
  f_MBC_wrongConver_matched->SetParLimits(5, 2.125, 2.135); // Gaussian center
  f_MBC_wrongConver_matched->SetParLimits(6, 0.005, 0.01);
  f_MBC_wrongConver_matched->SetParLimits(9, 2.110, 2.112);
  f_MBC_wrongConver_matched->SetParLimits(10, 0.001, 0.1);
  h_MBC_wrongConver_matched->Fit(f_MBC_wrongConver_matched, "REFM");
  
  // Matched wrong sign Ds + unmatched photon
  TCanvas *c_wrongConver_unmatched = new TCanvas("c_wrongConver_unmatched");
  h_MBC_wrongConver_unmatched->Draw();
  TF1 *f_MBC_wrongConver_unmatched=new TF1("f_MBC_wrongConver_unmatched", wrongDsUnmatchedPhoton, 2.08, 2.150, 6);
  f_MBC_wrongConver_unmatched->SetLineWidth(0);
  f_MBC_wrongConver_unmatched->SetParLimits(0, 2.125, 2.130);
  f_MBC_wrongConver_unmatched->SetParLimits(1, 0.001, 0.01);
  f_MBC_wrongConver_unmatched->SetParLimits(2, 0.1, 2.0);
  f_MBC_wrongConver_unmatched->SetParLimits(3, 1.0, 5.0);
  f_MBC_wrongConver_unmatched->SetParLimits(4, -0.1, 0.0);
  f_MBC_wrongConver_unmatched->SetParLimits(5, 0., 1.);
  h_MBC_wrongConver_unmatched->Fit(f_MBC_wrongConver_unmatched, "REFM");
  
  /*
  // Second, fit the wrongConver sample
  title="m_{BC} Distribution in Wrong-Sign D_{s} #rightarrow ";
  title+=decay;
  TCanvas *c_wrongConver = new TCanvas("c_wrongConver", "c_wrongConver");
  c_wrongConver->Divide(1,2);
  c_wrongConver->cd(1);
  h_MBC_wrongConver->SetTitle(title.c_str());
  h_MBC_wrongConver->GetYaxis()->SetTitle("# Events / MeV");
  h_MBC_wrongConver->GetYaxis()->CenterTitle();
  h_MBC_wrongConver->GetYaxis()->SetTitleOffset(1.2);
  h_MBC_wrongConver->Draw();
  TF1* f_wrongConverFit=new TF1("f_wrongConverFit", wrongConverFit, xmin, xmax, 14);
  f_wrongConverFit->SetParLimits(0, 2.112, 2.118); // 2.114 low
  f_wrongConverFit->SetParLimits(1, 0.001, 0.1);
  f_wrongConverFit->SetParLimits(2, 1.0, 5.0);
  f_wrongConverFit->SetParLimits(3, 0.1, 2.0);
  f_wrongConverFit->SetParLimits(4, -3.0, -0.1);
  f_wrongConverFit->SetParLimits(5, 0.5, 5.0); // 1.0 low
  f_wrongConverFit->SetParLimits(8, 2.127, 2.135);
  f_wrongConverFit->SetParLimits(9, 0.001, 0.01);
  f_wrongConverFit->SetParLimits(10, 2.05, 2.07);
  f_wrongConverFit->SetParLimits(11, 0.005, 0.19);
  f_wrongConverFit->SetParLimits(12, 0.5, 2.0);
  f_wrongConverFit->SetLineWidth(0);
  h_MBC_wrongConver->Fit(f_wrongConverFit, "R");
  c_wrongConver->cd(2);
  h_MBC_generic_wrong->SetLineColor(kCyan);
  h_MBC_wrongConver->DrawNormalized();
  h_MBC_generic_wrong->DrawNormalized("SAME");
  filename_eps=decay;
  filename_png=decay;
  filename_eps+="_wrongDsGammaEff_MBC.eps";
  filename_png+="_wrongDsGammaEff_MBC.png";
  c_wrongConver->Print(filename_eps.c_str());
  c_wrongConver->Print(filename_png.c_str());
  double wrongSign_eff=(f_wrongConverFit->Integral(xmin, xmax))*(120.0/0.12);
  std::cout<<"Wrong sign efficiency = "<<wrongSign_eff<<std::endl;
  */
  
  // Stack Backgrounds
  THStack *s_MBC_Background=new THStack("s_MBC_Background", "");
  s_MBC_Background->Add(h_MBC_continu, "hist");
  s_MBC_Background->Add(h_MBC_generic_veto, "hist");
  s_MBC_Background->Add(h_MBC_generic_wrong, "hist");
  
  // Now fit generic MC
  TCanvas *c_genericFit = new TCanvas("genericFit", "genericFit");
  title="m_{BC} Distribution in generic for ";
  title+=texDsstgamma;
  title+=", D_{s} #rightarrow ";
  title+=decay;
  h_MBC_generic->SetTitle(title.c_str());
  h_MBC_generic->GetXaxis()->SetTitle("m_{BC} (GeV)");
  h_MBC_generic->GetYaxis()->SetTitle("# Events / MeV");
  h_MBC_generic->GetYaxis()->CenterTitle();
  h_MBC_generic->GetYaxis()->SetTitleOffset(1.2);
  TF1 *f_genericFit=new TF1("f_genericFit", dataFit, xmin, xmax, 32);
  f_genericFit->SetParLimits(0, 2.140, 2.155); // right limit of Argus function
  f_genericFit->SetParLimits(1, 0.1, 1.0);
  // par 2 movable - argus
  // par 3 movable - argus
  // par 4 movable - argus
  // par 5 movable - argus
  f_genericFit->FixParameter(6, f_MBC_wrongConver_unmatched->GetParameter(0));
  f_genericFit->FixParameter(7, f_MBC_wrongConver_unmatched->GetParameter(1));
  f_genericFit->FixParameter(8, f_MBC_wrongConver_unmatched->GetParameter(2));
  f_genericFit->FixParameter(9, f_MBC_wrongConver_unmatched->GetParameter(3));
  f_genericFit->FixParameter(10, f_MBC_wrongConver_unmatched->GetParameter(4));
  // par 11 movable - scale of wrong sign Ds, unmatched photon
  f_genericFit->FixParameter(12, f_MBC_wrongConver_matched->GetParameter(0));
  f_genericFit->FixParameter(13, f_MBC_wrongConver_matched->GetParameter(1));
  f_genericFit->FixParameter(14, f_MBC_wrongConver_matched->GetParameter(2));
  f_genericFit->FixParameter(15, f_MBC_wrongConver_matched->GetParameter(3));
  f_genericFit->FixParameter(16, f_MBC_wrongConver_matched->GetParameter(4));
  f_genericFit->FixParameter(17, f_MBC_wrongConver_matched->GetParameter(5));
  f_genericFit->FixParameter(18, f_MBC_wrongConver_matched->GetParameter(6));
  f_genericFit->FixParameter(19, f_MBC_wrongConver_matched->GetParameter(7));
  f_genericFit->FixParameter(20, f_MBC_wrongConver_matched->GetParameter(8));
  f_genericFit->FixParameter(21, f_MBC_wrongConver_matched->GetParameter(9));
  f_genericFit->FixParameter(22, f_MBC_wrongConver_matched->GetParameter(10));
  // par 23 movable - scale of wrong sign Ds, matched photon
  f_genericFit->SetParLimits(24, 2.111, 2.113);
  f_genericFit->FixParameter(25, f_converFit->GetParameter(1));
  f_genericFit->FixParameter(26, f_converFit->GetParameter(2));
  f_genericFit->FixParameter(27, f_converFit->GetParameter(3));
  // par 28 movable - scale of signal shape
  f_genericFit->FixParameter(29, f_converFit->GetParameter(5));
  f_genericFit->FixParameter(30, f_converFit->GetParameter(6));
  f_genericFit->FixParameter(31, f_converFit->GetParameter(7));
  f_genericFit->SetLineWidth(0);
  f_genericFit->SetNpx(1000);
  h_MBC_generic->Fit(f_genericFit, "REFM");
  s_MBC_Background->SetMaximum(h_MBC_generic->GetMaximum());
  // s_MBC_Background->Draw();
  h_MBC_generic->Draw("SAME");
  Double_t genericFitPar[32];
  f_genericFit->GetParameters(genericFitPar);
  TF1 *f_genericFit_argus=new TF1("f_genericFit_argus", argus, xmin, xmax, 6);
  f_genericFit_argus->SetParameters(genericFitPar);
  f_genericFit_argus->SetLineWidth(0);
  f_genericFit_argus->Draw("SAME");
  TF1 *f_genericFit_wrongDsUnmatchedPhoton=new TF1("f_genericFit_wrongDsUnmatchedPhoton", dataFit_wrongDsUnmatchedPhoton, xmin, xmax, 12);
  f_genericFit_wrongDsUnmatchedPhoton->SetParameters(genericFitPar);
  f_genericFit_wrongDsUnmatchedPhoton->SetLineWidth(0);
  f_genericFit_wrongDsUnmatchedPhoton->Draw("SAME");
  TF1 *f_genericFit_wrongDsMatchedPhoton=new TF1("f_genericFit_wrongDsMatchedPhoton", dataFit_wrongDsMatchedPhoton, xmin, xmax, 24);
  f_genericFit_wrongDsMatchedPhoton->SetParameters(genericFitPar);
  f_genericFit_wrongDsMatchedPhoton->SetLineWidth(0);
  f_genericFit_wrongDsMatchedPhoton->Draw("SAME");
  ymax=(h_MBC_generic->GetMaximum())*0.95;
  line=new TLine(xmin, 0, xmin, ymax); line->Draw();
  line=new TLine(xmax, 0, xmax, ymax); line->Draw();
  /*
  filename_eps=decay;
  filename_png=decay;
  filename_eps+="_DsGammaGeneric_MBC.eps";
  filename_png+="_DsGammaGeneric_MBC.png";
  c_genericFit->SaveAs(filename_eps.c_str());
  c_genericFit->SaveAs(filename_png.c_str()); 
  */
  constant_ = genericFitPar[23]/genericFitPar[28];
  double nRegion_generic=(f_genericFit->Integral(xmin, xmax))*(120.0/0.12);
  double nBackground_generic=(f_genericFit_wrongDsMatchedPhoton->Integral(xmin, xmax))*(120.0/0.12);
  double nArgusBackground_generic=(f_genericFit_argus->Integral(xmin, xmax))*(120.0/0.12);
  double nSignal_generic=nRegion_generic-nBackground_generic;
  double scaleRegion_errorFr_generic=(f_genericFit->GetParError(28)/f_genericFit->GetParameter(28));
  double nSignal_generic_error=nSignal_generic*scaleRegion_errorFr_generic;
  double branchingFraction_generic=nSignal_generic/(luminosity*prodCrossSection_DsDss*branchingFr_mode_generic*nSignal_MC);
  double branchingFraction_error_generic=branchingFraction_generic*pow(
                                                       pow(1./1052., 2)
                                                      +pow(nSignal_MC_error/nSignal_MC, 2) 
                                                      +pow(scaleRegion_errorFr_generic, 2), 0.5);
  double n_wrongDsUnmatchedPhoton=(f_genericFit_wrongDsUnmatchedPhoton->Integral(xmin, xmax))*(120./0.12);
  std::cout<<"ratio = "<<(nBackground_generic-n_wrongDsUnmatchedPhoton)/(nRegion_generic-nBackground_generic)<<std::endl;
  std::cout<<"constant_ = "<<constant_<<std::endl;
  /*
  double nWrongSignBackground=nBackground_generic-nArgusBackground_generic;
  double nWrongSignScale_errorFr=(f_genericFit_wrongDsMatchedPhoton->GetParError(18)/f_genericFit_wrongDsMatchedPhoton->GetParameter(18));
  double nWrongSignBackground_error=nWrongSignBackground*nWrongSignScale_errorFr;
  double wrongSignProd=nWrongSignBackground/wrongSign_eff;
  double branchingFraction_wrongSign=nWrongSignBackground/(luminosity*prodCrossSection_DsDss*wrongSign_eff);
  double branchingFraction_wrongSign_error=branchingFraction_wrongSign*pow(
                                                        pow(1./1052., 2)
                                                       +pow(nWrongSignBackground_error/nWrongSignBackground, 2), 0.5);
  */
  
  // Now fit the data
  TCanvas *c_dataFit = new TCanvas("dataFit", "dataFit");
  title="m_{BC} Distribution in Data for ";
  title+=texDsstgamma;
  title+=", D_{s} #rightarrow ";
  title+=decay;
  h_MBC_physics->SetTitle(title.c_str());
  h_MBC_physics->GetXaxis()->SetTitle("m_{BC} (GeV)");
  h_MBC_physics->GetYaxis()->SetTitle("# Events / MeV");
  h_MBC_physics->GetYaxis()->CenterTitle();
  h_MBC_physics->GetYaxis()->SetTitleOffset(1.2);
  TF1 *f_dataFit=new TF1("f_dataFit", dataFit1, xmin, xmax, 32);
  f_dataFit->SetParLimits(0, 2.140, 2.155); // right limit of Argus function
  f_dataFit->SetParLimits(1, 0.1, 1.0);
  // par 1 movable - argus
  // par 2 movable - argus
  // par 3 movable - argus
  // par 4 movable - argus
  // par 5 movable - argus
  f_dataFit->FixParameter(6, f_MBC_wrongConver_unmatched->GetParameter(0));
  f_dataFit->FixParameter(7, f_MBC_wrongConver_unmatched->GetParameter(1));
  f_dataFit->FixParameter(8, f_MBC_wrongConver_unmatched->GetParameter(2));
  f_dataFit->FixParameter(9, f_MBC_wrongConver_unmatched->GetParameter(3));
  f_dataFit->FixParameter(10, f_MBC_wrongConver_unmatched->GetParameter(4));
  // par 11 movable - scale of wrong sign Ds, unmatched photon
  f_dataFit->FixParameter(12, f_MBC_wrongConver_matched->GetParameter(0));
  f_dataFit->FixParameter(13, f_MBC_wrongConver_matched->GetParameter(1));
  f_dataFit->FixParameter(14, f_MBC_wrongConver_matched->GetParameter(2));
  f_dataFit->FixParameter(15, f_MBC_wrongConver_matched->GetParameter(3));
  f_dataFit->FixParameter(16, f_MBC_wrongConver_matched->GetParameter(4));
  f_dataFit->FixParameter(17, f_MBC_wrongConver_matched->GetParameter(5));
  f_dataFit->FixParameter(18, f_MBC_wrongConver_matched->GetParameter(6));
  f_dataFit->FixParameter(19, f_MBC_wrongConver_matched->GetParameter(7));
  f_dataFit->FixParameter(20, f_MBC_wrongConver_matched->GetParameter(8));
  f_dataFit->FixParameter(21, f_MBC_wrongConver_matched->GetParameter(9));
  f_dataFit->FixParameter(22, f_MBC_wrongConver_matched->GetParameter(10));
  // par 23 movable - scale of wrong sign Ds, matched photon
  f_dataFit->SetParLimits(24, 2.111, 2.113);
  f_dataFit->FixParameter(25, f_converFit->GetParameter(1));
  f_dataFit->FixParameter(26, f_converFit->GetParameter(2));
  f_dataFit->FixParameter(27, f_converFit->GetParameter(3));
  // par 28 movable - scale of signal shape
  f_dataFit->FixParameter(29, f_converFit->GetParameter(5));
  f_dataFit->FixParameter(30, f_converFit->GetParameter(6));
  f_dataFit->FixParameter(31, f_converFit->GetParameter(7));
  f_dataFit->SetLineWidth(0);
  h_MBC_physics->Fit(f_dataFit, "R");
  Double_t dataFitPar[32];
  f_dataFit->GetParameters(dataFitPar);
  TF1 *f_dataFit_argus=new TF1("f_dataFit_argus", argus, xmin, xmax, 6);
  f_dataFit_argus->SetParameters(dataFitPar);
  f_dataFit_argus->SetLineWidth(0);
  f_dataFit_argus->Draw("SAME");
  TF1 *f_dataFit_wrongDsUnmatchedPhoton=new TF1("f_dataFit_wrongDsUnmatchedPhoton", dataFit_wrongDsUnmatchedPhoton, xmin, xmax, 12);
  f_dataFit_wrongDsUnmatchedPhoton->SetParameters(dataFitPar);
  f_dataFit_wrongDsUnmatchedPhoton->SetLineWidth(0);
  f_dataFit_wrongDsUnmatchedPhoton->Draw("SAME");
  TF1 *f_dataFit_wrongDsMatchedPhoton=new TF1("f_dataFit_wrongDsMatchedPhoton", dataFit_wrongDsMatchedPhoton, xmin, xmax, 24);
  f_dataFit_wrongDsMatchedPhoton->SetParameters(dataFitPar);
  f_dataFit_wrongDsMatchedPhoton->SetLineWidth(0);
  f_dataFit_wrongDsMatchedPhoton->Draw("SAME");
  ymax=(h_MBC_physics->GetMaximum())*0.95;
  line=new TLine(xmin, 0, xmin, ymax); line->Draw();
  line=new TLine(xmax, 0, xmax, ymax); line->Draw();
  filename_eps=decay;
  filename_png=decay;
  filename_eps+="_DsGammaData_MBC.eps";
  filename_png+="_DsGammaData_MBC.png";
  c_dataFit->SaveAs(filename_eps.c_str());
  c_dataFit->SaveAs(filename_png.c_str());
  double constant_data=dataFitPar[23]/dataFitPar[28];
  std::cout<<"constant_data = "<<constant_data<<std::endl;
  double nRegion=(f_dataFit->Integral(xmin, xmax))*(120.0/0.12);
  double nBackground=(f_dataFit_wrongDsMatchedPhoton->Integral(xmin, xmax))*(120.0/0.12);
  double nSignal=nRegion-nBackground;
  //double scaleRegion_errorFr=(f_dataFit->GetParError(28))/(f_dataFit->GetParameter(28));
  double scaleRegion_errorFr=(f_dataFit->GetParError(23))/(f_dataFit->GetParameter(23));
  double nSignal_error=nSignal*scaleRegion_errorFr;
  double branchingFraction=nSignal/(luminosity*prodCrossSection_DsDss*branchingFr_mode*nSignal_MC);
  double branchingFraction_error=branchingFraction*pow(pow(luminosity_error/luminosity, 2)
                                                      +pow(prodCrossSection_DsDss_error/prodCrossSection_DsDss, 2)
                                                      +pow(branchingFr_mode_error/branchingFr_mode, 2)
                                                      +pow(nSignal_MC_error/nSignal_MC, 2) 
                                                      +pow(scaleRegion_errorFr, 2), 0.5);
  double n_wrongDsUnmatchedPhoton_data=(f_dataFit_wrongDsUnmatchedPhoton->Integral(xmin, xmax))*(120./0.12);
  std::cout<<"ratio = "<<(nBackground-n_wrongDsUnmatchedPhoton_data)/(nRegion-nBackground)<<std::endl;
  
  std::cout<<"- Signal Efficiency -"<<std::endl;
  std::cout<<"Number of signal MC events in Region = "<<nRegion_MC<<std::endl;
  std::cout<<"Number of signal MC events as Combinatorics in the Region = "<<nCombo_MC<<std::endl;;
  std::cout<<"Number of signal MC events identified as signal = "<<nSignal_MC<<"+-"<<nSignal_MC_error<<std::endl;
  std::cout<<"---"<<std::endl;
  
  std::cout<<"- Data -"<<std::endl;
  std::cout<<"B(Ds->i) = "<<branchingFr_mode<<" +- "<<branchingFr_mode_error<<std::endl;
  std::cout<<"Number of events in Region = "<<nRegion<<std::endl;
  std::cout<<"Number of events as Combinatorics in the Region = "<<nBackground<<std::endl;;
  std::cout<<"Number of events identified as signal = "<<nSignal<<" +- "<<nSignal_error<<std::endl;
  std::cout<<"Scale error fraction = "<<scaleRegion_errorFr<<std::endl;
  std::cout<<"Branching fraction inferred = "<<branchingFraction<<" +- "<<branchingFraction_error<<std::endl;
  std::cout<<"TeX output: "<<std::endl;
  std::cout<<branchingFr_mode<<" $\pm$ "<<branchingFr_mode_error<<" & ";
  std::cout<<nSignal_MC<<" $\pm$ "<<nSignal_MC_error<<" & ";
  std::cout<<nSignal<<" $\pm$ "<<nSignal_error<<" & ";
  std::cout<<branchingFraction<<" $\pm$ "<<branchingFraction_error<<" \\\\ "<<std::endl;
  std::cout<<"---"<<std::endl;
  std::cout<<"- Generic MC -"<<std::endl;
  std::cout<<"B(Ds->i) generic = "<<branchingFr_mode_generic<<std::endl;
  std::cout<<"Number of events in Region = "<<nRegion_generic<<std::endl;
  std::cout<<"Number of events as Combinatorics in the Region = "<<nBackground_generic<<std::endl;
  std::cout<<"Number of events identified as signal = "<<nSignal_generic<<" +- "<<nSignal_generic_error<<std::endl;
  std::cout<<"Branching fraction inferred = "<<branchingFraction_generic<<" +- "<<branchingFraction_error_generic<<std::endl;
  std::cout<<"TeX output: "<<std::endl;
  std::cout<<branchingFr_mode_generic<<" & ";
  std::cout<<nSignal_MC<<" $\pm$ "<<nSignal_MC_error<<" & ";
  std::cout<<nSignal_generic<<" $\pm$ "<<nSignal_generic_error<<" & ";
  std::cout<<branchingFraction_generic<<" $\pm$ "<<branchingFraction_error_generic<<" \\\\ "<<std::endl;
  /*
  std::cout<<"---"<<std::endl;
  std::cout<<"nWrongSignBackground = "<<nWrongSignBackground<<std::endl;
  std::cout<<"wrongSign_eff = "<<wrongSign_eff<<std::endl;
  std::cout<<"Wrong sign prod = "<<wrongSignProd<<std::endl;
  std::cout<<"Branching fraction of generic from wrong sign = "<<branchingFraction_wrongSign<<" +- "<<branchingFraction_wrongSign_error<<std::endl;
  */
  
}
