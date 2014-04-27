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

std::string decay="KKpi";
// options: KKpi, KsK, pieta, pietaprime, KKpipi0, pipipi, KsKmpipi, pipi0eta, pietaprimerho

std::string texDsstgamma="D_{s}^{*} #to D_{s} #gamma";
double mbcCut_center, mbcCut_range;
double mbc_x1, mbc_x2, mbc_y1, mbc_y2;
double conver_xmin, conver_xmax, physics_xmin, physics_xmax;
double xmin, xmax;

double luminosity=586; // /pb
double luminosity_error=6;
double prodCrossSection_DsDss=948;
double prodCrossSection_DsDss_error=36;
double branchingFr_mode;
double branchingFr_mode_error=0;

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

Double_t argus(Double_t *x, Double_t *par) // 3 Parameters
{
  Double_t result=0;
  if (par[0]>x[0]) result=(par[1]+par[2]*x[0])*pow(par[0]-x[0], 0.5);
  return result;
}

Double_t crystalBall(Double_t *x, Double_t *par) // 7 Parameters
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

Double_t converFit(Double_t *x, Double_t *par) // 10 Parameters
{
  Double_t result;
  
  result=crystalBall(x, par);
  result=result+argus(x, par+7);
  
  return result;
}

Double_t wrongConverFit(Double_t *x, Double_t *par) // 14 Parameters // 18 Parameters now
{
  Double_t result;
  
  result=crystalBall(x, par); // 7 parameters
  
  /*
  // Now add the second Gaussian
  result+=par[7]*exp(-0.5*pow((x[0]-par[8])/par[9], 2));
  */
  // Add a second CB
  result+=crystalBall(x, par+7); // 14 parameters
  
  // Now add the background function
  /*
  Double_t width=(x[0]-par[10])/par[11];
  if (x[0]>par[10]) result=result+exp(-pow(width, 2))*pow(width, par[12])*par[13];
  */
  Double_t width=(x[0]-par[14])/par[15];
  if (x[0]>par[14]) result=result+exp(-pow(width, 2))*pow(width, par[16])*par[17];
  
  return result;
}

Double_t dataFit_wrongConver(Double_t *x, Double_t *par) // 18 Parameters // now 22
{
  Double_t result;
  result=argus(x, par);
  result+=wrongConverFit(x, par+3)*par[21];
  
  return result;
}

Double_t dataFit(Double_t *x, Double_t *par) // 25 Parameters
{
  Double_t result;
  result=dataFit_wrongConver(x, par);
  result+=crystalBall(x, par+22);
  return result;
}

void setValues()
{
  if (decay=="KKpi")
  {
    mbcCut_center=2.112; mbcCut_range=0.004;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
    conver_xmin=2.05; conver_xmax=2.150;
    physics_xmin=2.06; physics_xmax=2.155;
    
    branchingFr_mode=0.055;
    branchingFr_mode_error=0.0028;
  }
  else if (decay=="KsK")
  {
    mbcCut_center=2.112; mbcCut_range=0.007;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.35; mbc_y2=0.6;
    
    conver_xmin=2.06; conver_xmax=2.155;
    physics_xmin=2.02; physics_xmax=2.160;
    
    branchingFr_mode=0.0149;
    branchingFr_mode_error=0.0009;
  }
  else if (decay=="pieta")
  {
    mbcCut_center=2.112; mbcCut_range=0.008;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
    conver_xmin=2.06; conver_xmax=2.155;
    physics_xmin=2.02; physics_xmax=2.160;
    
    branchingFr_mode=0.0158*0.3931;
    branchingFr_mode_error=branchingFr_mode*pow(pow(0.0021/0.0158, 2)+pow(0.0020/0.3931, 2), 0.5);
  }
  else if (decay=="pietaprime")
  {
    mbcCut_center=2.112; mbcCut_range=0.011;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
    conver_xmin=2.06; conver_xmax=2.155;
    physics_xmin=2.02; physics_xmax=2.160;
    
    branchingFr_mode=0.038*0.446*0.3931;
    branchingFr_mode=branchingFr_mode*pow(pow(0.004/0.038, 2)+pow(0.0014/0.446, 2)+pow(0.0020/0.3931, 2), 0.5);
  }
  else if (decay=="KKpipi0")
  {
    mbcCut_center=2.112; mbcCut_range=0.004;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
    conver_xmin=2.06; conver_xmax=2.155;
    physics_xmin=2.02; physics_xmax=2.160;
    
    branchingFr_mode=0.056;
    branchingFr_mode_error=0.005;
  }
  else if (decay=="pipipi")
  {
    mbcCut_center=2.112; mbcCut_range=0.004;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
    conver_xmin=2.06; conver_xmax=2.155;
    physics_xmin=2.02; physics_xmax=2.160;
    
    branchingFr_mode=0.0111;
    branchingFr_mode_error=0.0008;
  }
  else if (decay=="KsKmpipi")
  {
    mbcCut_center=2.112; mbcCut_range=0.005;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
    conver_xmin=2.06; conver_xmax=2.155;
    physics_xmin=2.02; physics_xmax=2.160;
    
    branchingFr_mode=0.0164;
    branchingFr_mode_error=0.0012;
  }
  else if (decay=="pipi0eta")
  {
    mbcCut_center=2.112; mbcCut_range=0.004;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
    conver_xmin=2.06; conver_xmax=2.155;
    physics_xmin=2.02; physics_xmax=2.160;
    
    branchingFr_mode=0.130*1.*0.3931;
    branchingFr_mode_error=branchingFr_mode*pow(pow(0.022/0.130, 2)+pow(0.0020/0.3931, 2), 0.5);
  }
  else if (decay=="pietaprimerho")
  {
    mbcCut_center=2.112; mbcCut_range=0.004;
    
    mbc_x1=0.15; mbc_y1=0.9;
    mbc_x2=0.42; mbc_y2=0.6;
    
    conver_xmin=2.06; conver_xmax=2.155;
    physics_xmin=2.02; physics_xmax=2.160;
    
    branchingFr_mode=0.038*0.294;
    branchingFr_mode_error=branchingFr_mode*pow(pow(0.004/0.038, 2)+pow(0.009/0.294, 2), 0.5);
  }
  xmin=mbcCut_center-mbcCut_range;
  xmax=mbcCut_center+mbcCut_range;
}

void DsGammaAnalysis_Fitter()
{

  setValues();
  
  std::string filename;
  filename=decay;
  filename+="_DsGamma_MBC.root";
  TFile *file=new TFile("KKpi_DsGamma_MBC.root");
  
  TH1D *h_MBC_conver = (TH1D*)gDirectory->Get("h_MBC_conver");
  TH1D *h_MBC_generic = (TH1D*)gDirectory->Get("h_MBC_generic");
  TH1D *h_MBC_generic_veto = (TH1D*)gDirectory->Get("h_MBC_generic_veto");
  TH1D *h_MBC_continu = (TH1D*)gDirectory->Get("h_MBC_continu");
  TH1D *h_MBC_physics = (TH1D*)gDirectory->Get("h_MBC_physics");
  TH1D *h_MBC_wrongConver = (TH1D*)gDirectory->Get("h_MBC_wrongConver");
  
  h_MBC_generic->SetFillColor(kGreen);
  h_MBC_generic_veto->SetFillColor(kGreen);
  h_MBC_continu->SetFillColor(kBlue);
  
  double ymax;
  TLine *line;
  
  gROOT->SetStyle("Plain");
  
  // First fit the conver sample
  std::string title;
  title="m_{BC} Distribution in Signal Sample of ";
  title+=texDsstgamma;
  title+=", D_{s} #to ";
  title+=decay;
  TCanvas *c_MBC_conver=new TCanvas("c_MBC_conver", "c_MBC_conver");
  h_MBC_conver->SetTitle(title.c_str());
  h_MBC_conver->Draw();
  TF1 *f_converFit=new TF1("f_converFit", converFit, conver_xmin, conver_xmax, 10);
  f_converFit->SetParLimits(0, 2.111, 2.113);
  f_converFit->SetParLimits(1, 0.001, 0.1);
  f_converFit->SetParLimits(2, 1.0, 1.5);
  f_converFit->SetParLimits(3, 0.5, 3.0);
  f_converFit->SetParLimits(4, -3.0, -1.5);
  f_converFit->SetParLimits(5, 0.5, 5.0);
  f_converFit->SetParLimits(7, 2.150, 2.155);
  f_converFit->SetLineWidth(0);
  h_MBC_conver->Fit(f_converFit, "REFM");
  TF1 *f_converFit_bg=new TF1("f_converFit_bg", argus, conver_xmin, conver_xmax, 3);
  Double_t converPar[10];
  f_converFit->GetParameters(converPar);
  f_converFit_bg->SetParameters(converPar+7);
  f_converFit_bg->SetLineWidth(0);
  f_converFit_bg->Draw("SAME");
  ymax=(h_MBC_conver->GetMaximum())*0.95;
  line=new TLine(xmin, 0, xmin, ymax); line->Draw();
  line=new TLine(xmax, 0, xmax, ymax); line->Draw();
  double nRegion_MC=(f_converFit->Integral(xmin, xmax))*(120.0/0.12);
  double nCombo_MC=(f_converFit_bg->Integral(xmin, xmax))*(120.0/0.12);
  double nSignal_MC=nRegion_MC-nCombo_MC;
  std::cout<<"- MC -"<<std::endl;
  std::cout<<"Number of signal MC events in Region = "<<nRegion_MC<<std::endl;
  std::cout<<"Number of signal MC events as Combinatorics in the Region = "<<nCombo_MC<<std::endl;;
  std::cout<<"Number of signal MC events identified as signal = "<<nSignal_MC<<std::endl;
  std::cout<<"---"<<std::endl;
  
  // Second, fit the wrongConver sample
  TCanvas *c_wrongConver = new TCanvas("c_wrongConver", "c_wrongConver");
  h_MBC_wrongConver->Draw();
  TF1* f_wrongConverFit=new TF1("f_wrongConverFit", wrongConverFit, conver_xmin, conver_xmax, 18); // 14
  f_wrongConverFit->SetParLimits(0, 2.110, 2.114);
  f_wrongConverFit->SetParLimits(1, 0.001, 0.1);
  f_wrongConverFit->SetParLimits(2, 1.0, 5.0);
  f_wrongConverFit->SetParLimits(3, 0.1, 2.0);
  f_wrongConverFit->SetParLimits(4, -3.0, -0.1);
  f_wrongConverFit->SetParLimits(5, 1.0, 5.0);
  f_wrongConverFit->SetParLimits(7, 2.127, 2.135); // 2nd CB
  f_wrongConverFit->SetParLimits(8, 0.001, 0.01);
  f_wrongConverFit->SetParLimits(9, 1.0, 5.0);
  f_wrongConverFit->SetParLimits(10, 0.1, 2.0);
  f_wrongConverFit->SetParLimits(11, -3.0, -0.05);
  f_wrongConverFit->SetParLimits(12, 0.01, 5.0);
  f_wrongConverFit->SetParLimits(14, 2.05, 2.07);
  f_wrongConverFit->SetParLimits(15, 0.005, 0.19);
  f_wrongConverFit->SetParLimits(16, 0.5, 1.5);
  /*
  f_wrongConverFit->SetParLimits(0, 2.110, 2.114);
  f_wrongConverFit->SetParLimits(1, 0.001, 0.1);
  f_wrongConverFit->SetParLimits(2, 1.0, 5.0);
  f_wrongConverFit->SetParLimits(3, 0.1, 2.0);
  f_wrongConverFit->SetParLimits(4, -3.0, -0.1);
  f_wrongConverFit->SetParLimits(5, 1.0, 5.0);
  f_wrongConverFit->SetParLimits(8, 2.127, 2.135);
  f_wrongConverFit->SetParLimits(9, 0.001, 0.01);
  f_wrongConverFit->SetParLimits(10, 2.05, 2.07);
  f_wrongConverFit->SetParLimits(11, 0.005, 0.19);
  f_wrongConverFit->SetParLimits(12, 0.5, 1.5);
  */
  
  f_wrongConverFit->SetLineWidth(0);
  h_MBC_wrongConver->Fit(f_wrongConverFit, "R");
  
  
  // Now fit the data
  TCanvas *c_dataFit = new TCanvas("dataFit", "dataFit");
  TF1 *f_dataFit=new TF1("f_dataFit", dataFit, physics_xmin, physics_xmax, 29); // 25
  f_dataFit->SetParLimits(0, 2.153, 2.156); // right limit of Argus function
  // par 1 movable
  // par 2 movable
  f_dataFit->FixParameter(3, f_wrongConverFit->GetParameter(0)); // Could be moved //f_dataFit->SetParLimits(3, 2.110, 2.116);
  f_dataFit->FixParameter(4, f_wrongConverFit->GetParameter(1));
  f_dataFit->FixParameter(5, f_wrongConverFit->GetParameter(2));
  f_dataFit->FixParameter(6, f_wrongConverFit->GetParameter(3));
  f_dataFit->FixParameter(7, f_wrongConverFit->GetParameter(4));
  f_dataFit->FixParameter(8, f_wrongConverFit->GetParameter(5));
  f_dataFit->FixParameter(9, f_wrongConverFit->GetParameter(6));
  f_dataFit->FixParameter(10, f_wrongConverFit->GetParameter(7));
  f_dataFit->FixParameter(11, f_wrongConverFit->GetParameter(8));
  f_dataFit->FixParameter(12, f_wrongConverFit->GetParameter(9));
  f_dataFit->FixParameter(13, f_wrongConverFit->GetParameter(10));
  f_dataFit->FixParameter(14, f_wrongConverFit->GetParameter(11));
  f_dataFit->FixParameter(15, f_wrongConverFit->GetParameter(12));
  f_dataFit->FixParameter(16, f_wrongConverFit->GetParameter(13));
  f_dataFit->FixParameter(17, f_wrongConverFit->GetParameter(14));
  f_dataFit->FixParameter(18, f_wrongConverFit->GetParameter(15));
  f_dataFit->FixParameter(19, f_wrongConverFit->GetParameter(16));
  f_dataFit->FixParameter(20, f_wrongConverFit->GetParameter(17));
  // par 21 is movable - ratio of wrongConver to argus function
  f_dataFit->SetParLimits(22, 2.110, 2.113);
  f_dataFit->FixParameter(23, f_converFit->GetParameter(1));
  f_dataFit->FixParameter(24, f_converFit->GetParameter(2));
  f_dataFit->FixParameter(25, f_converFit->GetParameter(3));
  f_dataFit->FixParameter(26, f_converFit->GetParameter(4));
  f_dataFit->FixParameter(27, f_converFit->GetParameter(5));
  // par 28 is movable - ratio of Crystal Ball shape to preceding backgrounds
  f_dataFit->SetLineWidth(0);
  h_MBC_physics->Fit(f_dataFit, "R");
  /*
  Double_t dataFitPar[25];
  f_dataFit->GetParameters(dataFitPar);
  TF1 *f_dataFit_argus=new TF1("f_dataFit_argus", argus, physics_xmin, physics_xmax, 3);
  f_dataFit_argus->SetParameters(dataFitPar);
  f_dataFit_argus->SetLineWidth(0);
  f_dataFit_argus->Draw("SAME");
  TF1 *f_dataFit_wrongConver=new TF1("f_dataFit_wrongConver", dataFit_wrongConver, physics_xmin, physics_xmax, 18);
  f_dataFit_wrongConver->SetParameters(dataFitPar);
  f_dataFit_wrongConver->SetLineWidth(0);
  f_dataFit_wrongConver->Draw("SAME");
  ymax=(h_MBC_physics->GetMaximum())*0.95;
  line=new TLine(xmin, 0, xmin, ymax); line->Draw();
  line=new TLine(xmax, 0, xmax, ymax); line->Draw();
  double nRegion=(f_dataFit->Integral(xmin, xmax))*(120.0/0.12);
  double nBackground=(f_dataFit_wrongConver->Integral(xmin, xmax))*(120.0/0.12);
  double nSignal=nRegion-nBackground;
  double branchingFraction=nSignal/(luminosity*prodCrossSection_DsDss*branchingFr_mode*nSignal_MC);
  double branchingFraction_error=branchingFraction*pow(pow(luminosity_error/luminosity, 2)
                                                      +pow(prodCrossSection_DsDss_error/prodCrossSection_DsDss, 2)
                                                      +pow(branchingFr_mode_error/branchingFr_mode, 2) ,0.5); // Ignoring error of efficiency
  std::cout<<"- MC -"<<std::endl;
  std::cout<<"Number of events in Region = "<<nRegion<<std::endl;
  std::cout<<"Number of events as Combinatorics in the Region = "<<nBackground<<std::endl;;
  std::cout<<"Number of events identified as signal = "<<nSignal<<std::endl;
  std::cout<<"Branching fraction inferred = "<<branchingFraction<<" +- "<<branchingFraction_error<<std::endl;
  std::cout<<"---"<<std::endl;
  
  
  // Now draw the Background and Fitted Data
  title="m_{BC} Distribution in Monte Carlo Backgrounds and Data in ";
  title+=decay;
  TCanvas *c_MBC_Stacked = new TCanvas("c_MBC_Stacked");
  THStack *s_MBC_Background=new THStack("s_MBC_Background", "");
  s_MBC_Background->Add(h_MBC_continu, "hist");
  s_MBC_Background->Add(h_MBC_generic_veto, "hist");
  double stack_ymax=h_MBC_generic->GetMaximum();
  double physics_ymax=h_MBC_physics->GetMaximum();
  if (physics_ymax>stack_ymax)
  {
    s_MBC_Background->SetMaximum(physics_ymax*1.);
    ymax=physics_ymax;
  }
  else
  {
    s_MBC_Background->SetMaximum(stack_ymax*1.);
    ymax=stack_ymax;
  }
  s_MBC_Background->Draw();
  h_MBC_physics->Draw("SAME");
  f_dataFit_argus->Draw("SAME");
  f_dataFit_wrongConver->Draw("SAME");
  s_MBC_Background->SetTitle(title.c_str());
  s_MBC_Background->GetXaxis()->SetTitle("m_{BC} (GeV)");
  s_MBC_Background->GetYaxis()->SetTitle("Number of Events");
  
  line=new TLine(xmin, 0, xmin, ymax); line->Draw();
  line=new TLine(xmax, 0, xmax, ymax); line->Draw();
  //TLegend *legendMBC=new TLegend(mbc_x1, mbc_y1, mbc_x2, mbc_y2);
  //std::string continu_string="Continuum MC: "; continu_string+=itoa(continuScale*h_MBC_continu->GetEntries()); continu_string+=" Events";
  //std::string generic_string="Generic MC: "; generic_string+=itoa(genericScale*h_MBC_generic->GetEntries()); generic_string+=" Events";
  //std::string physics_string="Data: "; physics_string+=itoa(h_MBC_physics->GetEntries()); physics_string+=" Events";
  //legendMBC->AddEntry(h_MBC_continu, continu_string.c_str());
  //legendMBC->AddEntry(h_MBC_generic, generic_string.c_str());
  //legendMBC->AddEntry(h_MBC_physics, physics_string.c_str());
  //legendMBC->SetFillColor(kWhite);
  //legendMBC->Draw();
  
  */
}
