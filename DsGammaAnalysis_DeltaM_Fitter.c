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

std::string decay="pietaprimerho";
// options: KKpi, KsK, pieta, pietaprime, KKpipi0, pipipi, KsKmpipi, pipi0eta, pietaprimerho
std::string decay_tex;

std::string texDsstgamma="D_{s}^{*#pm}#rightarrow D_{s}^{#pm}#gamma";
double deltaMCut_center, deltaMCut_range;
double deltam_x1, deltam_x2, deltam_y1, deltam_y2;
double conver_xmin=0.04, conver_xmax=0.2;
double xmin, xmax;

double luminosity=586; // /pb
double luminosity_error=6;
double prodCrossSection_DsDss=948;
double prodCrossSection_DsDss_error=36;
double branchingFr_mode;
double branchingFr_mode_error=0;

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

Double_t combo(Double_t *x, Double_t *par) // 2 Parameters
{
  Double_t result=(par[0]+par[1]*x[0]);
  return result;
}

Double_t combo2(Double_t *x, Double_t *par) // 3 Parameters
{
  Double_t result=(par[0]+par[1]*T_1(x[0])+par[2]*T_2(x[0]));
  //result+=par[3]*exp(-0.5*pow((x[0]-par[4])/par[5], 2));      // Now add the second Gaussian
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

Double_t converFit(Double_t *x, Double_t *par) // 9 Parameters
{
  Double_t result;
  
  result=crystalBall(x, par);
  result=result+combo(x, par+7);
  
  return result;
}

Double_t wrongConverFit(Double_t *x, Double_t *par) // 10 Parameters
{
  Double_t result;
  
  result=crystalBall(x, par);
  result=result+combo2(x, par+7);
  
  return result;
}

Double_t dataFit_wrongConver(Double_t *x, Double_t *par) // 2+10+1=13 Parameters
{
  Double_t result;
  result=combo(x, par); // 2 parameters
  result+=wrongConverFit(x, par+2)*par[12]; // 
  
  return result;
}

Double_t dataFit(Double_t *x, Double_t *par) // 7+13=20 Parameters
{
  Double_t result;
  result=dataFit_wrongConver(x, par);
  result+=crystalBall(x, par+13);
  return result;
}

void setValues()
{
  if (decay=="KKpi")
  {
    decay_tex="K^{+} K^{-} #pi^{#pm}";
    deltaMCut_center=0.14; deltaMCut_range=0.02; // widened.
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    branchingFr_mode=0.055;
    branchingFr_mode_error=0.0028;
  }
  else if (decay=="KsK")
  {
    decay_tex="K_{S}^{0} K^{#pm}";
    deltaMCut_center=0.14; deltaMCut_range=0.02;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    branchingFr_mode=0.0149;
    branchingFr_mode_error=0.0009;
  }
  else if (decay=="pieta")
  {
    decay_tex="#pi^{#pm} #eta, #eta #rightarrow #gamma #gamma";
    deltaMCut_center=0.14; deltaMCut_range=0.02;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    branchingFr_mode=0.0158*0.3931;
    branchingFr_mode_error=branchingFr_mode*pow(pow(0.0021/0.0158, 2)+pow(0.0020/0.3931, 2), 0.5);
  }
  else if (decay=="pietaprime")
  {
    deltaMCut_center=0.14; deltaMCut_range=0.02;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    branchingFr_mode=0.038*0.446*0.3931;
    branchingFr_mode=branchingFr_mode*pow(pow(0.004/0.038, 2)+pow(0.0014/0.446, 2)+pow(0.0020/0.3931, 2), 0.5);
  }
  else if (decay=="KKpipi0")
  {
    deltaMCut_center=0.14; deltaMCut_range=0.02;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    branchingFr_mode=0.056;
    branchingFr_mode_error=0.005;
  }
  else if (decay=="pipipi")
  {
    deltaMCut_center=0.14; deltaMCut_range=0.02;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    branchingFr_mode=0.0111;
    branchingFr_mode_error=0.0008;
  }
  else if (decay=="KsKmpipi")
  {
    deltaMCut_center=0.14; deltaMCut_range=0.02;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    branchingFr_mode=0.0164;
    branchingFr_mode_error=0.0012;
  }
  else if (decay=="pipi0eta")
  {
    deltaMCut_center=0.14; deltaMCut_range=0.02;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    branchingFr_mode=0.130*1.*0.3931;
    branchingFr_mode_error=branchingFr_mode*pow(pow(0.022/0.130, 2)+pow(0.0020/0.3931, 2), 0.5);
  }
  else if (decay=="pietaprimerho")
  {
    deltaMCut_center=0.14; deltaMCut_range=0.02;
    
    deltam_x1=0.65; deltam_y1=0.9;
    deltam_x2=0.9; deltam_y2=0.6;
    
    branchingFr_mode=0.038*0.294;
    branchingFr_mode_error=branchingFr_mode*pow(pow(0.004/0.038, 2)+pow(0.009/0.294, 2), 0.5);
  }
  xmin=deltaMCut_center-deltaMCut_range;
  xmax=deltaMCut_center+deltaMCut_range;
}

void DsGammaAnalysis_DeltaM_Fitter()
{

  setValues();
  
  std::string filename_root=decay;
  filename_root+="_DsGamma_DeltaM.root";
  
  TFile *file=new TFile(filename_root.c_str());
  
  TH1D *h_DeltaM_conver = (TH1D*)gDirectory->Get("h_DeltaM_conver");
  TH1D *h_DeltaM_generic = (TH1D*)gDirectory->Get("h_DeltaM_generic");
  TH1D *h_DeltaM_generic_veto = (TH1D*)gDirectory->Get("h_DeltaM_generic_veto");
  TH1D *h_DeltaM_continu = (TH1D*)gDirectory->Get("h_DeltaM_continu");
  TH1D *h_DeltaM_physics = (TH1D*)gDirectory->Get("h_DeltaM_physics");
  TH1D *h_DeltaM_wrongConver = (TH1D*)gDirectory->Get("h_DeltaM_wrongConver");
  
  h_DeltaM_generic->SetFillColor(kCyan);
  h_DeltaM_generic_veto->SetFillColor(kGreen);
  h_DeltaM_continu->SetFillColor(kBlue);
  
  h_DeltaM_generic_veto->Scale(1./19.2);
  
  double ymax;
  TLine *line;
  
  gROOT->SetStyle("Plain");
  
  // First fit the conver sample
  std::string title;
  title="#deltam Distribution in Signal Sample of ";
  title+=texDsstgamma;
  title+=", D_{s}^{#pm} #rightarrow ";
  title+=decay_tex;
  TCanvas *c_DeltaM_conver=new TCanvas("c_DeltaM_conver", "c_DeltaM_conver");
  h_DeltaM_conver->SetTitle(title.c_str());
  h_DeltaM_conver->GetYaxis()->SetTitle("Efficiency / MeV");
  h_DeltaM_conver->GetYaxis()->CenterTitle();
  h_DeltaM_conver->GetYaxis()->SetTitleOffset(1.3);
  h_DeltaM_conver->Draw();
  TF1 *f_converFit=new TF1("f_converFit", converFit, conver_xmin, conver_xmax, 9);
  f_converFit->SetParLimits(0, 0.143, 0.145);
  f_converFit->SetParLimits(1, 0.001, 0.1);
  f_converFit->SetParLimits(2, 1.0, 2.0);
  f_converFit->SetParLimits(3, 0.5, 5.0);
  f_converFit->SetParLimits(4, -3.0, -1.0);
  f_converFit->SetParLimits(5, 0.5, 5.0);
  f_converFit->SetLineWidth(0);
  h_DeltaM_conver->Fit(f_converFit, "R");
  TF1 *f_converFit_bg=new TF1("f_converFit_bg", combo, conver_xmin, conver_xmax, 2);
  Double_t converPar[9];
  f_converFit->GetParameters(converPar);
  f_converFit_bg->SetParameters(converPar+7);
  f_converFit_bg->SetLineWidth(0);
  f_converFit_bg->Draw("SAME");
  ymax=(h_DeltaM_conver->GetMaximum())*0.95;
  line=new TLine(xmin, 0, xmin, ymax); line->Draw();
  line=new TLine(xmax, 0, xmax, ymax); line->Draw();
  std::string filename_eps=decay;
  std::string filename_png=decay;
  filename_eps+="_DsGammaEff_DeltaM.eps";
  filename_png+="_DsGammaEff_DeltaM.png";
  c_DeltaM_conver->Print(filename_eps.c_str());
  c_DeltaM_conver->Print(filename_png.c_str());
  double nRegion_MC=(f_converFit->Integral(xmin, xmax))*(150.0/0.15);
  double nCombo_MC=(f_converFit_bg->Integral(xmin, xmax))*(150.0/0.15);
  double nSignal_MC=nRegion_MC-nCombo_MC;
  double nSignal_MC_error=pow(nSignal_MC/(nConversionSample_Dsp+nConversionSample_Dsm), 0.5);
  std::cout<<"- MC -"<<std::endl;
  std::cout<<"Number of signal MC events in Region = "<<nRegion_MC<<std::endl;
  std::cout<<"Number of signal MC events as Combinatorics in the Region = "<<nCombo_MC<<std::endl;;
  std::cout<<"Number of signal MC events identified as signal = "<<nSignal_MC<<"+-"<<nSignal_MC_error<<std::endl;
  std::cout<<"---"<<std::endl;
  
  // Second, fit the wrongConver sample
  title="#deltam Distribution in Wrong-Sign D_{s} #rightarrow ";
  title+=decay;
  TCanvas *c_wrongConver = new TCanvas("c_wrongConver", "c_wrongConver");
  h_DeltaM_wrongConver->SetTitle(title.c_str());
  h_DeltaM_wrongConver->GetYaxis()->SetTitle("# Events / MeV");
  h_DeltaM_wrongConver->GetXaxis()->SetTitle("#deltam (GeV)");
  h_DeltaM_wrongConver->GetYaxis()->CenterTitle();
  h_DeltaM_wrongConver->GetYaxis()->SetTitleOffset(1.2);
  h_DeltaM_wrongConver->Draw();
  TF1* f_wrongConverFit=new TF1("f_wrongConverFit", wrongConverFit, conver_xmin, conver_xmax, 10);
  f_wrongConverFit->SetParLimits(0, 0.139, 0.142);
  f_wrongConverFit->SetParLimits(1, 0.001, 0.1);
  f_wrongConverFit->SetParLimits(2, 1.0, 6.0);
  f_wrongConverFit->SetParLimits(3, 0.001, 2.0);
  f_wrongConverFit->SetParLimits(4, -3.0, -0.5);
  f_wrongConverFit->SetParLimits(5, 0.5, 7.0);
  //f_wrongConverFit->SetParLimits(11, 0.103, 0.113);
  f_wrongConverFit->SetLineWidth(0);
  h_DeltaM_wrongConver->Fit(f_wrongConverFit, "R");
  
  // Stack and Fit the Continuum, Generic veto Dsstgamma
  //TCanvas *c_DeltaM_Stacked = new TCanvas("c_DeltaM_Stacked");
  THStack *s_DeltaM_Background=new THStack("s_DeltaM_Background", "");
  s_DeltaM_Background->Add(h_DeltaM_continu, "hist");
  s_DeltaM_Background->Add(h_DeltaM_generic_veto, "hist");
  //s_DeltaM_Background->Draw();
  
  // Now fit the data
  TCanvas *c_dataFit = new TCanvas("dataFit", "dataFit");
  title="#deltam Distribution in Data for ";
  title+=texDsstgamma;
  title+=", D_{s} #rightarrow ";
  title+=decay;
  h_DeltaM_physics->SetTitle(title.c_str());
  h_DeltaM_physics->GetXaxis()->SetTitle("#deltam (GeV)");
  h_DeltaM_physics->GetYaxis()->SetTitle("# Events / MeV");
  h_DeltaM_physics->GetYaxis()->CenterTitle();
  h_DeltaM_physics->GetYaxis()->SetTitleOffset(1.2);
  //h_DeltaM_physics->Draw("SAME");
  TF1 *f_dataFit=new TF1("f_dataFit", dataFit, conver_xmin, conver_xmax, 20);
  // par 0 free
  // par 1 free
  f_dataFit->FixParameter(2, f_wrongConverFit->GetParameter(0)); // Could be moved //f_dataFit->SetParLimits(3, 2.110, 2.116);
  f_dataFit->FixParameter(3, f_wrongConverFit->GetParameter(1));
  f_dataFit->FixParameter(4, f_wrongConverFit->GetParameter(2));
  f_dataFit->FixParameter(5, f_wrongConverFit->GetParameter(3));
  f_dataFit->FixParameter(6, f_wrongConverFit->GetParameter(4));
  f_dataFit->FixParameter(7, f_wrongConverFit->GetParameter(5));
  f_dataFit->FixParameter(8, f_wrongConverFit->GetParameter(6));
  f_dataFit->FixParameter(9, f_wrongConverFit->GetParameter(7));
  f_dataFit->FixParameter(10, f_wrongConverFit->GetParameter(8));
  f_dataFit->FixParameter(11, f_wrongConverFit->GetParameter(9));
  // par 12 free - ratio of wrongConver to combo
  f_dataFit->SetParLimits(13, 0.143, 0.147);
  f_dataFit->FixParameter(14, f_converFit->GetParameter(1));
  f_dataFit->FixParameter(15, f_converFit->GetParameter(2));
  f_dataFit->FixParameter(16, f_converFit->GetParameter(3));
  f_dataFit->FixParameter(17, f_converFit->GetParameter(4));
  f_dataFit->FixParameter(18, f_converFit->GetParameter(5));
  // par 19 is movable - ratio of Crystal Ball shape to preceding backgrounds
  f_dataFit->SetLineWidth(0);
  h_DeltaM_physics->Fit(f_dataFit, "R");
  //h_DeltaM_generic->Draw("SAME");
  s_DeltaM_Background->Draw("SAME");
  Double_t dataFitPar[20];
  f_dataFit->GetParameters(dataFitPar);
  TF1 *f_dataFit_combo=new TF1("f_dataFit_combo", combo, conver_xmin, conver_xmax, 2);
  f_dataFit_combo->SetParameters(dataFitPar);
  f_dataFit_combo->SetLineWidth(0);
  f_dataFit_combo->Draw("SAME");
  TF1 *f_dataFit_wrongConver=new TF1("f_dataFit_wrongConver", dataFit_wrongConver, conver_xmin, conver_xmax, 13);
  f_dataFit_wrongConver->SetParameters(dataFitPar);
  f_dataFit_wrongConver->SetLineWidth(0);
  f_dataFit_wrongConver->Draw("SAME");
  ymax=(h_DeltaM_physics->GetMaximum())*0.95;
  line=new TLine(xmin, 0, xmin, ymax); line->Draw();
  line=new TLine(xmax, 0, xmax, ymax); line->Draw();
  double nRegion=(f_dataFit->Integral(xmin, xmax))*(120.0/0.12);
  double nBackground=(f_dataFit_wrongConver->Integral(xmin, xmax))*(120.0/0.12);
  double nSignal=nRegion-nBackground;
  double branchingFraction=nSignal/(luminosity*prodCrossSection_DsDss*branchingFr_mode*nSignal_MC);
  double branchingFraction_error=branchingFraction*pow(pow(luminosity_error/luminosity, 2)
                                                      +pow(prodCrossSection_DsDss_error/prodCrossSection_DsDss, 2)
                                                      +pow(branchingFr_mode_error/branchingFr_mode, 2)
                                                      +pow(nSignal_MC_error/nSignal_MC, 2) ,0.5);
  std::cout<<"- Data -"<<std::endl;
  std::cout<<"Number of events in Region = "<<nRegion<<std::endl;
  std::cout<<"Number of events as Combinatorics in the Region = "<<nBackground<<std::endl;;
  std::cout<<"Number of events identified as signal = "<<nSignal<<std::endl;
  std::cout<<"Branching fraction inferred = "<<branchingFraction<<" +- "<<branchingFraction_error<<std::endl;
  std::cout<<"---"<<std::endl;
  
  // Now fit the generic MC
  TCanvas *c_genericFit = new TCanvas("genericFit", "genericFit");
  title="#deltam Distribution in Generic MC for ";
  title+=texDsstgamma;
  title+=", D_{s} #rightarrow ";
  title+=decay;
  h_DeltaM_generic->SetTitle(title.c_str());
  h_DeltaM_generic->GetXaxis()->SetTitle("#deltam (GeV)");
  h_DeltaM_generic->GetYaxis()->SetTitle("# Events / MeV");
  h_DeltaM_generic->GetYaxis()->CenterTitle();
  h_DeltaM_generic->GetYaxis()->SetTitleOffset(1.2);
  TF1 *f_genericFit=new TF1("f_genericFit", dataFit, conver_xmin, conver_xmax, 20);
  // par 0 free
  // par 1 free
  f_genericFit->FixParameter(2, f_wrongConverFit->GetParameter(0)); // Could be moved //f_dataFit->SetParLimits(3, 2.110, 2.116);
  f_genericFit->FixParameter(3, f_wrongConverFit->GetParameter(1));
  f_genericFit->FixParameter(4, f_wrongConverFit->GetParameter(2));
  f_genericFit->FixParameter(5, f_wrongConverFit->GetParameter(3));
  f_genericFit->FixParameter(6, f_wrongConverFit->GetParameter(4));
  f_genericFit->FixParameter(7, f_wrongConverFit->GetParameter(5));
  f_genericFit->FixParameter(8, f_wrongConverFit->GetParameter(6));
  f_genericFit->FixParameter(9, f_wrongConverFit->GetParameter(7));
  f_genericFit->FixParameter(10, f_wrongConverFit->GetParameter(8));
  f_genericFit->FixParameter(11, f_wrongConverFit->GetParameter(9));
  // par 12 free - ratio of wrongConver to combo
  f_genericFit->SetParLimits(13, 0.143, 0.147);
  f_genericFit->FixParameter(14, f_converFit->GetParameter(1));
  f_genericFit->FixParameter(15, f_converFit->GetParameter(2));
  f_genericFit->FixParameter(16, f_converFit->GetParameter(3));
  f_genericFit->FixParameter(17, f_converFit->GetParameter(4));
  f_genericFit->FixParameter(18, f_converFit->GetParameter(5));
  // par 19 is movable - ratio of Crystal Ball shape to preceding backgrounds
  f_genericFit->SetLineWidth(0);
  h_DeltaM_generic->Fit(f_genericFit, "R");
  s_DeltaM_Background->Draw("SAME");
  Double_t genericFitPar[20];
  f_genericFit->GetParameters(genericFitPar);
  TF1 *f_genericFit_combo=new TF1("f_genericFit_combo", combo, conver_xmin, conver_xmax, 2);
  f_genericFit_combo->SetParameters(genericFitPar);
  f_genericFit_combo->SetLineWidth(0);
  f_genericFit_combo->Draw("SAME");
  TF1 *f_genericFit_wrongConver=new TF1("f_genericFit_wrongConver", dataFit_wrongConver, conver_xmin, conver_xmax, 13);
  f_genericFit_wrongConver->SetParameters(genericFitPar);
  f_genericFit_wrongConver->SetLineWidth(0);
  f_genericFit_wrongConver->Draw("SAME");
  ymax=(h_DeltaM_generic->GetMaximum())*0.95;
  line=new TLine(xmin, 0, xmin, ymax); line->Draw();
  line=new TLine(xmax, 0, xmax, ymax); line->Draw();
  double nRegion_generic=(f_genericFit->Integral(xmin, xmax))*(120.0/0.12);
  double nBackground_generic=(f_genericFit_wrongConver->Integral(xmin, xmax))*(120.0/0.12);
  double nSignal_generic=nRegion_generic-nBackground_generic;
  double branchingFraction_generic=nSignal_generic/(luminosity*prodCrossSection_DsDss*branchingFr_mode*nSignal_MC);
  double branchingFraction_error_generic=branchingFraction_generic*pow(pow(luminosity_error/luminosity, 2)
                                                      +pow(prodCrossSection_DsDss_error/prodCrossSection_DsDss, 2)
                                                      +pow(branchingFr_mode_error/branchingFr_mode, 2)
                                                      +pow(nSignal_MC_error/nSignal_MC, 2) ,0.5);
  std::cout<<"- Generic MC -"<<std::endl;
  std::cout<<"Number of events in Region_generic = "<<nRegion_generic<<std::endl;
  std::cout<<"Number of events as Combinatorics in the Region_generic = "<<nBackground_generic<<std::endl;;
  std::cout<<"Number of events identified as signal_generic = "<<nSignal_generic<<std::endl;
  std::cout<<"Branching fraction inferred_generic = "<<branchingFraction_generic<<" +- "<<branchingFraction_error_generic<<std::endl;
  std::cout<<"---"<<std::endl;
  
}
