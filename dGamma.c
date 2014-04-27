#include <iostream>
#include <math.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TText.h"
#include "TLegend.h"
#include "TStyle.h"

using namespace std;

const double m_e=0.000511;
const double m_e2=m_e*m_e;
const double m_Ds=1.969;
const double m_Ds2=m_Ds*m_Ds;
const double m_Dsstr=2.112;
const double m_Dsstr2=m_Dsstr*m_Dsstr;
const double m_rho=0.77549;
const double m_omega=0.782;
const double m_phi=1.019455;

const int nq2=1000000;

const double q2min=4*m_e*m_e;
const double q2max=(m_Dsstr-m_Ds)*(m_Dsstr-m_Ds);
const double deltaq2=(q2max-q2min)/nq2;

const double alpha=1.0/137.0;
const double pi=4.0*atan(1.0);

const double Egamma=(m_Dsstr2-m_Ds2)/(2*m_Dsstr);

double C(double q2)
{
  //return 1./(1.+q2/(m_omega*m_omega));
  return 1./(1.-q2/(m_phi*m_phi));
  // return 1.;
}

double f(double q2){

  double EDs=(m_Dsstr2+m_Ds2-q2)/(2.0*m_Dsstr);
  double PDs=sqrt(EDs*EDs-m_Ds2);
  double A=m_Dsstr2-m_Ds2+q2;
  
  //double prefact=PDs*alpha*C*C/(48*pi*pi*m_Dsstr2);
  double prefact=PDs*pow(1-4*m_e2/q2, 0.5)*alpha*C(q2)*C(q2)/(48*pi*pi*m_Dsstr2);
	/*
  return prefact*(
		  (2.0/3.0)*A*A*m_e2/(q2*q2)
		  +(2.0*A*A/6.0+4*m_e2*m_Dsstr2)/q2
                  -(4.0/3.0)*m_Dsstr2
		  );
  */
  
  return prefact*(
		  (2.0/3.0)*A*A*m_e2/(q2*q2)
		  +(A*A/3.0-(8.0/3.0)*m_e2*m_Dsstr2)/q2
                  -(4.0/3.0)*m_Dsstr2
	);
  
}
/*
double fsimple(double q2){

  static bool first=true;

  double EDs=(m_Dsstr2+m_Ds2-q2)/(2.0*m_Dsstr);
  double PDs=sqrt(EDs*EDs-m_Ds2);
  double A=m_Dsstr2-m_Ds2+q2;
  
  if (first) {
    cout << "PDs="<<PDs<<endl;
    first=false;
  }

  double prefact=PDs*alpha*C*C/(3*pi*pi*16*m_Dsstr2*q2*q2);
	    
  return prefact*(0.75*A*A*q2);

}
*/

void dGamma()
{
  bool zoom=true;
  
  double xmin=0; // 2/nq2;
  double xmax;
  
  if (zoom) xmax=5*q2min;
  else xmax=0.0002;

 TChain *signalTree=new TChain("DsTaggedDecaysProc/nt8");
 signalTree->Add("/nfs/cor/an2/souvik/MC_vtosll_Dsp_KKpi/DsTaggedDecaysProc_MC_vtosll_Dsp_KKpi.root");
 signalTree->Add("/nfs/cor/an2/souvik/MC_vtosll_Dsm_KKpi/DsTaggedDecaysProc_MC_vtosll_Dsm_KKpi.root");
 TH1D *h_meeMC = new TH1D("h_meeMC", "k^{2} Distribution; m_{ee}^{2} (GeV^{2}); d#Gamma/dk^{2}", (xmax-xmin)/(2*deltaq2), xmin, xmax); h_meeMC->SetLineColor(kRed); h_meeMC->Sumw2();
 TH1D *h_mee = new TH1D("h_mee", "k^{2} Distribution; m_{ee}^{2} (GeV^{2}); d#Gamma/dk^{2}", (xmax-xmin)/(2*deltaq2), xmin, xmax); h_mee->SetLineColor(kGreen);
 TH1D *h_meeMC_recalc = new TH1D("h_meeMC_recalc", "k^{2} Distribution; m_{ee}^{2} (GeV^{2}); d#Gamma/dk^{2}", (xmax-xmin)/(2*deltaq2), xmin, xmax); h_meeMC_recalc->SetLineColor(kBlue); h_meeMC_recalc->Sumw2();
 TH1D *h_ratio = new TH1D("h_ratio", "h_ratio", 100, xmin, xmax);
 TH1D *h_me = new TH1D("h_me", "Invariant Mass Distribution of the Electron from Generator Level Monte Carlo; m_{e} (GeV)", 100, 0., 5*.000511);
 TH2D *h_Ee = new TH2D("h_Ee", "h_Ee", 100, 0., 0.140, 100, 0., 0.140);
 TH2D *h_Ee_me = new TH2D("h_Ee_me", "Invariant Mass vs Energy of Electron from Generator Level Monte Carlo; E (GeV); m_{e} (GeV)", 100, 0., 0.140, 100, 0., 5*.000511);
 
 float mee_signalMC;
 float electron_E, electron_Px, electron_Py, electron_Pz;
 float positron_E, positron_Px, positron_Py, positron_Pz;
 
 signalTree->SetBranchAddress("keeMass_MC", &(mee_signalMC));
 signalTree->SetBranchAddress("kElectron1E_MC", &(electron_E));
 signalTree->SetBranchAddress("kElectron1Px_MC", &(electron_Px));
 signalTree->SetBranchAddress("kElectron1Py_MC", &(electron_Py));
 signalTree->SetBranchAddress("kElectron1Pz_MC", &(electron_Pz));
 signalTree->SetBranchAddress("kElectron2E_MC", &(positron_E));
 signalTree->SetBranchAddress("kElectron2Px_MC", &(positron_Px));
 signalTree->SetBranchAddress("kElectron2Py_MC", &(positron_Py));
 signalTree->SetBranchAddress("kElectron2Pz_MC", &(positron_Pz));
 
 int nEntries=signalTree->GetEntries();
 for (int i=0; i<nEntries; ++i)
 {
  signalTree->GetEvent(i);
  h_meeMC->Fill(mee_signalMC);
  
  double me=pow(electron_E*electron_E
               -electron_Px*electron_Px
               -electron_Py*electron_Py
               -electron_Pz*electron_Pz, 0.5);
               
  double mp=pow(positron_E*positron_E
               -positron_Px*positron_Px
               -positron_Py*positron_Py
               -positron_Pz*positron_Pz, 0.5);
               
  h_me->Fill(me);
  h_me->Fill(mp);
  
  double Ee=pow(electron_Px*electron_Px
               +electron_Py*electron_Py
               +electron_Pz*electron_Pz
               +.000511*.000511, 0.5);
               
  double Ep=pow(positron_Px*positron_Px
               +positron_Py*positron_Py
               +positron_Pz*positron_Pz
               +.000511*.000511, 0.5);
               
  h_Ee->Fill(Ee, electron_E);
  h_Ee->Fill(Ep, positron_E);
  
  double mee_recalc=(Ee+Ep)*(Ee+Ep)
                   -(electron_Px+positron_Px)*(electron_Px+positron_Px)
                   -(electron_Py+positron_Py)*(electron_Py+positron_Py)
                   -(electron_Pz+positron_Pz)*(electron_Pz+positron_Pz);
                       
  h_meeMC_recalc->Fill(mee_recalc);
  //std::cout<<mee_recalc<<", ";
  
  h_Ee_me->Fill(electron_E, me);
  h_Ee_me->Fill(positron_E, mp);
  
 }
 
 std::cout<<"q2min = "<<q2min<<std::endl;
 
  double sum=0.0;

  for(int iq2=0; iq2<nq2 ; iq2++ ) {
    double q2=q2min+(iq2+0.5)*deltaq2;
    sum+=f(q2);
    h_mee->Fill(q2, f(q2));
  }

  double Gammaee=sum*deltaq2;

  double Gammagamma=pow(Egamma,3)/(12.0*pi); // C=1

  cout << "pi :" << pi << endl;
  cout << "C  :" << C << endl;
  cout << "m_e:" << m_e << endl;
  cout << "Egamma:" << Egamma << endl;

  cout << "nq2:" << nq2 << endl;

  cout << "sum:" << sum << endl;

  cout << "Gammagamma:"<<Gammagamma<<endl;
  cout << "Gammaee:"   <<Gammaee<<endl;
  
  cout << "Gammaee/Gammagamma:"<<Gammaee/Gammagamma<<endl;
  
  h_mee->Scale(1/h_mee->Integral(h_mee->FindBin(q2min), h_mee->FindBin(q2max)));
  h_meeMC->Scale(1/h_meeMC->GetEntries());
  h_meeMC_recalc->Scale(1/h_meeMC_recalc->GetEntries());
  
  std::cout<<"FindBin(q2min) = "<<h_mee->FindBin(q2min)<<std::endl;
  std::cout<<"FindBin(q2max) = "<<h_mee->FindBin(q2max)<<std::endl;
  
  h_ratio->Divide(h_mee, h_meeMC_recalc);
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  TCanvas *c_mee_problem = new TCanvas("c_mee_problem", "c_mee_problem", 500, 500);
  if (!zoom) 
  {
    gPad->SetLogx(1);
    gPad->SetLeftMargin(.14);
    h_mee->GetYaxis()->SetTitleOffset(1.7);
  }
  else 
  {
    gPad->SetLeftMargin(0.16);
    h_mee->GetYaxis()->SetTitleOffset(2.1);
  }
  h_mee->GetYaxis()->CenterTitle();
  h_mee->GetXaxis()->SetTitleOffset(1.1);
  h_mee->GetXaxis()->CenterTitle();
  if (!zoom) h_mee->SetMaximum(0.09);
  else h_mee->SetMaximum(0.0035);
  h_meeMC->SetMarkerStyle(7);
  h_mee->Draw("C");
  h_meeMC->Draw("SAME");
  TLegend *legend;
  if (!zoom) legend=new TLegend(.5, .9, .9, .7);
  else legend=new TLegend(.5, .9, .9, .7);
  legend->AddEntry(h_mee, "Analytical Distribution");
  legend->AddEntry(h_meeMC, "Distribution from Monte Carlo");
  legend->SetFillColor(kWhite);
  legend->Draw();
  TText *figNum;
  if (!zoom) figNum=new TText(.7, .6, "(a)");
  else figNum=new TText(.7, .6, "(b)");
  figNum->SetNDC(kTRUE);
  figNum->Draw("SAME");
  
  TCanvas *c_mee = new TCanvas("c_mee", "c_mee", 500, 500);
  if (!zoom) 
  {
    gPad->SetLogx(1);
    gPad->SetLeftMargin(.14);
    h_mee->GetYaxis()->SetTitleOffset(1.7);
  }
  else 
  {
    gPad->SetLeftMargin(0.16);
    h_mee->GetYaxis()->SetTitleOffset(2.1);
  }
  h_mee->GetYaxis()->CenterTitle();
  h_mee->GetXaxis()->SetTitleOffset(1.1);
  h_mee->GetXaxis()->CenterTitle();
  if (!zoom) h_mee->SetMaximum(0.09);
  else h_mee->SetMaximum(0.0035);
  h_meeMC_recalc->SetMarkerStyle(7);
  h_mee->Draw("C");
  h_meeMC_recalc->Draw("SAME");
  legend=new TLegend(.5, .9, .9, .7);
  legend->AddEntry(h_mee, "Analytical Distribution");
  legend->AddEntry(h_meeMC_recalc, "Distribution from Corrected Monte Carlo");
  legend->SetFillColor(kWhite);
  legend->Draw();
  if (!zoom) figNum=new TText(.7, .6, "(a)");
  else figNum=new TText(.7, .6, "(b)");
  figNum->SetNDC(kTRUE);
  figNum->Draw("SAME");
  
  TCanvas *c_me = new TCanvas("c_me", "c_me", 500, 1000);
  c_me->Divide(1,3);
  c_me->cd(1);
  h_me->Draw();
  c_me->cd(2);
  h_Ee->Draw("box");
  c_me->cd(3);
  h_Ee_me->Draw("box");
  
  TCanvas *c_eMass=new TCanvas("c_eMass", "c_eMass", 500, 500);
  h_me->Draw();
  figNum=new TText(.7, .6, "(a)");
  figNum->SetNDC(kTRUE);
  figNum->Draw("SAME");
  
  TCanvas *c_eMass_E=new TCanvas("c_eMass_E", "c_eMass_E", 500, 500);
  gPad->SetLeftMargin(.18);
  h_Ee_me->GetYaxis()->SetTitleOffset(2.3);
  h_Ee_me->GetYaxis()->CenterTitle();
  h_Ee_me->GetXaxis()->CenterTitle();
  h_Ee_me->Draw("box");
  figNum=new TText(.3, .7, "(b)");
  figNum->SetNDC(kTRUE);
  figNum->Draw("SAME");
  
  
  std::cout<<"At q2=0.00002, calc says "<<h_mee->GetBinContent(h_mee->FindBin(0.00002));
  std::cout<<", MC says "<<h_meeMC->GetBinContent(h_meeMC->FindBin(0.00002))<<std::endl;
  
  std::cout<<"At q2=0.00006, calc says "<<h_mee->GetBinContent(h_mee->FindBin(0.00006));
  std::cout<<", MC says "<<h_meeMC->GetBinContent(h_meeMC->FindBin(0.00006))<<std::endl;
  
  //h_mee->KolmogorovTest(h_meeMC);


} 
