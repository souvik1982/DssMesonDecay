#include <iostream>
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TFile.h"
#include "TPaveStats.h"
#include "TLorentzVector.h"
#include <set>
#include <map>

double pi=3.14159265358979;
double jPsiMass_e_=3.092;
double jPsiMass_e_cut_=0.030;
double jPsiMass_mu_=3.09692; //3.09692;
double jPsiMass_mu_cut_=0.030;
double psi2SMass_pdg_=3.686;
double psi2SMass_center=3.67;
double psi2SMass_range=0.054;
double pi0Mass_pdg_=0.1349766;
double pi0PullMass_cut_=2.5;
double pi0Mass_reco_range=0.018;
double lowerShowerCut_=0.07;
double x925Unf_min_=1.05;
double x925Unf_max_=1.17;
double psi2S_pxCut=0.04;
double psi2S_pyCut=0.04;
double psi2S_pzCut=0.04;
double psi2S_eCenter=-0.005, psi2S_eCut=0.022;
double mee_center=0;
double mee_range=0.009;
double pi0_2_missingMass_max=0.058; // 0.058
double diffPsi_center=psi2SMass_pdg_-jPsiMass_mu_;
double diffPsi_range=0.03;
double diffD0Cut=-0.005;
double dPhiCutLess=0.12;
double missingElectronEnergy_min=0.1;
double missingElectronEnergy_max=0.15;
double cosThetaMax_=0.9;
double missingElectronMomentum_min=0.05;
double missingElectronMomentum_max=0.5;
double electronLowerLimit=0.05;
double photonLowerLimit=0.05;
double electronEnergyThreshold=0.144;

bool data=true;

std::string ftoa (float value)
{
	char buffer[50];
	sprintf(buffer,"%1.2f",value);
	std::string str(buffer);
	return str;
}

bool isGoodPsi2S(TLorentzVector psi2S_vector, TLorentzVector lab4momentum)
{
  return (fabs(psi2S_vector.Px()-lab4momentum.Px())<psi2S_pxCut &&
          fabs(psi2S_vector.Py()-lab4momentum.Py())<psi2S_pyCut &&
          fabs(psi2S_vector.Pz()-lab4momentum.Pz())<psi2S_pzCut); // &&
          // fabs(psi2S_vector.E()-lab4momentum.E())<psi2S_eCut);
}

bool isGoodJPsi(double jPsiMass, int lepton)
{
  if (lepton==0) return fabs((jPsiMass-jPsiMass_mu_)<jPsiMass_mu_cut_);
  if (lepton==1) return fabs((jPsiMass-jPsiMass_e_)<jPsiMass_e_cut_);
}

bool isGoodPi0(double pi0PullMass)
{
 return fabs(pi0PullMass)<pi0PullMass_cut_;
}

bool isGoodLowerShower(double loShower)
{
  return (loShower>lowerShowerCut_);
}

bool isGoodPi0_2_Missing(double pi0MissMassSq)
{
  return (pi0MissMassSq>0 && pi0MissMassSq<pi0_2_missingMass_max);
}

bool isGoodPi0_2(double pi0Mass)
{
  return (fabs(pi0Mass-pi0Mass_pdg_)<pi0Mass_reco_range);
}

bool dD0(float diff)
{
  return ((diff)>diffD0Cut);
}

double deltaPhi(double phi1, double phi2)
{
  double dphi=phi1-phi2;
  if (dphi<-pi) dphi=2*pi+dphi;
  if (dphi>pi) dphi=2*pi-dphi;
  return dphi;
}

bool dPhiCut(float dPhi)
{
  return (dPhi<dPhiCutLess);
}

bool electronEnergyCut(double electronEnergy, double positronEnergy)
{
  return (electronEnergy<electronEnergyThreshold && positronEnergy<electronEnergyThreshold);
}

bool isGoodDiffPsi(double diffPsi)
{
  return fabs(diffPsi-diffPsi_center)<diffPsi_range;
}

bool missingElectronEnergy(double energy)
{
  if (energy>missingElectronEnergy_min && energy<missingElectronEnergy_max) return true;
  else return false;
}

bool missingElectronMomentum(double p)
{
  if (p>missingElectronMomentum_min && p<missingElectronMomentum_max) return true;
  else return false;
}

bool missingElectronAngle(double px, double py, double pz, TH1F *h_e2Missing_cosTheta)
{
  double pt=pow(px*px+py*py, 0.5);
  double costheta=pz/pow(pt*pt+pz*pz, 0.5);
  h_e2Missing_cosTheta->Fill(costheta);
  // std::cout<<"costheta = "<<costheta<<std::endl;
  return (fabs(costheta)<cosThetaMax_);
}

Double_t e_eff_Fit(Double_t *x, Double_t *par)
{
  return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2], 2))
        +par[3]*exp(-0.5*pow((x[0]-par[4])/par[5], 2))
        +par[6]*exp(-0.5*pow((x[0]-par[7])/par[8], 2));
}

Double_t e_ineff_Fit(Double_t *x, Double_t *par)
{
  return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2], 2))
        +par[3]*exp(-0.5*pow((x[0]-par[4])/par[5], 2))
        +par[6]*exp(-0.5*pow((x[0]-par[7])/par[8], 2))
        +par[9]*exp(-0.5*pow((x[0]-par[10])/par[11], 2));
}

Double_t e_eff_Peak(Double_t *x, Double_t *par)
{
  return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2], 2))+par[3]*exp(-0.5*pow((x[0]-par[4])/par[5], 2));
}

Double_t e_eff_Back(Double_t *x, Double_t *par)
{
  return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2], 2));
}

Double_t e_ineff_Peak(Double_t *x, Double_t *par)
{
  return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2], 2))
        +par[3]*exp(-0.5*pow((x[0]-par[4])/par[5], 2));
}

Double_t e_ineff_Back(Double_t *x, Double_t *par)
{
  return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2], 2))
        +par[3]*exp(-0.5*pow((x[0]-par[4])/par[5], 2));
}

int LowEnergyElectronEfficiency_Reco()
{
  int runNumber, eventNumber;
  float lab4momentum_e, lab4momentum_px, lab4momentum_py, lab4momentum_pz;
  float jPsi_e, jPsi_px, jPsi_py, jPsi_pz, jPsi_m;
  float jPsi_e_mc, jPsi_px_mc, jPsi_py_mc, jPsi_pz_mc;
  int jPsi_matched;
  int jPsi_lepton; // 0 if electron, 1 if muon
  float pi0_e, pi0_px, pi0_py, pi0_pz, pi0_RawMass, pi0_PullMass, pi0_hiShower_e, pi0_loShower_e;
  float pi0_e_mc, pi0_px_mc, pi0_py_mc, pi0_pz_mc;
  int pi0_matched;
  float e1_e, e1_px, e1_py, e1_pz, e1_d0;
  float e1_e_mc, e1_px_mc, e1_py_mc, e1_pz_mc; // positron
  int e1_charge;
  int e1_matched;
  float e2_e, e2_px, e2_py, e2_pz, e2_d0;
  float e2_e_mc, e2_px_mc, e2_py_mc, e2_pz_mc; // electron
  float ph_e, ph_px, ph_py, ph_pz, x925Unf;
  float ph_e_mc, ph_px_mc, ph_py_mc, ph_pz_mc;
  int ph_matched;
  
  TFile *file;
  
  if (data) file=new TFile("/nfs/cor/an2/souvik/LowEnergyElectronEfficiency/LowEnergyElectronEfficiency_ReTaggedData_219739_220729.root");
  else 
  {
    // file=new TFile("/nfs/cor/an2/souvik/LowEnergyElectronEfficiency_MC/LowEnergyElectronEfficiencyMC_219739.root");
    file=new TFile("/nfs/cor/an2/souvik/LowEnergyElectronEfficiency_gamma_MC/LowEnergyElectronEfficiency_gamma_MC_219739.root");
  }
  
  TTree* nt1=(TTree*)gDirectory->Get("recoTree");
  
  nt1->SetBranchAddress("kRun_", &(runNumber));
  nt1->SetBranchAddress("kEvent_", &(eventNumber));
  nt1->SetBranchAddress("kLab4Momentum_E_", &(lab4momentum_e));
  nt1->SetBranchAddress("kLab4Momentum_Px_", &(lab4momentum_px));
  nt1->SetBranchAddress("kLab4Momentum_Py_", &(lab4momentum_py));
  nt1->SetBranchAddress("kLab4Momentum_Pz_", &(lab4momentum_pz));
  nt1->SetBranchAddress("kJPsi_E_", &(jPsi_e));
  nt1->SetBranchAddress("kJPsi_Px_", &(jPsi_px));
  nt1->SetBranchAddress("kJPsi_Py_", &(jPsi_py));
  nt1->SetBranchAddress("kJPsi_Pz_", &(jPsi_pz));
  nt1->SetBranchAddress("kJPsi_Matched_", &(jPsi_matched));
  nt1->SetBranchAddress("kJPsiDecayMode_", &(jPsi_lepton));
  nt1->SetBranchAddress("kPi0_1_E_", &(pi0_e));
  nt1->SetBranchAddress("kPi0_1_Px_", &(pi0_px));
  nt1->SetBranchAddress("kPi0_1_Py_", &(pi0_py));
  nt1->SetBranchAddress("kPi0_1_Pz_", &(pi0_pz));
  nt1->SetBranchAddress("kPi0_1_RawMass_", &(pi0_RawMass));
  nt1->SetBranchAddress("kPi0_1_PullMass_", &(pi0_PullMass));
  nt1->SetBranchAddress("kPi0_1_HiShower_E_", &(pi0_hiShower_e));
  nt1->SetBranchAddress("kPi0_1_LoShower_E_", &(pi0_loShower_e));
  nt1->SetBranchAddress("kPi0_1_Matched_", &(pi0_matched));
  nt1->SetBranchAddress("kPi0_2_Electron1_E_", &(e1_e));
  nt1->SetBranchAddress("kPi0_2_Electron1_Px_", &(e1_px));
  nt1->SetBranchAddress("kPi0_2_Electron1_Py_", &(e1_py));
  nt1->SetBranchAddress("kPi0_2_Electron1_Pz_", &(e1_pz));
  nt1->SetBranchAddress("kPi0_2_Electron1_D0_", &(e1_d0));
  nt1->SetBranchAddress("kPi0_2_Electron1_Charge_", &(e1_charge));
  nt1->SetBranchAddress("kPi0_2_Electron1_Matched_", &(e1_matched));
  nt1->SetBranchAddress("kPi0_2_Electron2_E_", &(e2_e));
  nt1->SetBranchAddress("kPi0_2_Electron2_Px_", &(e2_px));
  nt1->SetBranchAddress("kPi0_2_Electron2_Py_", &(e2_py));
  nt1->SetBranchAddress("kPi0_2_Electron2_Pz_", &(e2_pz));
  nt1->SetBranchAddress("kPi0_2_Electron2_D0_", &(e2_d0));
  nt1->SetBranchAddress("kPi0_2_Photon_E_", &(ph_e));
  nt1->SetBranchAddress("kPi0_2_Photon_Px_", &(ph_px));
  nt1->SetBranchAddress("kPi0_2_Photon_Py_", &(ph_py));
  nt1->SetBranchAddress("kPi0_2_Photon_Pz_", &(ph_pz));
  nt1->SetBranchAddress("kx925Unf_", &(x925Unf));
  nt1->SetBranchAddress("kPi0_2_Photon_Matched_", &(ph_matched));
  nt1->SetBranchStatus("*", 1);
  
  TH1F *h_Psi2S_m = new TH1F("h_Psi2S_m", "#psi(2S) Mass; m (GeV); Events / 15 MeV", 100, 3.0, 4.5);
  TH1F *h_Psi2S_px = new TH1F("h_Psi2S_px", "#psi(2S) p_{x}; p_{x} (GeV); Events / 12 MeV", 100, -.2, .2);
  TH1F *h_Psi2S_py = new TH1F("h_Psi2S_py", "#psi(2S) p_{y}; p_{y} (GeV); Events / 12 MeV", 100, -.2, .2);
  TH1F *h_Psi2S_pz = new TH1F("h_Psi2S_pz", "#psi(2S) p_{z}; p_{z} (GeV); Events / 12 MeV", 100, -.2, .2);
  TH1F *h_diffPsi = new TH1F("h_diffPsi", "#psi(2S)-J/#psi Mass; m_{#psi(2S)}-m_{J/#psi} (GeV)", 100, 0.2, 0.9);
  TH1F *h_diffPsi_px = new TH1F("h_diffPsi_px", "px_{#Psi(2S)}-px_{Beam}; #Deltap_{x} (GeV)", 100, -0.2, 0.2);
  TH1F *h_diffPsi_py = new TH1F("h_diffPsi_py", "py_{#Psi(2S)}-py_{Beam}; #Deltap_{y} (GeV)", 100, -0.2, 0.2);
  TH1F *h_diffPsi_pz = new TH1F("h_diffPsi_pz", "pz_{#Psi(2S)}-pz_{Beam}; #Deltap_{z} (GeV)", 100, -0.2, 0.2);
  TH1F *h_diffPsi_e = new TH1F("h_diffPsi_e", "e_{#Psi(2S)}-e_{Beam}", 100, -0.2, 0.2);
  TH1F *h_Psi2S_e = new TH1F("h_Psi2S_e", "#Psi(2S) Energy; E (GeV)", 100, 0., 5.);
  TH1F *h_jPsi_m_e = new TH1F("h_jPsi_m_e", "J/#psi Mass from e^{+}e^{-}; m (GeV)", 50, 3.0492, 3.14692);
  TH1F *h_jPsi_m_m = new TH1F("h_jPsi_m_m", "J/#psi Mass from #mu^{+}#mu^{-}; m (GeV)", 50, 3.04692, 3.14692);
  TH2F *h_jPsi_Reco_Truth_p = new TH2F("h_jPsi_Reco_Truth_p", "h_jPsi_Reco_Truth_p", 100, 0., .5, 100, 0., .5);
  TH2F *h_pi0Mass_RawPull = new TH2F("h_pi0Mass_RawPull", "h_pi0Mass_RawPull", 100, 0.11, 0.16, 100, -3., 3.);
  TH1F *h_pi0Mass_Pull = new TH1F("h_pi0Mass_Pull", "#pi^{0} Pull Mass; #sigma", 50, -3., 3.);
  TH2F *h_pi0_Reco_Truth_p = new TH2F("h_pi0_Reco_Truth_p", "h_pi0_Reco_Truth_p", 100, 0., 0.5, 100, 0., 0.5);
  TH1F *h_pi0_hiShower_E = new TH1F("h_pi0_hiShower_E", "h_pi0_hiShower_E", 100, 0., 0.5);
  TH1F *h_pi0_loShower_E = new TH1F("h_pi0_loShower_E", "h_pi0_loShower_E", 100, 0., 0.5);
  TH2F *h_pi0_RecoMC_hiShower = new TH2F("h_pi0_RecoMC_hiShower", "h_pi0_RecoMC_hiShower", 100, 0, 0.5, 100, -0.2, 0.2);
  TH2F *h_pi0_RecoMC_loShower = new TH2F("h_pi0_RecoMC_loShower", "h_pi0_RecoMC_loShower", 100, 0, 0.5, 100, -0.2, 0.2);
  TH1F *h_pi0_2_missing = new TH1F("h_pi0_2_missing", "Missing Mass Squared of 2nd #pi^{0}; m^{2} (GeV^{2})", 50, -0.2, 0.2);
  TH1F *h_pi0_2_missing_p = new TH1F("h_pi0_2_missing_p", "h_pi0_2_missing_p", 100, 0., 0.5);
  TH1F *h_pi0_2_Reco_m = new TH1F("h_pi0_2_Reco_m", "Invariant Mass of Second #pi^{0}; m (GeV)", 100, 0.05, 0.4);
  TH2F *h_pi0_2_Reco_m_e1diff = new TH2F("h_pi0_2_Reco_m_e1diff", "h_pi0_2_Reco_m_e1diff", 100, -0.1, 0.1, 100, 0.05, 0.2);
  TH2F *h_pi0_2_Reco_m_phdiff = new TH2F("h_pi0_2_Reco_m_phdiff", "h_pi0_2_Reco_m_phdiff", 100, -0.1, 0.1, 100, 0.05, 0.2);
  TH2F *h_e1_Reco_Truth_p = new TH2F("h_e1_Reco_Truth_p", "h_e1_Reco_Truth_p", 100, 0., 0.5, 100, 0., 0.5);
  TH2F *h_ph_Reco_Truth_p = new TH2F("h_ph_Reco_Truth_p", "h_ph_Reco_Truth_p", 100, 0., 0.5, 100, 0., 0.5);
  TH1F *h_x925 = new TH1F("h_x925", "h_x925", 100, 0.5, 1.5);
  TH2F *h_ph_RecoTruth_x925 = new TH2F("h_ph_RecoTruth_x925", "h_ph_RecoTruth_x925", 100, 0.5, 1.5, 100, -0.5, 0.5);
  TH1D *h_diffD0 = new TH1D("h_diffD0", "#Deltad_{0}; #Deltad_{0} m", 50, -0.01, 0.01);
  TH1D *h_dPhi = new TH1D("h_dPhi", "#Delta#phi_{0}; #Delta#phi_{0}", 50, -pi/2, pi/2);
  TH1D *h_n_dPhi = new TH1D("h_n_dPhi", "Conversion Type #Delta#phi_{0}", 50, -pi/2, pi/2);
  
  TH1D *h_e_Reco_e = new TH1D("h_e_Reco_e", "Energy of the Electron; E (GeV)", 100, 0., 0.3);
  
  unsigned int yield=0;
  unsigned int sidebandYield=0;
  unsigned int oldRun=0, oldEvent=0;
  unsigned int oldRunSideband=0, oldEventSideband=0;
  
  unsigned int nEvents=nt1->GetEntries();
  for (int i=0; i<nEvents; ++i)
  {
    nt1->GetEvent(i);
    
    TLorentzVector lab4momentum(lab4momentum_px, lab4momentum_py, lab4momentum_pz, lab4momentum_e);
    
    TLorentzVector jPsi_vector(jPsi_px, jPsi_py, jPsi_pz, jPsi_e);
    if (jPsi_lepton==0) h_jPsi_m_e->Fill(jPsi_vector.M());
    else if (jPsi_lepton==1) h_jPsi_m_m->Fill(jPsi_vector.M());
    // if (jPsi_lepton==1)
    {
    if (isGoodJPsi(jPsi_vector.M(), jPsi_lepton))
    {
      TLorentzVector pi0_vector(pi0_px, pi0_py, pi0_pz, pi0_e);
      h_pi0Mass_RawPull->Fill(pi0_RawMass, pi0_PullMass);
      h_pi0Mass_Pull->Fill(pi0_PullMass);
      if (isGoodPi0(pi0_PullMass))
      {
        TLorentzVector e1_vector(e1_px, e1_py, e1_pz, e1_e);
        TLorentzVector e2_vector(e2_px, e2_py, e2_pz, e2_e);
        TLorentzVector ph_vector(ph_px, ph_py, ph_pz, ph_e);
        TLorentzVector pi0_2_vector=e1_vector+e2_vector+ph_vector;
        h_pi0_2_Reco_m->Fill(pi0_2_vector.M());
        h_e_Reco_e->Fill(e1_e);
        h_e_Reco_e->Fill(e2_e);
        if (isGoodPi0_2(pi0_2_vector.M()) && electronEnergyCut(e1_e, e2_e))
        {
          double diffD0;
          if (e1_charge>0) diffD0=e2_d0-e1_d0;
          else diffD0=e1_d0-e2_d0;
          double phi_e1=atan2(e1_py, e1_px);
          double phi_e2=atan2(e2_py, e2_px);
          double dPhi;
          if (e1_charge>0) dPhi=deltaPhi(phi_e2, phi_e1);
          else dPhi=deltaPhi(phi_e1, phi_e2);
          
          TLorentzVector psi2S_vector=jPsi_vector+pi0_vector+pi0_2_vector;
          h_diffPsi_px->Fill(psi2S_vector.Px()-lab4momentum.Px());
          h_diffPsi_py->Fill(psi2S_vector.Py()-lab4momentum.Py());
          h_diffPsi_pz->Fill(psi2S_vector.Pz()-lab4momentum.Pz());
          h_diffPsi_e->Fill(psi2S_vector.E()-lab4momentum.E());
          if (isGoodPsi2S(psi2S_vector, lab4momentum))
          {
            float diffPsi=psi2S_vector.M()-jPsi_vector.M();
            h_diffPsi->Fill(diffPsi);
            if (isGoodDiffPsi(diffPsi))
            {
              // h_e_Reco_e->Fill(e1_e);
              // h_e_Reco_e->Fill(e2_e);
              h_diffD0->Fill(diffD0);
              if (dD0(diffD0))
              {
                h_dPhi->Fill(dPhi);
                if (dPhiCut(dPhi))
                {
                  if (oldRun!=runNumber || oldEvent!=eventNumber)
                  {
                    ++yield;
                    // h_e_Reco_e->Fill(e1_e);
                    // h_e_Reco_e->Fill(e2_e);
                  }
                  /* {
                    std::cout<<"Run number = "<<runNumber<<", Event number = "<<eventNumber<<std::endl;
                    std::cout<<"e1, e = "<<e1_e<<", e2, e = "<<e2_e<<std::endl;
                  } */
                  oldRun=runNumber;
                  oldEvent=eventNumber;
                }
              }
              else
              {
                h_n_dPhi->Fill(dPhi);
                if (!dPhiCut(dPhi))
                {
                  if (oldRunSideband!=runNumber || oldEventSideband!=eventNumber)
                  {
                    ++sidebandYield;
                  }
                  oldRunSideband=runNumber;
                  oldEventSideband=eventNumber;
                }
              }
            }
          }
        }
      }
    }
    }
  } // Iterate over events
  
  std::cout<<"Yield = "<<yield<<std::endl;
  std::cout<<"Sideband yield = "<<sidebandYield<<std::endl;
  
  TLine *line;
  float xmin, xmax, ymax;
  gROOT->SetStyle("Plain");
  
  TCanvas *jPsi_mass = new TCanvas("jPsi_mass", "jPsi_mass", 500, 1000);
  jPsi_mass->Divide(1,2);
  jPsi_mass->cd(1);
  h_jPsi_m_e->Draw();  
  xmin=jPsiMass_e_-jPsiMass_e_cut_;
  xmax=jPsiMass_e_+jPsiMass_e_cut_;
  ymax=(h_jPsi_m_e->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  jPsi_mass->cd(2);
  h_jPsi_m_m->Draw();
  xmin=jPsiMass_mu_-jPsiMass_mu_cut_;
  xmax=jPsiMass_mu_+jPsiMass_mu_cut_;
  ymax=(h_jPsi_m_m->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  
  TCanvas *pi01_mass = new TCanvas("pi01_mass", "pi01_mass", 500, 1000);
  pi01_mass->Divide(1,2);
  pi01_mass->cd(1);
  h_pi0Mass_RawPull->Draw("box");
  pi01_mass->cd(2);
  h_pi0Mass_Pull->Draw();
  line = new TLine(pi0PullMass_cut_, h_pi0Mass_Pull->GetMinimum(), pi0PullMass_cut_, h_pi0Mass_Pull->GetMaximum()*0.90); line->Draw();
  line = new TLine(-pi0PullMass_cut_, h_pi0Mass_Pull->GetMinimum(), -pi0PullMass_cut_, h_pi0Mass_Pull->GetMaximum()*0.90); line->Draw();
  
  TCanvas *pi02_mass = new TCanvas("pi02_mass", "pi02_mass", 500, 500);
  h_pi0_2_Reco_m->Draw();
  line = new TLine(pi0Mass_pdg_-pi0Mass_reco_range, 0, pi0Mass_pdg_-pi0Mass_reco_range, h_pi0_2_Reco_m->GetMaximum()*0.90); line->Draw();
  line = new TLine(pi0Mass_pdg_+pi0Mass_reco_range, 0, pi0Mass_pdg_+pi0Mass_reco_range, h_pi0_2_Reco_m->GetMaximum()*0.90); line->Draw();
  
  TCanvas *c_diffD0 = new TCanvas("c_diffD0", "c_diffD0", 500, 500);
  h_diffD0->Draw();
  line = new TLine(diffD0Cut, 0, diffD0Cut, h_diffD0->GetMaximum()); line->Draw();
  
  TCanvas *c_dPhi = new TCanvas("c_dPhi", "c_dPhi", 500, 500);
  h_dPhi->Draw();
  line = new TLine(dPhiCutLess, 0, dPhiCutLess, h_dPhi->GetMaximum()); line->Draw();
  
  TCanvas *c_n_dPhi = new TCanvas("c_n_dPhi", "c_n_dPhi", 500, 500);
  h_n_dPhi->Draw();
  line = new TLine(dPhiCutLess, 0, dPhiCutLess, h_n_dPhi->GetMaximum()); line->Draw();
  
  TCanvas *c_RelativePsi = new TCanvas("c_RelativePsi", "c_RelativePsi", 500, 1500);
  c_RelativePsi->Divide(1,3);
  c_RelativePsi->cd(1);
  h_diffPsi_px->Draw();
  xmin=-psi2S_pxCut;
  xmax=psi2S_pxCut;
  ymax=(h_diffPsi_px->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  c_RelativePsi->cd(2);
  h_diffPsi_py->Draw();
  xmin=-psi2S_pyCut;
  xmax=psi2S_pyCut;
  ymax=(h_diffPsi_py->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  c_RelativePsi->cd(3);
  h_diffPsi_pz->Draw();
  xmin=-psi2S_pzCut;
  xmax=psi2S_pzCut;
  ymax=(h_diffPsi_pz->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  
  TCanvas *c_RelativePsi_e = new TCanvas("c_RelativePsi_e");
  h_diffPsi_e->Draw();
  xmin=-psi2S_eCut;
  xmax=psi2S_eCut;
  ymax=(h_diffPsi_e->GetMaximum())*0.75;
  //line = new TLine(xmin,ymax,xmin,0); line->Draw();
  //line = new TLine(xmax,ymax,xmax,0); line->Draw();
  
  TCanvas *c_diffPsi = new TCanvas("c_diffPsi", "c_diffPsi", 500, 500);
  c_diffPsi->Divide(1,1);
  c_diffPsi->cd(1);
  h_diffPsi->Draw();
  line = new TLine(diffPsi_center-diffPsi_range, 0, diffPsi_center-diffPsi_range, h_diffPsi->GetMaximum()*0.90); line->Draw();
  line = new TLine(diffPsi_center+diffPsi_range, 0, diffPsi_center+diffPsi_range, h_diffPsi->GetMaximum()*0.90); line->Draw();
  
  TCanvas *c_e_Reco_e = new TCanvas("c_e_Reco_e", "c_e_Reco_e", 500, 500);
  h_e_Reco_e->Draw();
  line = new TLine(electronEnergyThreshold, 0, electronEnergyThreshold, h_e_Reco_e->GetMaximum()*0.90); line->Draw();
  
  TFile *tfile = new TFile("eEnergy.root", "RECREATE");
  h_e_Reco_e->Write();
  tfile->Close();
  
  return 0;
}
