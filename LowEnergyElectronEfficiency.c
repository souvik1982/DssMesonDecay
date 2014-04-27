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

double pi_=3.14159265358979;
double jPsiMass_e_=3.092;
double jPsiMass_e_cut_=0.015;
double jPsiMass_mu_=3.09692; //3.09692;
double jPsiMass_mu_cut_=0.015;
double psi2SMass_pdg_=3.686;
double psi2SMass_center=3.67;
double psi2SMass_range=0.054;
double pi0Mass_pdg_=0.1349766;
double pi0PullMass_cut_=1.5;
double pi0Mass_reco_range=0.012;
double lowerShowerCut_=0.07;
double x925Unf_min_=1.05;
double x925Unf_max_=1.17;
double psi2S_pxCut=0.022;
double psi2S_pyCut=0.022;
double psi2S_pzCut=0.022;
double psi2S_eCut=0.022;
double mee_center=0;
double mee_range=0.009;
double pi0_2_missingMass_max=0.058; // 0.058
double diffPsi_center=psi2SMass_pdg_-jPsiMass_mu_;
double diffPsi_range=0.03;
double missingElectronEnergy_min=0.1;
double missingElectronEnergy_max=0.15;
double cosThetaMax_=0.9;
double missingElectronMomentum_min=0.05;
double missingElectronMomentum_max=0.5;
double electronLowerLimit=0.05;
double photonLowerLimit=0.05;

bool data=false;

std::string ftoa (float value)
{
	char buffer[50];
	sprintf(buffer,"%1.2f",value);
	std::string str(buffer);
	return str;
}

bool isGoodPsi2S(double px, double py, double pz, double e)
{
  return (fabs(px)<psi2S_pxCut &&
          fabs(py)<psi2S_pyCut &&
          fabs(pz)<psi2S_pzCut);
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

int LowEnergyElectronEfficiency()
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
  float e1_e, e1_px, e1_py, e1_pz;
  float e1_e_mc, e1_px_mc, e1_py_mc, e1_pz_mc; // positron
  int e1_charge;
  int e1_matched;
  int e2_n;
  float e2_e[100], e2_px[100], e2_py[100], e2_pz[100];
  float e2_e_mc, e2_px_mc, e2_py_mc, e2_pz_mc; // electron
  float ph_e, ph_px, ph_py, ph_pz, x925Unf;
  float ph_e_mc, ph_px_mc, ph_py_mc, ph_pz_mc;
  int ph_matched;
  
  TFile *file;
  
  if (data) file=new TFile("/nfs/cor/an2/souvik/LowEnergyElectronEfficiency/LowEnergyElectronEfficiency_ReTaggedData_219739_220729.root");
  else file=new TFile("/nfs/cor/an2/souvik/LowEnergyElectronEfficiency_MC/LowEnergyElectronEfficiencyMC_219739.root");
  
  TTree* nt1=(TTree*)gDirectory->Get("LowEnergyElectronEfficiency");
  
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
  nt1->SetBranchAddress("kJPsi_E_MC_", &(jPsi_e_mc));
  nt1->SetBranchAddress("kJPsi_Px_MC_", &(jPsi_px_mc));
  nt1->SetBranchAddress("kJPsi_Py_MC_", &(jPsi_py_mc));
  nt1->SetBranchAddress("kJPsi_Pz_MC_", &(jPsi_pz_mc));
  nt1->SetBranchAddress("kPi0_1_E_", &(pi0_e));
  nt1->SetBranchAddress("kPi0_1_Px_", &(pi0_px));
  nt1->SetBranchAddress("kPi0_1_Py_", &(pi0_py));
  nt1->SetBranchAddress("kPi0_1_Pz_", &(pi0_pz));
  nt1->SetBranchAddress("kPi0_1_RawMass_", &(pi0_RawMass));
  nt1->SetBranchAddress("kPi0_1_PullMass_", &(pi0_PullMass));
  nt1->SetBranchAddress("kPi0_1_HiShower_E_", &(pi0_hiShower_e));
  nt1->SetBranchAddress("kPi0_1_LoShower_E_", &(pi0_loShower_e));
  nt1->SetBranchAddress("kFirstPi0_E_MC_", &(pi0_e_mc));
  nt1->SetBranchAddress("kFirstPi0_Px_MC_", &(pi0_px_mc));
  nt1->SetBranchAddress("kFirstPi0_Py_MC_", &(pi0_py_mc));
  nt1->SetBranchAddress("kFirstPi0_Pz_MC_", &(pi0_pz_mc));
  nt1->SetBranchAddress("kPi0_1_Matched_", &(pi0_matched));
  nt1->SetBranchAddress("kPi0_2_Electron1_E_", &(e1_e));
  nt1->SetBranchAddress("kPi0_2_Electron1_Px_", &(e1_px));
  nt1->SetBranchAddress("kPi0_2_Electron1_Py_", &(e1_py));
  nt1->SetBranchAddress("kPi0_2_Electron1_Pz_", &(e1_pz));
  nt1->SetBranchAddress("kPi0_2_Electron1_Charge_", &(e1_charge));
  nt1->SetBranchAddress("kPi0_2_Electron1_Matched_", &(e1_matched));
  nt1->SetBranchAddress("kSecondPi0Positron_E_MC_", &(e1_e_mc));
  nt1->SetBranchAddress("kSecondPi0Positron_Px_MC_", &(e1_px_mc));
  nt1->SetBranchAddress("kSecondPi0Positron_Py_MC_", &(e1_py_mc));
  nt1->SetBranchAddress("kSecondPi0Positron_Pz_MC_", &(e1_pz_mc));
  nt1->SetBranchAddress("kPi0_2_Electron2_N_", &(e2_n));
  nt1->SetBranchAddress("kPi0_2_Electron2_E_", &(e2_e));
  nt1->SetBranchAddress("kPi0_2_Electron2_Px_", &(e2_px));
  nt1->SetBranchAddress("kPi0_2_Electron2_Py_", &(e2_py));
  nt1->SetBranchAddress("kPi0_2_Electron2_Pz_", &(e2_pz));
  nt1->SetBranchAddress("kSecondPi0Electron_E_MC_", &(e2_e_mc));
  nt1->SetBranchAddress("kSecondPi0Electron_Px_MC_", &(e2_px_mc));
  nt1->SetBranchAddress("kSecondPi0Electron_Py_MC_", &(e2_py_mc));
  nt1->SetBranchAddress("kSecondPi0Electron_Pz_MC_", &(e2_pz_mc));
  nt1->SetBranchAddress("kPi0_2_Photon_E_", &(ph_e));
  nt1->SetBranchAddress("kPi0_2_Photon_Px_", &(ph_px));
  nt1->SetBranchAddress("kPi0_2_Photon_Py_", &(ph_py));
  nt1->SetBranchAddress("kPi0_2_Photon_Pz_", &(ph_pz));
  nt1->SetBranchAddress("kx925Unf_", &(x925Unf));
  nt1->SetBranchAddress("kPi0_2_Photon_Matched_", &(ph_matched));
  nt1->SetBranchAddress("kSecondPi0Photon_E_MC_", &(ph_e_mc));
  nt1->SetBranchAddress("kSecondPi0Photon_Px_MC_", &(ph_px_mc));
  nt1->SetBranchAddress("kSecondPi0Photon_Py_MC_", &(ph_py_mc));
  nt1->SetBranchAddress("kSecondPi0Photon_Pz_MC_", &(ph_pz_mc));
  nt1->SetBranchStatus("*", 1);
  
  TH1F *h_Psi2S_m = new TH1F("h_Psi2S_m", "#Psi(2S) Mass; m (GeV); Events / 15 MeV", 100, 3.0, 4.5);
  TH1F *h_Psi2S_px = new TH1F("h_Psi2S_px", "#Psi(2S) p_{x}; p_{x} (GeV); Events / 12 MeV", 100, -.2, .2);
  TH1F *h_Psi2S_py = new TH1F("h_Psi2S_py", "#Psi(2S) p_{y}; p_{y} (GeV); Events / 12 MeV", 100, -.2, .2);
  TH1F *h_Psi2S_pz = new TH1F("h_Psi2S_pz", "#Psi(2S) p_{z}; p_{z} (GeV); Events / 12 MeV", 100, -.2, .2);
  TH1F *h_diffPsi = new TH1F("h_diffPsi", "#Psi(2S)-J/#psi Mass", 100, 0.2, 0.9);
  TH1F *h_diffPsi_px = new TH1F("h_diffPsi_px", "px_{#Psi(2S)}-px_{Beam}", 100, -0.2, 0.2);
  TH1F *h_diffPsi_py = new TH1F("h_diffPsi_py", "py_{#Psi(2S)}-py_{Beam}", 100, -0.2, 0.2);
  TH1F *h_diffPsi_pz = new TH1F("h_diffPsi_pz", "pz_{#Psi(2S)}-pz_{Beam}", 100, -0.2, 0.2);
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
  TH1F *h_pi0_2_Reco_m = new TH1F("h_pi0_2_Reco_m", "h_pi0_2_Reco_m", 100, 0.05, 0.4);
  TH2F *h_pi0_2_Reco_m_e1diff = new TH2F("h_pi0_2_Reco_m_e1diff", "h_pi0_2_Reco_m_e1diff", 100, -0.1, 0.1, 100, 0.05, 0.2);
  TH2F *h_pi0_2_Reco_m_phdiff = new TH2F("h_pi0_2_Reco_m_phdiff", "h_pi0_2_Reco_m_phdiff", 100, -0.1, 0.1, 100, 0.05, 0.2);
  TH2F *h_e1_Reco_Truth_p = new TH2F("h_e1_Reco_Truth_p", "h_e1_Reco_Truth_p", 100, 0., 0.5, 100, 0., 0.5);
  TH2F *h_ph_Reco_Truth_p = new TH2F("h_ph_Reco_Truth_p", "h_ph_Reco_Truth_p", 100, 0., 0.5, 100, 0., 0.5);
  TH1F *h_x925 = new TH1F("h_x925", "h_x925", 100, 0.5, 1.5);
  TH2F *h_ph_RecoTruth_x925 = new TH2F("h_ph_RecoTruth_x925", "h_ph_RecoTruth_x925", 100, 0.5, 1.5, 100, -0.5, 0.5);
  
  TH1F *h_e_eff_Mass2 = new TH1F("h_e_eff_Mass2", "Efficiency Plot; m_{e}^{2} (GeV^{2})", 50, -0.1, 0.05);
  TH1F *h_e_ineff_Mass2 = new TH1F("h_e_ineff_Mass2", "Inefficiency Plot; m_{e}^{2} (GeV^{2})", 50, -0.1, 0.05);
  std::string title_eff, title_ineff;
  if (data)
  {
    title_eff="Efficiency Plot for Data with ";  title_ineff="Inefficiency Plot for Data with ";
    title_eff+=ftoa(missingElectronEnergy_min);  title_ineff+=ftoa(missingElectronEnergy_min);
    title_eff+=" < E < ";                        title_ineff+=" < E < ";
    title_eff+=ftoa(missingElectronEnergy_max);  title_ineff+=ftoa(missingElectronEnergy_max);
    title_eff+=" GeV";                           title_ineff+=" GeV";
  }
  else
  {
    title_eff="Efficiency Plot for MC with ";    title_ineff="Inefficiency Plot for MC with ";
    title_eff+=ftoa(missingElectronEnergy_min);  title_ineff+=ftoa(missingElectronEnergy_min);
    title_eff+=" < E < ";                        title_ineff+=" < E < ";
    title_eff+=ftoa(missingElectronEnergy_max);  title_ineff+=ftoa(missingElectronEnergy_max);
    title_eff+=" GeV";                           title_ineff+=" GeV";
  }
  h_e_eff_Mass2->SetTitle(title_eff.c_str());
  h_e_ineff_Mass2->SetTitle(title_ineff.c_str());
  
  TH2F *h_e_eff_E_p = new TH2F("h_e_eff_E_p", "h_e_eff_E_p", 100, 0., 0.03, 100, 0., 0.03);
  TH2F *h_e_ineff_E_p = new TH2F("h_e_ineff_E_p", "h_e_ineff_E_p", 100, 0., 0.03, 100, 0., 0.03);
  TH1F *h_eff_Angle = new TH1F("h_eff_Angle", "h_eff_Angle", 100, 0., pi_);
  TH1F *h_ineff_Angle = new TH1F("h_ineff_Angle", "h_ineff_Angle", 100, 0., pi_);
  TH1F *h_e2Missing_p = new TH1F("h_e2Missing_p", "h_e2Missing_p", 100, 0., 0.5);
  TH1F *h_e2Missing_e = new TH1F("h_e2Missing_e", "h_e2Missing_e", 100, 0., 0.5);
  TH1F *h_e2Missing_cosTheta = new TH1F("h_e2Missing_cosTheta", "h_e2Missing_cosTheta", 100, -1, 1);
  TH1I *h_nMissingElectrons = new TH1I("h_nMissingElectrons", "h_nMissingElectrons", 6, 0, 5);
  TH1F *h_e2Candidate_p = new TH1F("h_e2Candidate_p", "h_e2Candidate_p", 100, 0., 0.5);
  TH1F *h_e2Candidate_e = new TH1F("h_e2Candidate_e", "h_e2Candidate_e", 100, 0., 0.5);
  TH1F *h_e2Missing_p_MC = new TH1F("h_e2Missing_p_MC", "h_e2Missing_p_MC", 100, 0., 0.5);
  TH1F *h_e2_p_MC = new TH1F("h_e2_p_MC", "h_e2_p_MC", 100, 0., 0.5);
  TH2F *h_e2_RecoMissing_MC_p = new TH2F("h_e2_RecoMissing_MC_p", "h_e2_RecoMissing_MC_p", 100, 0., 0.5, 100, 0., .5);
  TH2F *h_e2_RecoMissing_MCMissing_p = new TH2F("h_e2_RecoMissing_MCMissing_p", "h_e2_RecoMissing_MCMissing_p", 100, 0., 0.5, 100, 0., 0.5);
  TH2F *h_e2_Missing_MC_p = new TH2F("h_e2_Missing_MC_p", "h_e2_Missing_MC_p", 100, 0., 0.5, 100, 0., 0.5);
  TH2F *h_pi0_2_missing_p_e2_missing_m = new TH2F("h_pi0_2_missing_p_e2_missing_m", "h_pi0_2_missing_p_e2_missing_m", 50, -0.1, 0.1, 100, 0., 0.5);
  TH2F *h_diffPsi_diffpi0 = new TH2F("h_diffPsi_diffpi0", "h_diffPsi_diffpi0", 100, -0.02, 0.02, 100, 0.2, 0.9);
  TH2F *h_diffPsi_pi0_2_Reco_m = new TH2F("h_diffPsi_pi0_2_Reco_m", "h_diffPsi_pi0_2_Reco_m", 100, 0.05, 0.4, 100, 0.2, 0.9);
  TH2F *h_diffPsi_hiShower = new TH2F("h_diffPsi_hiShower", "h_diffPsi_hiShower", 100, 0., 0.5, 100, 0.2, 0.9);
  TH2F *h_diffPsi_loShower = new TH2F("h_diffPsi_loShower", "h_diffPsi_loShower", 100, 0., 0.5, 100, 0.2, 0.9);
  TH1F *h_e2_MissingMass2_Resolution = new TH1F("h_e2_MissingMass2_Resolution", "h_e2_MissingMass2_Resolution", 200, -0.1, 0.1);
  TH2F *h_e2_MissingMass2_MissingP = new TH2F("h_e2_MissingMass2_MissingP", "h_e2_MissingMass2_MissingP", 100, 0., 0.5, 100, -0.1, 0.1);
  TH2F *h_e2_MissingMass2_MissingE = new TH2F("h_e2_MissingMass2_MissingE", "h_e2_MissingMass2_MissingE", 100, 0., 0.3, 100, -0.2, 0.05);
  TH2F *h_Missingpi0_pi0 = new TH2F("h_Missingpi0_pi0", "h_Missingpi0_pi0", 100, 0., 0.2, 100, 0., 0.2);
  
  
  unsigned int nEvents=nt1->GetEntries();
  for (int i=0; i<nEvents; ++i)
  {
    nt1->GetEvent(i);
    
    TLorentzVector psi2S_vector_rest(lab4momentum_px, lab4momentum_py, lab4momentum_pz, lab4momentum_e);
    TLorentzVector psi2S_vector_rest_MC(lab4momentum_px, lab4momentum_py, lab4momentum_pz, lab4momentum_e);
    
    // if (jPsi_lepton==1)
    {
      // if (jPsi_matched==1)
      {
        TLorentzVector jPsi_vector(jPsi_px, jPsi_py, jPsi_pz, jPsi_e);
        TLorentzVector jPsi_vector_MC(jPsi_px_mc, jPsi_py_mc, jPsi_pz_mc, jPsi_e_mc);
        if (jPsi_lepton==0) h_jPsi_m_e->Fill(jPsi_vector.M());
        else if (jPsi_lepton==1) h_jPsi_m_m->Fill(jPsi_vector.M());
 
        if (isGoodJPsi(jPsi_vector.M(), jPsi_lepton))
        {
          // if (pi0_matched==1)
          {
            TLorentzVector pi0_vector(pi0_px, pi0_py, pi0_pz, pi0_e);
            TLorentzVector pi0_vector_MC(pi0_px_mc, pi0_py_mc, pi0_pz_mc, pi0_e_mc);
            h_pi0Mass_RawPull->Fill(pi0_RawMass, pi0_PullMass);
            h_pi0Mass_Pull->Fill(pi0_PullMass);
 
            if (isGoodPi0(pi0_PullMass) && isGoodLowerShower(pi0_loShower_e))
            {
              TLorentzVector pi0_2_missing=psi2S_vector_rest-jPsi_vector-pi0_vector;
              h_pi0_2_missing->Fill(pi0_2_missing.M2());
              
              // h_Missingpi0_pi0->Fill(pi0_vector.P()*pi0_vector.P(), pi0_2_missing.P()*pi0_2_missing.P());
 
              if (isGoodPi0_2_Missing(pi0_2_missing.M2()))
              {
              
                h_Missingpi0_pi0->Fill(pi0_vector.P()*pi0_vector.P(), pi0_2_missing.P()*pi0_2_missing.P());
                
                // if (e1_matched==1 && ph_matched==1) // e1_matched==1 && 
                {
                  TLorentzVector e1_vector(e1_px, e1_py, e1_pz, e1_e);
                  TLorentzVector ph_vector(ph_px, ph_py, ph_pz, ph_e);
                  TLorentzVector e2_vector_exp=psi2S_vector_rest-jPsi_vector-pi0_vector-e1_vector-ph_vector;
                  TLorentzVector e1_vector_MC;
                  if (e1_charge==-1) // electron
                  {
                    e1_vector_MC=TLorentzVector(e2_px_mc, e2_py_mc, e2_pz_mc, e2_e_mc);
                  }
                  else
                  {
                    e1_vector_MC=TLorentzVector(e1_px_mc, e1_py_mc, e1_pz_mc, e1_e_mc);
                  }
                  TLorentzVector ph_vector_MC(ph_px_mc, ph_py_mc, ph_pz_mc, ph_e_mc);
                  TLorentzVector e2_vector_exp_MC=psi2S_vector_rest_MC-jPsi_vector_MC-pi0_vector_MC-e1_vector_MC-ph_vector_MC;
                  
                  h_x925->Fill(x925Unf);
                  h_ph_RecoTruth_x925->Fill(x925Unf, ph_vector.P()-ph_vector_MC.P());
                  // if (fabs(ph_vector.P()-ph_vector_MC.P())<0.02 && fabs(e1_vector.P()-e1_vector_MC.P())<0.02)
                  if (ph_vector.P()>photonLowerLimit && e1_vector.P()>electronLowerLimit && x925Unf>x925Unf_min_ && x925Unf<x925Unf_max_)
                  {
                    TLorentzVector e2_vector_MC;
                    if (e1_charge==-1) // electron, so take the positron
                    {
                      e2_vector_MC=*(new TLorentzVector(e1_px_mc, e1_py_mc, e1_pz_mc, e1_e_mc));
                    }
                    else               // positron, so take the electron
                    {
                      e2_vector_MC=*(new TLorentzVector(e2_px_mc, e2_py_mc, e2_pz_mc, e2_e_mc));
                    }
                    
                    h_e2Missing_p->Fill(e2_vector_exp.P());
                    h_e2Missing_e->Fill(e2_vector_exp.E());
                    // if (missingElectronEnergy(e2_vector_exp.P()) && missingElectronAngle(e2_vector_exp.Px(), e2_vector_exp.Py(), e2_vector_exp.Pz(), h_e2Missing_cosTheta))
                    if (missingElectronEnergy(e2_vector_exp.E()) && missingElectronAngle(e2_vector_exp.Px(), e2_vector_exp.Py(), e2_vector_exp.Pz(), h_e2Missing_cosTheta))
                    {
                      
                      h_jPsi_Reco_Truth_p->Fill(jPsi_vector_MC.P(), jPsi_vector.P());
                      h_pi0_Reco_Truth_p->Fill(pi0_vector_MC.P(), pi0_vector.P());
                      h_e1_Reco_Truth_p->Fill(e1_vector_MC.P(), e1_vector.P());
                      h_ph_Reco_Truth_p->Fill(ph_vector_MC.P(), ph_vector.P());
                      h_e2_RecoMissing_MC_p->Fill(e2_vector_MC.P(), e2_vector_exp.P());
                      h_e2_Missing_MC_p->Fill(e2_vector_MC.P(), e2_vector_exp_MC.P());
                      
                      
                      h_e2_MissingMass2_Resolution->Fill(e2_vector_exp.M2()-e2_vector_exp_MC.M2());
                      h_e2_MissingMass2_MissingP->Fill(e2_vector_exp.P(), e2_vector_exp.M2()-e2_vector_exp_MC.M2());
                      h_e2_MissingMass2_MissingE->Fill(e2_vector_exp.E(), e2_vector_exp.M2()-e2_vector_exp_MC.M2());
                      
                      h_e2_p_MC->Fill(e2_vector_MC.P());
                      
                      h_pi0_2_missing_p->Fill(pi0_2_missing.P());
                      
                      h_pi0_hiShower_E->Fill(pi0_hiShower_e);
                      h_pi0_loShower_E->Fill(pi0_loShower_e);
                      h_pi0_RecoMC_hiShower->Fill(pi0_hiShower_e, pi0_vector.P()-pi0_vector_MC.P());
                      h_pi0_RecoMC_loShower->Fill(pi0_loShower_e, pi0_vector.P()-pi0_vector_MC.P());
                      
                      bool efficientElectron=false;
                      h_nMissingElectrons->Fill(e2_n);
                      for (unsigned int i_lastElectron=0; i_lastElectron<e2_n; ++i_lastElectron)
                      {
                        TLorentzVector e2_vector(e2_px[i_lastElectron], e2_py[i_lastElectron], e2_pz[i_lastElectron], e2_e[i_lastElectron]);
                        TLorentzVector pi0_2_vector=e1_vector+e2_vector+ph_vector;
                        TLorentzVector psi2S_vector=jPsi_vector+pi0_vector+pi0_2_vector;
                        float diffPsi=psi2S_vector.M()-jPsi_vector.M();
                        h_diffPsi->Fill(diffPsi);
                        h_e2Candidate_p->Fill(e2_vector.P());
                        h_e2Candidate_e->Fill(e2_vector.E());
                        h_diffPsi_diffpi0->Fill(pi0_vector.P()-pi0_vector_MC.P(), diffPsi);
                        h_pi0_2_Reco_m_e1diff->Fill(e1_vector.P()-e1_vector_MC.P(), pi0_2_vector.M());
                        h_pi0_2_Reco_m_phdiff->Fill(ph_vector.P()-ph_vector_MC.P(), pi0_2_vector.M());
                        h_diffPsi_hiShower->Fill(pi0_hiShower_e, diffPsi);
                        h_diffPsi_loShower->Fill(pi0_loShower_e, diffPsi);
                        
                        //if (isGoodPi0_2(pi0_2_vector.M()))
                        {
                          h_diffPsi_pi0_2_Reco_m->Fill(pi0_2_vector.M(), diffPsi);
                          if (isGoodDiffPsi(diffPsi))
                          {
                            h_pi0_2_Reco_m->Fill(pi0_2_vector.M());
                            
                            h_Psi2S_m->Fill(psi2S_vector.M());
                            h_Psi2S_px->Fill(psi2S_vector.Px());
                            h_Psi2S_py->Fill(psi2S_vector.Py());
                            h_Psi2S_pz->Fill(psi2S_vector.Pz());
                            h_Psi2S_e->Fill(psi2S_vector.E());
 
                            double diffPx=psi2S_vector.Px()-lab4momentum_px;
                            double diffPy=psi2S_vector.Py()-lab4momentum_py;
                            double diffPz=psi2S_vector.Pz()-lab4momentum_pz;
                            double diffE=psi2S_vector.E()-lab4momentum_e;
                            h_diffPsi_px->Fill(diffPx);
                            h_diffPsi_py->Fill(diffPy);
                            h_diffPsi_pz->Fill(diffPz);
                            h_diffPsi_e->Fill(diffE);
                            //if (isGoodPsi2S(diffPx, diffPy, diffPz, diffE)) 
                            {
                              efficientElectron=true;
                            }
                          }
                        }
                      } // Loop over last electron
 
                      if (efficientElectron)
                      {
                        h_e_eff_E_p->Fill(e2_vector_exp.P()*e2_vector_exp.P(), e2_vector_exp.E()*e2_vector_exp.E());
                        h_e_eff_Mass2->Fill(e2_vector_exp.M2());
                        /*
                        h_jPsi_Reco_Truth_p->Fill(jPsi_vector_MC.P(), jPsi_vector.P());
                        h_pi0_Reco_Truth_p->Fill(pi0_vector_MC.P(), pi0_vector.P());
                        h_e1_Reco_Truth_p->Fill(e1_vector_MC.P(), e1_vector.P());
                        h_ph_Reco_Truth_p->Fill(ph_vector_MC.P(), ph_vector.P());
                        h_e2_RecoMissing_MC_p->Fill(e2_vector_MC.P(), e2_vector_exp.P());
                        h_e2_Missing_MC_p->Fill(e2_vector_MC.P(), e2_vector_exp_MC.P());
                        */
                      }
                      else
                      {
                        h_e_ineff_E_p->Fill(e2_vector_exp.P()*e2_vector_exp.P(), e2_vector_exp.E()*e2_vector_exp.E());
                        h_e_ineff_Mass2->Fill(e2_vector_exp.M2());
                        h_pi0_2_missing_p_e2_missing_m->Fill(e2_vector_exp.M2(), pi0_2_missing.P());
                        if (e2_vector_exp.M2()>0.002 && e2_vector_exp.M2()<0.007)
                        {
                          std::cout<<runNumber<<" "<<eventNumber<<std::endl;
                          std::cout<<" J/Psi = "<<jPsi_vector.Px()<<", "<<jPsi_vector.Py()<<", "<<jPsi_vector.Pz()<<", "<<jPsi_vector.E()<<std::endl;
                          std::cout<<" 1st pi0 = "<<pi0_vector.Px()<<", "<<pi0_vector.Py()<<", "<<pi0_vector.Pz()<<", "<<pi0_vector.E()<<std::endl;
                          std::cout<<"  photon = "<<ph_vector.Px()<<", "<<ph_vector.Py()<<", "<<ph_vector.Pz()<<", "<<ph_vector.E()<<std::endl;
                          std::cout<<"  1st electron = "<<e1_vector.Px()<<", "<<e1_vector.Py()<<", "<<e1_vector.Pz()<<", "<<e1_vector.E()<<std::endl;
                          for (unsigned int i_lastElectron=0; i_lastElectron<e2_n; ++i_lastElectron)
                          {
                            TLorentzVector e2_vector(e2_px[i_lastElectron], e2_py[i_lastElectron], e2_pz[i_lastElectron], e2_e[i_lastElectron]);
                            TLorentzVector pi0_2_vector=e1_vector+e2_vector+ph_vector;
                            TLorentzVector psi2S_vector=jPsi_vector+pi0_vector+pi0_2_vector;
                            std::cout<<"  2nd electron = "<<e2_vector.Px()<<", "<<e2_vector.Py()<<", "<<e2_vector.Pz()<<", "<<e2_vector.E()<<std::endl;
                            std::cout<<" 2nd pi0 = "<<pi0_2_vector.Px()<<", "<<pi0_2_vector.Py()<<", "<<pi0_2_vector.Pz()<<", "<<pi0_2_vector.E()<<std::endl;
                            std::cout<<" psi(2S) = "<<psi2S_vector.Px()<<", "<<psi2S_vector.Py()<<", "<<psi2S_vector.Pz()<<", "<<psi2S_vector.E()<<std::endl;
                            std::cout<<" psi(2S) - J/psi mass = "<<psi2S_vector.M()-jPsi_vector.M()<<std::endl;
                          }
                        }
                        /*
                        h_jPsi_Reco_Truth_p->Fill(jPsi_vector_MC.P(), jPsi_vector.P());
                        h_pi0_Reco_Truth_p->Fill(pi0_vector_MC.P(), pi0_vector.P());
                        h_e1_Reco_Truth_p->Fill(e1_vector_MC.P(), e1_vector.P());
                        h_ph_Reco_Truth_p->Fill(ph_vector_MC.P(), ph_vector.P());
                        h_e2_RecoMissing_MC_p->Fill(e2_vector_MC.P(), e2_vector_exp.P());
                        h_e2_Missing_MC_p->Fill(e2_vector_MC.P(), e2_vector_exp_MC.P());
                        */
                      }
                    } // Missing Electron Energy and Angle Cut
                  } // Cleaning Last Photon
                } // Second pi0's Electron and Photon MC Matched
              }   // Missing mass of 2nd pi0 was acceptable
            } // pi0 is good
          } // First pi0 is MC Matched
        } // J/Psi is good
      } // J/Psi is MC Matched
    } // J/Psi is composed of electrons / muons
  } // Iterate over events
  
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
  pi01_mass->Divide(1,6);
  pi01_mass->cd(1);
  h_pi0Mass_RawPull->Draw("box");
  pi01_mass->cd(2);
  h_pi0Mass_Pull->Draw();
  line = new TLine(pi0PullMass_cut_, h_pi0Mass_Pull->GetMinimum(), pi0PullMass_cut_, h_pi0Mass_Pull->GetMaximum()*0.90); line->Draw();
  line = new TLine(-pi0PullMass_cut_, h_pi0Mass_Pull->GetMinimum(), -pi0PullMass_cut_, h_pi0Mass_Pull->GetMaximum()*0.90); line->Draw();
  pi01_mass->cd(3);
  h_pi0_hiShower_E->Draw();
  pi01_mass->cd(4);
  h_pi0_loShower_E->Draw();
  pi01_mass->cd(5);
  h_pi0_RecoMC_hiShower->Draw("box");
  pi01_mass->cd(6);
  h_pi0_RecoMC_loShower->Draw("box");
  
  TCanvas *c_pi0_Missingpi0 = new TCanvas("c_pi0_Missingpi0", "c_pi0_Missingpi0");
  h_Missingpi0_pi0->Draw("box");
  
  TCanvas *pi02_mass = new TCanvas("pi02_mass", "pi02_mass", 500, 1000);
  pi02_mass->Divide(1,5);
  pi02_mass->cd(1);
  h_pi0_2_missing->Draw();
  line = new TLine(0, 0, 0, h_pi0_2_missing->GetMaximum()*0.90); line->Draw();
  line = new TLine(pi0_2_missingMass_max, 0, pi0_2_missingMass_max, h_pi0_2_missing->GetMaximum()*0.90); line->Draw();
  pi02_mass->cd(2);
  h_pi0_2_Reco_m->Draw();
  line = new TLine(pi0Mass_pdg_-pi0Mass_reco_range, 0, pi0Mass_pdg_-pi0Mass_reco_range, h_pi0_2_Reco_m->GetMaximum()*0.90); line->Draw();
  line = new TLine(pi0Mass_pdg_+pi0Mass_reco_range, 0, pi0Mass_pdg_+pi0Mass_reco_range, h_pi0_2_Reco_m->GetMaximum()*0.90); line->Draw();
  pi02_mass->cd(3);
  h_pi0_2_Reco_m_e1diff->Draw("box");
  pi02_mass->cd(4);
  h_pi0_2_Reco_m_phdiff->Draw("box");
  pi02_mass->cd(5);
  h_pi0_2_missing_p->Draw();
  
  TCanvas *c_diffPsi = new TCanvas("c_diffPsi", "c_diffPsi", 500, 1000);
  c_diffPsi->Divide(1,5);
  c_diffPsi->cd(1);
  h_diffPsi->Draw();
  line = new TLine(diffPsi_center-diffPsi_range, 0, diffPsi_center-diffPsi_range, h_diffPsi->GetMaximum()*0.90); line->Draw();
  line = new TLine(diffPsi_center+diffPsi_range, 0, diffPsi_center+diffPsi_range, h_diffPsi->GetMaximum()*0.90); line->Draw();
  c_diffPsi->cd(2);
  h_diffPsi_diffpi0->Draw("box");
  c_diffPsi->cd(3);
  h_diffPsi_pi0_2_Reco_m->Draw("box");
  c_diffPsi->cd(4);
  h_diffPsi_hiShower->Draw("box");
  c_diffPsi->cd(5);
  h_diffPsi_loShower->Draw("box");
  
  TCanvas *Psi2S = new TCanvas("Psi2S", "Psi2S", 1000, 1000);
  Psi2S->Divide(2,2);
  Psi2S->cd(1);
  h_Psi2S_m->Draw();
  xmin=psi2SMass_center-psi2SMass_range;
  xmax=psi2SMass_center+psi2SMass_range;
  ymax=(h_Psi2S_m->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  Psi2S->cd(2);
  h_Psi2S_px->Draw();
  xmin=-psi2S_pxCut;
  xmax=psi2S_pxCut;
  ymax=(h_Psi2S_px->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  Psi2S->cd(3);
  h_Psi2S_py->Draw();
  xmin=-psi2S_pyCut;
  xmax=psi2S_pyCut;
  ymax=(h_Psi2S_py->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  Psi2S->cd(4);
  h_Psi2S_pz->Draw();
  xmin=-psi2S_pzCut;
  xmax=psi2S_pzCut;
  ymax=(h_Psi2S_pz->GetMaximum())*0.75;
  line = new TLine(xmin,ymax,xmin,0); line->Draw();
  line = new TLine(xmax,ymax,xmax,0); line->Draw();
  
  TCanvas *c_RelativePsi = new TCanvas("c_RelativePsi", "c_RelativePsi", 1000, 1000);
  c_RelativePsi->Divide(2,2);
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
  c_RelativePsi->cd(4);
  h_diffPsi_e->Draw();
  
  TCanvas *c_e2Missing_p=new TCanvas("c_e2Missing_p", "c_e2Missing_p");
  c_e2Missing_p->Divide(1,5);
  c_e2Missing_p->cd(1);
  h_e2Missing_e->Draw();
  ymax=(h_e2Missing_e->GetMaximum())*0.9;
  line = new TLine(missingElectronEnergy_min,ymax,missingElectronEnergy_min,0); line->Draw();
  line = new TLine(missingElectronEnergy_max,ymax,missingElectronEnergy_max,0); line->Draw();
  c_e2Missing_p->cd(2);
  h_e2Missing_p->Draw();
  ymax=(h_e2Missing_p->GetMaximum())*0.9;
  line = new TLine(missingElectronMomentum_min,ymax,missingElectronMomentum_min,0); line->Draw();
  line = new TLine(missingElectronMomentum_max,ymax,missingElectronMomentum_max,0); line->Draw();
  c_e2Missing_p->cd(3);
  h_e2Missing_cosTheta->Draw();
  ymax=(h_e2Missing_cosTheta->GetMaximum())*0.9;
  line = new TLine(-cosThetaMax_, ymax, -cosThetaMax_, 0); line->Draw();
  line = new TLine(cosThetaMax_, ymax, cosThetaMax_, 0); line->Draw();
  c_e2Missing_p->cd(4);
  h_e2_p_MC->Draw();
  c_e2Missing_p->cd(5);
  h_nMissingElectrons->Draw();
  
  TCanvas *c_Tags = new TCanvas("c_Tags", "c_Tags");
  c_Tags->Divide(1,6);
  c_Tags->cd(1);
  h_jPsi_Reco_Truth_p->Draw("box");
  c_Tags->cd(2);
  h_pi0_Reco_Truth_p->Draw("box");
  c_Tags->cd(3);
  h_e1_Reco_Truth_p->Draw("box");
  c_Tags->cd(4);
  h_ph_Reco_Truth_p->Draw("box");
  c_Tags->cd(5);
  h_e2_Missing_MC_p->Draw("box");
  c_Tags->cd(6);
  h_e2_RecoMissing_MC_p->Draw("box");
  
  TCanvas *c_x925 = new TCanvas("c_x925", "c_x925");
  c_x925->Divide(1,2);
  c_x925->cd(1);
  h_x925->Draw();
  ymax=h_x925->GetMaximum()*0.9;
  line = new TLine(x925Unf_min_,ymax,x925Unf_min_,0); line->Draw();
  line = new TLine(x925Unf_max_,ymax,x925Unf_max_,0); line->Draw();
  c_x925->cd(2);  
  h_ph_RecoTruth_x925->Draw("box");
  ymax=h_ph_RecoTruth_x925->GetMaximum();
  line = new TLine(x925Unf_min_,-0.5,x925Unf_min_,0.5); line->Draw();
  line = new TLine(x925Unf_max_,-0.5,x925Unf_max_,0.5); line->Draw();
  
  TCanvas *c_e2Candidate = new TCanvas("c_e2Candidate", "c_e2Candidate");
  c_e2Candidate->Divide(1,2);
  c_e2Candidate->cd(1);
  h_e2Candidate_p->Draw();
  c_e2Candidate->cd(2);
  h_e2Candidate_e->Draw();
  
  TCanvas *c_e2_MissingMass2_Resolution = new TCanvas("c_e2_MissingMass2_Resolution", "c_e2_MissingMass2_Resolution");
  c_e2_MissingMass2_Resolution->Divide(1,3);
  c_e2_MissingMass2_Resolution->cd(1);
  h_e2_MissingMass2_Resolution->Draw();
  c_e2_MissingMass2_Resolution->cd(2);
  h_e2_MissingMass2_MissingP->Draw("box");
  c_e2_MissingMass2_Resolution->cd(3);
  h_e2_MissingMass2_MissingE->Draw("box");
  
  TCanvas *c_e_E_p = new TCanvas("c_e_E_p", "c_e_E_p");
  c_e_E_p->Divide(1,2);
  c_e_E_p->cd(1);
  h_e_eff_E_p->Draw("box");
  c_e_E_p->cd(2);
  h_e_ineff_E_p->Draw("box");
  
  TCanvas *electronMass2 = new TCanvas("electronMass2", "electronMass2", 500, 1000);
  electronMass2->Divide(1,3);
  electronMass2->cd(1);
  h_e_eff_Mass2->Draw();
  /*
  TF1 *f_eff=new TF1("f_eff", e_eff_Fit, -0.2, 0.2, 9);
  f_eff->SetLineWidth(1);
  f_eff->SetLineColor(kRed);
  f_eff->SetParLimits(1, -0.004, 0.004);
  f_eff->SetParLimits(2, 0.001, 0.02);
  f_eff->SetParLimits(4, -0.004, 0.004);
  f_eff->SetParLimits(5, 0.001, 0.08);
  f_eff->SetParLimits(7, -0.004, 0.004);
  f_eff->SetParLimits(8, 0.0005, 0.01);
  h_e_eff_Mass2->Fit(f_eff, "RLEM");
  line = new TLine(mee_center-mee_range, 0, mee_center-mee_range, h_e_eff_Mass2->GetMaximum()*0.75); line->Draw();
  line = new TLine(mee_center+mee_range, 0, mee_center+mee_range, h_e_eff_Mass2->GetMaximum()*0.75); line->Draw();
  TF1 *f_e_eff_Peak=new TF1("f_e_eff_Peak", e_eff_Peak, -0.2, 0.2, 6);
  f_e_eff_Peak->SetParameter(0, f_eff->GetParameter(0));
  f_e_eff_Peak->SetParameter(1, f_eff->GetParameter(1));
  f_e_eff_Peak->SetParameter(2, f_eff->GetParameter(2));
  f_e_eff_Peak->SetParameter(3, f_eff->GetParameter(6));
  f_e_eff_Peak->SetParameter(4, f_eff->GetParameter(7));
  f_e_eff_Peak->SetParameter(5, f_eff->GetParameter(8));
  double effEvents=f_e_eff_Peak->Integral(mee_center-mee_range, mee_center+mee_range)*250;
  
  TF1 *f_e_eff_Back=new TF1("f_e_eff_Peak", e_eff_Back, -0.2, 0.2, 3);
  f_e_eff_Back->SetLineWidth(1);
  f_e_eff_Back->SetLineColor(kBlue);
  f_e_eff_Back->SetParameter(0, f_eff->GetParameter(3));
  f_e_eff_Back->SetParameter(1, f_eff->GetParameter(4));
  f_e_eff_Back->SetParameter(2, f_eff->GetParameter(5));
  f_e_eff_Back->Draw("SAME");
  */
  
  electronMass2->cd(2);
  h_e_ineff_Mass2->Draw();
  /*
  TF1 *f_e_ineff_Mass2=new TF1("f_e_ineff_Mass2", e_ineff_Fit, -0.2, 0.2, 12);
  f_e_ineff_Mass2->SetLineWidth(1);
  f_e_ineff_Mass2->SetLineColor(kRed);
  f_e_ineff_Mass2->SetParLimits(0, 0.0, 1e5);
  f_e_ineff_Mass2->SetParLimits(1, -0.004, 0.004);
  f_e_ineff_Mass2->SetParLimits(2, 0.001, 0.02);
  f_e_ineff_Mass2->SetParLimits(4, 0.01, 0.03);
  f_e_ineff_Mass2->SetParLimits(5, 0.01, 0.08);
  f_e_ineff_Mass2->SetParLimits(7, -0.01, 0.004);
  f_e_ineff_Mass2->SetParLimits(8, 0.01, 0.1);
  f_e_ineff_Mass2->SetParLimits(10, -0.004, 0.004);
  f_e_ineff_Mass2->SetParLimits(11, 0.0005, 0.01);
  h_e_ineff_Mass2->Fit(f_e_ineff_Mass2, "RLEM");
  line = new TLine(mee_center-mee_range, 0, mee_center-mee_range, h_e_ineff_Mass2->GetMaximum()*0.75); line->Draw();
  line = new TLine(mee_center+mee_range, 0, mee_center+mee_range, h_e_ineff_Mass2->GetMaximum()*0.75); line->Draw();
  TF1 *f_e_ineff_Peak=new TF1("f_e_ineff_Peak", e_ineff_Peak, -0.2, 0.2, 6);
  f_e_ineff_Peak->SetParameter(0, f_e_ineff_Mass2->GetParameter(0));
  f_e_ineff_Peak->SetParameter(1, f_e_ineff_Mass2->GetParameter(1));
  f_e_ineff_Peak->SetParameter(2, f_e_ineff_Mass2->GetParameter(2));
  f_e_ineff_Peak->SetParameter(3, f_e_ineff_Mass2->GetParameter(9));
  f_e_ineff_Peak->SetParameter(4, f_e_ineff_Mass2->GetParameter(10));
  f_e_ineff_Peak->SetParameter(5, f_e_ineff_Mass2->GetParameter(11));
  double ineffEvents=f_e_ineff_Peak->Integral(mee_center-mee_range, mee_center+mee_range)*250;
  TF1 *f_e_ineff_Back=new TF1("f_e_ineff_Back", e_ineff_Back, -0.2, 0.2, 6);
  f_e_ineff_Back->SetLineWidth(1);
  f_e_ineff_Back->SetLineColor(kBlue);
  f_e_ineff_Back->SetParameter(0, f_e_ineff_Mass2->GetParameter(3));
  f_e_ineff_Back->SetParameter(1, f_e_ineff_Mass2->GetParameter(4));
  f_e_ineff_Back->SetParameter(2, f_e_ineff_Mass2->GetParameter(5));
  f_e_ineff_Back->SetParameter(3, f_e_ineff_Mass2->GetParameter(6));
  f_e_ineff_Back->SetParameter(4, f_e_ineff_Mass2->GetParameter(7));
  f_e_ineff_Back->SetParameter(5, f_e_ineff_Mass2->GetParameter(8));
  f_e_ineff_Back->Draw("SAME");
  */
  
  electronMass2->cd(3);
  h_pi0_2_missing_p_e2_missing_m->Draw("box");
  
  /*
  std::cout<<"Number of events under efficient peak = "<<effEvents<<std::endl;
  std::cout<<"Number of events under inefficient peak = "<<ineffEvents<<std::endl;
  double efficiency=effEvents/(effEvents+ineffEvents);
  std::cout<<"Electron reconstruction efficiency = "<<efficiency<<std::endl;
  */
  
  /*
  TCanvas *angle = new TCanvas("angle");
  angle->Divide(1,2);
  angle->cd(1);
  h_eff_Angle->Draw();
  angle->cd(2);
  h_ineff_Angle->Draw();
  */
  return 0;
}
