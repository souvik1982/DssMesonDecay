#include <iostream>
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TFile.h"
#include "TChain.h"
#include "TPaveStats.h"
#include <set>
#include <map>

int DsTagCompare()
{
  //TFile *originFile = new TFile("/nfs/cor/an3/souvik/Dataset39/DsTaggedProc_Data_213586_214863.root");
  //TFile *originFile = new TFile("/nfs/cor/an3/souvik/Dataset40/DsTaggedProc_Data_215307_217385.root");
  //TFile *originFile = new TFile("/nfs/cor/an3/souvik/Dataset41/DsTaggedProc_Data_217687_219721.root");
  //TFile *originFile = new TFile("/nfs/cor/an3/souvik/Dataset47/DsTaggedProc_Data_230474_232255.root");
  //TFile *originFile = new TFile("/nfs/cor/an3/souvik/Dataset48/DsTaggedDecaysProc_OriginalData_data48.root");
  //originFile->cd("DsTaggedDecaysProc");
  //TTree *originTree = (TTree*)gDirectory->Get("nt4");
  
  TChain *originTree=new TChain("DsTaggedDecaysProc/nt4");
  originTree->Add("/nfs/cor/an3/souvik/Dataset39/DsTaggedProc_Data_213586_214863.root");
  originTree->Add("/nfs/cor/an3/souvik/Dataset40/DsTaggedProc_Data_215307_217385.root");
  originTree->Add("/nfs/cor/an3/souvik/Dataset41/DsTaggedProc_Data_217687_219721.root");
  originTree->Add("/nfs/cor/an3/souvik/Dataset47/DsTaggedProc_Data_230474_232255.root");
  originTree->Add("/nfs/cor/an3/souvik/Dataset48/DsTaggedDecaysProc_OriginalData_data48.root");
  
  //TFile *reProcFile = new TFile("/nfs/cor/an3/souvik/Dataset39/DsTaggedDecaysProc_ReTaggedData_213586_214863.root");
  //TFile *reProcFile = new TFile("/nfs/cor/an3/souvik/Dataset40/DsTaggedDecaysProc_ReTaggedData_215307_217385.root");
  //TFile *reProcFile = new TFile("/nfs/cor/an3/souvik/Dataset41/DsTaggedDecaysProc_ReTaggedData_217687_219721.root");
  //TFile *reProcFile = new TFile("/nfs/cor/an3/souvik/Dataset47/DsTaggedDecaysProc_ReTaggedData_230474_232255.root");
  //TFile *reProcFile = new TFile("/nfs/cor/an3/souvik/Dataset48/DsTaggedDecaysProc_ReTaggedData_232264_234607.root");
  //reProcFile->cd("DsTaggedDecaysProc");
  //TTree *reProcTree = (TTree*)gDirectory->Get("nt4");
  
  TChain *reProcTree=new TChain("DsTaggedDecaysProc/nt4");
  reProcTree->Add("/nfs/cor/an3/souvik/Dataset39/DsTaggedDecaysProc_ReTaggedData_213586_214863.root");
  reProcTree->Add("/nfs/cor/an3/souvik/Dataset40/DsTaggedDecaysProc_ReTaggedData_215307_217385.root");
  reProcTree->Add("/nfs/cor/an3/souvik/Dataset41/DsTaggedDecaysProc_ReTaggedData_217687_219721.root");
  reProcTree->Add("/nfs/cor/an3/souvik/Dataset47/DsTaggedDecaysProc_ReTaggedData_230474_232255.root");
  reProcTree->Add("/nfs/cor/an3/souvik/Dataset48/DsTaggedDecaysProc_ReTaggedData_232264_234607.root");
  
  float run_origin, event_origin, decayMode_origin, dsPlusM_origin, dsPlusE_origin, dsPlusPx_origin, dsPlusPy_origin, dsPlusPz_origin;
  float run_reProc, event_reProc, decayMode_reProc, dsPlusM_reProc, dsPlusE_reProc, dsPlusPx_reProc, dsPlusPy_reProc, dsPlusPz_reProc;
  float beamEnergy_origin;
  float beamEnergy_reProc;
  float pi0RawMass_origin, pi0PullMass_origin, etaRawMass_origin, etaPullMass_origin;
  float pi0RawMass_reProc, pi0PullMass_reProc, etaRawMass_reProc, etaPullMass_reProc;
  
  TH1F *h_dsPlusM = new TH1F("h_dsPlusM", "h_dsPlusM", 100, -0.000001, 0.000001);
  TH2F *h_dsPlusM_event = new TH2F("h_dsPlusM_event", "h_dsPlusM_event", 243377, 0, 243377, 100, -0.00001, 0.00001);
  TH1F *h_beamEnergy = new TH1F("h_beamEnergy", "h_beamEnergy", 100, -0.0001, 0.0001);
  
  // Original File
  originTree->SetBranchAddress("Run", &(run_origin));
  originTree->SetBranchAddress("Event", &(event_origin));
  originTree->SetBranchAddress("DecayMode", &(decayMode_origin));
  originTree->SetBranchAddress("dsPlusM", &(dsPlusM_origin));
  originTree->SetBranchAddress("BeamEnergy", &(beamEnergy_origin));
  originTree->SetBranchAddress("kPi0RawMass", &(pi0RawMass_origin));
  originTree->SetBranchAddress("kPi0PullMass", &(pi0PullMass_origin));
  originTree->SetBranchAddress("kEtaRawMass", &(etaRawMass_origin));
  originTree->SetBranchAddress("kEtaPullMass", &(etaPullMass_origin));
  
  // ReProcessed File
  reProcTree->SetBranchAddress("Run", &(run_reProc));
  reProcTree->SetBranchAddress("Event", &(event_reProc));
  reProcTree->SetBranchAddress("DecayMode", &(decayMode_reProc));
  reProcTree->SetBranchAddress("dsPlusM", &(dsPlusM_reProc));
  reProcTree->SetBranchAddress("BeamEnergy", &(beamEnergy_reProc));
  reProcTree->SetBranchAddress("kPi0RawMass", &(pi0RawMass_reProc));
  reProcTree->SetBranchAddress("kPi0PullMass", &(pi0PullMass_reProc));
  reProcTree->SetBranchAddress("kEtaRawMass", &(etaRawMass_reProc));
  reProcTree->SetBranchAddress("kEtaPullMass", &(etaPullMass_reProc));
  
  bool details=true;
  int i=0, k=0;
  int nMissedTags=0, nNewTags=0;
  int totMissedTags=0, totNewTags=0;
  int oldRun=0;
  int nEvents_origin=0, nEvents_reProc=0;
  int old_event_origin=-1, old_event_reProc=-1;
  int nDTags_origin=0, nDTags_reProc=0;
  int nEntries_origin=originTree->GetEntries();
  
  while (run_origin<=234607 && i<nEntries_origin)
  {
    originTree->GetEvent(i);
    reProcTree->GetEvent(i+k);
    // if (run_origin>234607) break;
    
    if (run_origin!=oldRun && run_reProc!=oldRun)
    {
      std::cout<<"# of events in original = "<<nEvents_origin<<" || # of events in reprocessed = "<<nEvents_reProc<<std::endl;
      std::cout<<"# of Ds tags in original = "<<i-nDTags_origin<<" || # of Ds tags in reprocessed = "<<i+k-nDTags_reProc<<std::endl;
      std::cout<<"Missed tags = "<<nMissedTags<<", New tags = "<<nNewTags<<std::endl<<std::endl;
      std::cout<<"New Run of original file "<<run_origin<<std::endl;
      std::cout<<"========================"<<std::endl;
      oldRun=run_origin;
      nEvents_origin=0; nEvents_reProc=0;
      old_event_origin=-1; old_event_reProc=-1;
      nDTags_origin=i; nDTags_reProc=i+k;
      nMissedTags=0; nNewTags=0;
      --i;
    }
    else
    {
      /*
      if (run_origin>214547)
      {
        std::cout<<"(run_origin = "<<run_origin<<", event_origin = "<<event_origin<<", i = "<<i<<") ||"
                 <<" (run_reProc = "<<run_reProc<<", event_reProc = "<<event_reProc<<", i+k = "<<i+k<<")"<<std::endl;
        std::cout<<" Original Ds mass = "<<dsPlusM_origin<<"Decay mode = "<<decayMode_origin<<" || Reprocessed Ds mass = "<<dsPlusM_reProc<<", Decay mode = "<<decayMode_reProc<<std::endl;
        std::cout<<" beamEnergy_origin = "<<beamEnergy_origin<<"|| beamEnergy_reProc = "<<beamEnergy_reProc<<std::endl;
      }
      */
      
      if (run_reProc>run_origin)
      {
        --k;
        std::cout<<"event_origin = "<<event_origin<<" Missed a DTag!"<<std::endl; ++nMissedTags; ++totMissedTags;
      }
      else if (run_origin>run_reProc)
      {
        ++k; --i;
        std::cout<<"event_reProc = "<<event_reProc<<" New DTag!"<<std::endl; ++nNewTags; ++totNewTags;
      }
      else if (event_reProc>event_origin) {
        --k;
        std::cout<<"event_origin = "<<event_origin<<" Missed a DTag!"<<std::endl; ++nMissedTags; ++totMissedTags;
      }
      else if (event_origin>event_reProc) {
        ++k; --i;
        std::cout<<"event_reProc = "<<event_reProc<<" New DTag!"<<std::endl; ++nNewTags; ++totNewTags;
      }
      
      /*
      else {
        h_dsPlusM->Fill(dsPlusM_origin-dsPlusM_reProc);
        h_dsPlusM_event->Fill(i, dsPlusM_origin-dsPlusM_reProc);
        h_beamEnergy->Fill(beamEnergy_origin-beamEnergy_reProc);
      }
      */
      
      if (event_origin>old_event_origin)
      {
        ++nEvents_origin;
        old_event_origin=event_origin;
      }
      if (event_reProc>old_event_reProc)
      {
        ++nEvents_reProc;
        old_event_reProc=event_reProc;
      }
      
    }
    
    ++i;
    
  }
  
  std::cout<<"# of events in original = "<<nEvents_origin<<" || # of events in reprocessed = "<<nEvents_reProc<<std::endl;
  std::cout<<"# of Ds tags in original = "<<i-nDTags_origin<<" || # of Ds tags in reprocessed = "<<i+k-nDTags_reProc<<std::endl;
  std::cout<<"Missed tags = "<<nMissedTags<<", New tags = "<<nNewTags<<std::endl<<std::endl;
  std::cout<<"========================"<<std::endl;
  std::cout<<"========================"<<std::endl;
  
  std::cout<<"Total DTags in original = "<<i<<", Total DTags in reprocessed = "<<i+k<<std::endl;
  std::cout<<"Total Missed DTags = "<<totMissedTags<<", Total New DTags = "<<totNewTags<<std::endl;
  
  /*
  TCanvas *c_dsPlusM = new TCanvas("c_dsPlusM");
  c_dsPlusM->Divide(1,2);
  c_dsPlusM->cd(1);
  h_dsPlusM->Draw();
  c_dsPlusM->cd(2);
  h_dsPlusM_event->Draw("box");
  
  TCanvas *c_beamEnergy = new TCanvas("c_beamEnergy");
  h_beamEnergy->Draw();
  */
  
  return 0;
}

/* if (decayMode_origin==422 || decayMode_origin==501)
    {
      std::cout<<"Origin pi0RawMass = "<<pi0RawMass_origin<<", pi0Pull Mass = "<<pi0PullMass_origin<<std::endl;
    }
    if (decayMode_reProc==422 || decayMode_reProc==501)
    {
      std::cout<<"reProc pi0RawMass = "<<pi0RawMass_reProc<<", pi0Pull Mass = "<<pi0PullMass_reProc<<std::endl;
    }
    if (decayMode_origin==442)
    {
      std::cout<<"Origin etaRawMass = "<<etaRawMass_origin<<", etaPull Mass = "<<etaPullMass_origin<<std::endl;
    }
    if (decayMode_reProc==442)
    {
      std::cout<<"reProc etaRawMass = "<<etaRawMass_reProc<<", etaPull Mass = "<<etaPullMass_reProc<<std::endl;
    } */
    /* if (beamEnergy_origin!=beamEnergy_reProc)
    {
      std::cout<<"run_origin = "<<run_origin<<", run_reProc = "<<run_reProc<<std::endl;
      std::cout<<"i = "<<i<<" event_origin = "<<event_origin<<", i+k = "<<i+k<<" event_reProc = "<<event_reProc<<std::endl;
      std::cout<<"beamEnergy_origin = "<<beamEnergy_origin<<", beamEnergy_reProc = "<<beamEnergy_reProc<<std::endl;
    } */
