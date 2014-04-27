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
#include "TPaveStats.h"
#include <set>
#include <map>

int CompareTracksShowers()
{

  TFile *originFile = new TFile("/nfs/cor/an2/souvik/Dataset48/DsTaggedDecaysProc_OriginalData_data48.root");
  originFile->cd("DsTaggedDecaysProc");
  TTree *originTree = (TTree*)gDirectory->Get("nt5");
  
  TFile *reProcFile = new TFile("/nfs/cor/an2/souvik/Dataset48/DsTaggedDecaysProc_ReTaggedData_232264_234607.root");
  reProcFile->cd("DsTaggedDecaysProc");
  TTree *reProcTree = (TTree*)gDirectory->Get("nt5");
  
  float run_origin, event_origin, trackNumber_origin, trackPx_origin, trackPy_origin, trackPz_origin, showerNumber_origin, showerE_origin;
  float run_reProc, event_reProc, trackNumber_reProc, trackPx_reProc, trackPy_reProc, trackPz_reProc, showerNumber_reProc, showerE_reProc;
  
  TH1F *h_dTrackPx = new TH1F("h_dTrackPx", "h_dTrackPx", 100, -0.0001, 0.0001);
  TH1F *h_dTrackPy = new TH1F("h_dTrackPy", "h_dTrackPy", 100, -0.0001, 0.0001);
  TH1F *h_dTrackPz = new TH1F("h_dTrackPz", "h_dTrackPz", 100, -0.0001, 0.0001);
  TH1F *h_dShowerE = new TH1F("h_dShowerE", "h_dShowerE", 100, -0.000001, 0.000001);
  
  // Original File
  originTree->SetBranchAddress("Run", &(run_origin));
  originTree->SetBranchAddress("Event", &(event_origin));
  originTree->SetBranchAddress("trackNumber", &(trackNumber_origin));
  originTree->SetBranchAddress("trackPx", &(trackPx_origin));
  originTree->SetBranchAddress("trackPy", &(trackPy_origin));
  originTree->SetBranchAddress("trackPz", &(trackPz_origin));
  originTree->SetBranchAddress("showerNumber", &(showerNumber_origin));
  originTree->SetBranchAddress("showerE", &(showerE_origin));
  
  // ReProcessed File
  reProcTree->SetBranchAddress("Run", &(run_reProc));
  reProcTree->SetBranchAddress("Event", &(event_reProc));
  reProcTree->SetBranchAddress("trackNumber", &(trackNumber_reProc));
  reProcTree->SetBranchAddress("trackPx", &(trackPx_reProc));
  reProcTree->SetBranchAddress("trackPy", &(trackPy_reProc));
  reProcTree->SetBranchAddress("trackPz", &(trackPz_reProc));
  reProcTree->SetBranchAddress("showerNumber", &(showerNumber_reProc));
  reProcTree->SetBranchAddress("showerE", &(showerE_reProc));
  
  cout.precision(6);
  int k=0;
  int nMissedTags=0, nNewTags=0;
  for (int i=0; i<originTree->GetEntries(); ++i)
  {
    originTree->GetEvent(i);
    reProcTree->GetEvent(i+k);
    if (run_origin>232269) break;
    //std::cout<<"run_origin = "<<run_origin<<", run_reProc = "<<run_reProc<<std::endl;
    //std::cout<<"i = "<<i<<" event_origin = "<<event_origin<<", i+k = "<<i+k<<" event_reProc = "<<event_reProc<<std::endl;
    
    if (event_reProc>event_origin) 
    {
      --k;
      //std::cout<<"Missed tag!"<<std::endl;
      ++nMissedTags;
      
    }
    else if (event_reProc<event_origin) 
    {
      ++k;
      //std::cout<<"New tag!"<<std::endl;
      ++nNewTags;
    }
    else if (run_origin==232269 && event_origin==88989 && run_reProc==232269 && event_reProc==88989)
    {
      if (trackNumber_origin!=-1 && trackNumber_reProc!=-1)
      {
        if (trackPx_reProc!=trackPx_origin) h_dTrackPx->Fill(trackPx_reProc-trackPx_origin);
        if (trackPy_reProc!=trackPy_origin) h_dTrackPy->Fill(trackPy_reProc-trackPy_origin);
        if (trackPz_reProc!=trackPz_origin) h_dTrackPz->Fill(trackPz_reProc-trackPz_origin);
        //std::cout<<(trackPz_origin-trackPz_reProc)<<std::endl;
      }
      if (showerNumber_origin!=-1 && showerNumber_reProc!=-1)
      {
         if (showerE_origin!=showerE_reProc) h_dShowerE->Fill(showerE_reProc-showerE_origin);
         //std::cout<<(showerE_reProc-showerE_origin)<<std::endl;
      }
    }
    
  }
  
  std::cout<<"Original # entries = "<<originTree->GetEntries()<<std::endl;
  std::cout<<"Reprocessed # entries = "<<reProcTree->GetEntries()<<std::endl;
  std::cout<<"missed tags = "<<nMissedTags<<", new tags = "<<nNewTags<<std::endl;
  
  TCanvas *c_dTrack = new TCanvas("c_dTrack");
  c_dTrack->Divide(1,3);
  c_dTrack->cd(1);
  h_dTrackPx->Draw();
  c_dTrack->cd(2);
  h_dTrackPy->Draw();
  c_dTrack->cd(3);
  h_dTrackPz->Draw();
  
  TCanvas *c_dShower = new TCanvas("c_dShower");
  h_dShowerE->Draw();
  
  return 0;
}
