int ComparePass2()
{
  TFile *orFile = new TFile("DssteeBIGFILES/DsTaggedProc_Data_230474_232255.root");
  orFile->cd("DsTaggedDecaysProc");
  TTree *orTree = (TTree*)gDirectory->Get("nt5");

  TFile *reFile = new TFile("DssteeBIGFILES/DsTaggedDecaysProc_ReTaggedData_230474_230617.root");
  reFile->cd("DsTaggedDecaysProc");
  TTree *reTree = (TTree*)gDirectory->Get("nt5");
  
  orTree->Draw("showerNumber", "Run==230474 && Event==1700");
  reTree->Draw("showerNumber", "Run==230474 && Event==1700");

  
}


