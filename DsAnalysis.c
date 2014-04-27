int DsAnalysis()
{
  
  bool vertexFitted=true;
  
  double pi=3.14159265358979;
  
  double signalScale=0.00936536;
  double converScale=0.0556227;
  
  TChain *signalTree;
  if (vertexFitted) signalTree=new TChain("DsTaggedDecaysProc/nt7");
  else signalTree=new TChain("DsTaggedDecaysProc/nt3");
  signalTree->Add("/nfs/cor/an2/souvik/MC_vtosll_Dsp_KKpi/DsTaggedDecaysProc_MC_vtosll_Dsp_KKpi.root");
  signalTree->Add("/nfs/cor/an2/souvik/MC_vtosll_Dsm_KKpi/DsTaggedDecaysProc_MC_vtosll_Dsm_KKpi.root");
  
  TChain *converTree;
  if (vertexFitted) converTree=new TChain("DsTaggedDecaysProc/nt7");
  else converTree=new TChain("DsTaggedDecaysProc/nt3");
  //converTree->Add("/nfs/cor/an2/souvik/MC_gamma_Dsp_KKpi/DsTaggedDecaysProc_MC_gamma_Dsp_KKpi.root");
  //converTree->Add("/nfs/cor/an2/souvik/MC_gamma_Dsm_KKpi/DsTaggedDecaysProc_MC_gamma_Dsm_KKpi.root");
  converTree->Add("/nfs/cor/an3/souvik/MC_gamma_Dsp_generic/DsTaggedDecaysProc_MC_gamma_Dsp_generic.root");
  converTree->Add("/nfs/cor/an3/souvik/MC_gamma_Dsm_generic/DsTaggedDecaysProc_MC_gamma_Dsm_generic.root");
  
  TCut decayMode       = "DecayMode==401";
  TCut electronQuality = "kElectron1E_reco<0.15 && kElectron2E_reco<0.15";
  if (vertexFitted)
  {
  TCut dsPlusMCut      = "abs(dsPlusM-1.96849)<0.011";
  TCut MBCCut          = "abs(MBC-2.112)<0.004";
  TCut deltaMCut       = "abs(DeltaM-0.1438)<0.004";
  TCut diffD0Cut       = "(kElectron1D0_reco-kElectron2D0_reco)>-0.005";
  TCut dPhiCut         = "(atan2(kElectron1Py_reco, kElectron1Px_reco)-atan2(kElectron2Py_reco, kElectron2Px_reco))<0.12";
  }
  else
  {
  TCut dsPlusMCut      = "abs(dsPlusM-1.96849)<0.011";
  TCut MBCCut          = "abs(MBC-2.112)<0.004";
  TCut deltaMCut       = "abs(DeltaM-0.1438)<0.006";
  TCut diffD0Cut       = "(kElectron1D0_reco-kElectron2D0_reco)>-0.006"; // -0.005
  TCut dPhiCut         = "(atan2(kElectron1Py_reco, kElectron1Px_reco)-atan2(kElectron2Py_reco, kElectron2Px_reco))<0.1"; // 0.12
  }
  TCut vtx             = "chisqVtx>0";
  TCut vtxCut          = "abs(0.013-pow(pow(kElectron1E_reco+kElectron2E_reco, 2)-(pow(kElectron1Px_reco+kElectron2Px_reco, 2)+pow(kElectron1Py_reco+kElectron2Py_reco, 2)+pow(kElectron1Pz_reco+kElectron2Pz_reco, 2)), 0.5))>0.00389";
  
  TH1D *h_mee_signal=new TH1D("h_mee_signal", "m_{ee} Signal MC; m_{ee} (GeV); Number of Events", 50, 0., 0.05); h_mee_signal->SetLineColor(kRed);
  TH1D *h_mee_conver=new TH1D("h_mee_conver", "m_{ee} Conversion MC; m_{ee} (GeV); Number of Events", 50, 0., 0.05); h_mee_conver->SetLineColor(kGreen);
  
  TH1D *h_diffD0_signal = new TH1D("h_diffD0_signal", "#Deltad_{0} Signal Sample; m", 50, -0.01, 0.01); h_diffD0_signal->SetLineColor(kRed);
  TH1D *h_diffD0_conver = new TH1D("h_diffD0_conver", "#Deltad_{0} Conversion Sample; m", 50, -0.01, 0.01); h_diffD0_conver->SetLineColor(kGreen);
  
  TH1D *h_dPhi_signal = new TH1D("h_dPhi_signal", "#Delta#Phi Signal Sample", 50, -pi/2, pi/2); h_dPhi_signal->SetLineColor(kRed);
  TH1D *h_dPhi_conver = new TH1D("h_dPhi_conver", "#Delta#Phi conver Sample", 50, -pi/2, pi/2); h_dPhi_conver->SetLineColor(kGreen);
  
  gROOT->SetStyle("Plain");
  
  if (vertexFitted)
  {
    signalTree->Draw("pow(pow(kElectron1E_reco+kElectron2E_reco, 2)-(pow(kElectron1Px_reco+kElectron2Px_reco, 2)+pow(kElectron1Py_reco+kElectron2Py_reco, 2)+pow(kElectron1Pz_reco+kElectron2Pz_reco, 2)), 0.5)>>h_mee_signal",
                      //"(kElectron1D0_reco-kElectron2D0_reco)>>h_diffD0_signal",
                      //"(atan2(kElectron1Py_reco, kElectron1Px_reco)-atan2(kElectron2Py_reco, kElectron2Px_reco))>>h_dPhi_signal",
                                          decayMode && electronQuality && dsPlusMCut && MBCCut && deltaMCut && vtxCut);
    converTree->Draw("pow(pow(kElectron1E_reco+kElectron2E_reco, 2)-(pow(kElectron1Px_reco+kElectron2Px_reco, 2)+pow(kElectron1Py_reco+kElectron2Py_reco, 2)+pow(kElectron1Pz_reco+kElectron2Pz_reco, 2)), 0.5)>>h_mee_conver",
                      //"(kElectron1D0_reco-kElectron2D0_reco)>>h_diffD0_conver",
                      //"(atan2(kElectron1Py_reco, kElectron1Px_reco)-atan2(kElectron2Py_reco, kElectron2Px_reco))>>h_dPhi_conver",
                                          decayMode && electronQuality && dsPlusMCut && MBCCut && deltaMCut && vtxCut);
    //                                    decayMode && electronQuality && dsPlusMCut && MBCCut && deltaMCut && diffD0Cut && dPhiCut && vtxCut);
    std::cout<<"Vertex Fitted"<<std::endl;
  }
  else
  {
    signalTree->Draw("pow(pow(kElectron1E_reco+kElectron2E_reco, 2)-(pow(kElectron1Px_reco+kElectron2Px_reco, 2)+pow(kElectron1Py_reco+kElectron2Py_reco, 2)+pow(kElectron1Pz_reco+kElectron2Pz_reco, 2)), 0.5)>>h_mee_signal",
                      //"(kElectron1D0_reco-kElectron2D0_reco)>>h_diffD0_signal",
                      //"(atan2(kElectron1Py_reco, kElectron1Px_reco)-atan2(kElectron2Py_reco, kElectron2Px_reco))>>h_dPhi_signal",
                                          decayMode && electronQuality && dsPlusMCut && MBCCut && deltaMCut && diffD0Cut && dPhiCut);
    converTree->Draw("pow(pow(kElectron1E_reco+kElectron2E_reco, 2)-(pow(kElectron1Px_reco+kElectron2Px_reco, 2)+pow(kElectron1Py_reco+kElectron2Py_reco, 2)+pow(kElectron1Pz_reco+kElectron2Pz_reco, 2)), 0.5)>>h_mee_conver",
                      //"(kElectron1D0_reco-kElectron2D0_reco)>>h_diffD0_conver",
                      //"(atan2(kElectron1Py_reco, kElectron1Px_reco)-atan2(kElectron2Py_reco, kElectron2Px_reco))>>h_dPhi_conver",
                                          decayMode && electronQuality && dsPlusMCut && MBCCut && deltaMCut && diffD0Cut && dPhiCut);
    //                                    decayMode && electronQuality && dsPlusMCut && MBCCut && deltaMCut && diffD0Cut && dPhiCut);
    std::cout<<"Without Vertex Fitting"<<std::endl;
  }
  
  h_mee_signal->Scale(signalScale);
  h_mee_conver->Scale(converScale);
  
  double nSignal=h_mee_signal->Integral(0, 51);
  double nConver=h_mee_conver->Integral(0, 51);
  std::cout<<"Signal within = "<<nSignal<<std::endl;
  std::cout<<"Conver within = "<<nConver<<std::endl;
  std::cout<<"Significance = "<<nSignal/pow(nConver, 0.5)<<std::endl;
  
  
  
  TCanvas *c_mee=new TCanvas("c_mee", "", 800, 350);
  c_mee->Divide(2,1);
  c_mee->cd(1);
  //c_mee_1->SetLogx();
  h_mee_signal->Draw();
  TLine *line;
  line=new TLine(0.013-0.00389, 0, 0.013-0.00389, 1.6); line->Draw();
  line=new TLine(0.013+0.00389, 0, 0.013+0.00389, 1.6); line->Draw();
  c_mee->cd(2);
  //c_mee_2->SetLogx();
  h_mee_conver->Draw();
  line=new TLine(0.013-0.00389, 0, 0.013-0.00389, 1.4); line->Draw();
  line=new TLine(0.013+0.00389, 0, 0.013+0.00389, 1.4); line->Draw();
  
  
  /*
  TCanvas *c_diffD0=new TCanvas("c_diffD0");
  c_diffD0->Divide(1,2);
  c_diffD0->cd(1);
  h_diffD0_signal->Draw();
  c_diffD0->cd(2);
  h_diffD0_conver->Draw();
  */
  
  /*
  TCanvas *c_dPhi=new TCanvas("c_dPhi");
  c_dPhi->Divide(1,2);
  c_dPhi->cd(1);
  h_dPhi_signal->Draw();
  c_dPhi->cd(2);
  h_dPhi_conver->Draw();
  */
  
/*
  
  TCut dsPlusMCut   = "abs(dsPlusM-1.96849)<0.02";
  TCut DeltaECut    = "abs(DeltaE)<0.05";
  TCut MBCCut       = "abs(MBC-2.112)<0.005";  
  
  TH1D *dsPlusM_signal = new TH1D("dsPlusM_signal", "dsPlusM_signal", 100, 1.9, 2.1); dsPlusM_signal->SetLineColor(kRed);
  TH1D *dsPlusM_bg = new TH1D("dsPlusM_bg", "dsPlusM_bg", 100, 1.9, 2.1); dsPlusM_bg->SetLineColor(kBlue);
  signalTree->Draw("dsPlusM>>dsPlusM_signal");
  bgTree->Draw("dsPlusM>>dsPlusM_bg");
  
  TH1D *DeltaE_signal = new TH1D("DeltaE_signal", "DeltaE", 100, -0.1, 0.2); DeltaE_signal->SetLineColor(kRed);
  TH1D *DeltaE_bg = new TH1D("DeltaE_bg", "DeltaE", 100, -0.1, 0.2); DeltaE_bg->SetLineColor(kBlue);
  signalTree->Draw("DeltaE>>DeltaE_signal", dsPlusMCut);
  bgTree->Draw("DeltaE>>DeltaE_bg", dsPlusMCut);
  
  TH1D *MBC_signal = new TH1D("MBC_signal", "MBC", 100, 2., 2.2); MBC_signal->SetLineColor(kRed);
  TH1D *MBC_bg = new TH1D("MBC_bg","MBC", 100, 2., 2.2); MBC_bg->SetLineColor(kBlue);
  signalTree->Draw("MBC>>MBC_signal", dsPlusMCut && DeltaECut);
  bgTree->Draw("MBC>>MBC_bg", dsPlusMCut && DeltaECut);
  
  TH1D *DeltaM_signal = new TH1D("DeltaM_signal", "DeltaM", 100, 0.0, 0.2); DeltaM_signal->SetLineColor(kRed);
  TH1D *DeltaM_bg = new TH1D("DeltaM_bg", "DeltaM", 100, 0.0, 0.2); DeltaM_bg->SetLineColor(kBlue);
  signalTree->Draw("DeltaM>>DeltaM_signal", dsPlusMCut && DeltaECut && MBCCut);
  bgTree->Draw("DeltaM>>DeltaM_bg", dsPlusMCut && DeltaECut && MBCCut);
  
  TH1D *d0_e_signal = new TH1D("d0_e_signal", "d0_e", 50, -0.01, 0.01); d0_e_signal->SetLineColor(kRed);
  TH1D *d0_e_bg = new TH1D("d0_e_bg", "d0_e", 50, -0.01, 0.01); d0_e_bg->SetLineColor(kBlue);
  signalTree->Draw("kElectron1D0_reco>>d0_e_signal", phiCut && dsPlusMCut && DeltaECut && MBCCut);
  bgTree->Draw("kElectron1D0_reco>>d0_e_bg", phiCut && dsPlusMCut && DeltaECut && MBCCut);
  
  TH1D *d0_p_signal = new TH1D("d0_p_signal", "d0_p", 50, -0.01, 0.01); d0_p_signal->SetLineColor(kRed);
  TH1D *d0_p_bg = new TH1D("d0_p_bg", "d0_p", 50, -0.01, 0.01); d0_p_bg->SetLineColor(kBlue);
  signalTree->Draw("kElectron2D0_reco>>d0_p_signal", phiCut && dsPlusMCut && DeltaECut && MBCCut);
  bgTree->Draw("kElectron2D0_reco>>d0_p_bg", phiCut && dsPlusMCut && DeltaECut && MBCCut);
  
  // d0 of e + d0 of p
  
  
  
  TH1D *RatioHitsToExpected_electron_signal = new TH1D("RatioHitsToExpected_electron_signal", "RatioHitsToExpected_electron_signal", 10, 0., 1.); RatioHitsToExpected_electron_signal->SetLineColor(kRed);
  TH1D *RatioHitsToExpected_electron_bg = new TH1D("RatioHitsToExpected_electron_bg", "RatioHitsToExpected_electron_bg", 10, 0., 1.); RatioHitsToExpected_electron_bg->SetLineColor(kRed);
  signalTree->Draw("kElectron1RatioHitsToExpected>>RatioHitsToExpected_electron_signal", phiCut && dsPlusMCut && DeltaECut && MBCCut);
  bgTree->Draw("kElectron1RatioHitsToExpected>>RatioHitsToExpected_electron_bg", phiCut && dsPlusMCut && DeltaECut && MBCCut);
  
  TH1D *RatioHitsToExpected_positron_signal = new TH1D("RatioHitsToExpected_positron_signal", "RatioHitsToExpected_positron_signal", 10, 0., 1.); RatioHitsToExpected_positron_signal->SetLineColor(kRed);
  TH1D *RatioHitsToExpected_positron_bg = new TH1D("RatioHitsToExpected_positron_bg", "RatioHitsToExpected_positron_bg", 10, 0., 1.); RatioHitsToExpected_positron_bg->SetLineColor(kRed);
  signalTree->Draw("kElectron2RatioHitsToExpected>>RatioHitsToExpected_positron_signal", phiCut && dsPlusMCut && DeltaECut && MBCCut);
  bgTree->Draw("kElectron2RatioHitsToExpected>>RatioHitsToExpected_positron_bg", phiCut && dsPlusMCut && DeltaECut && MBCCut);
  
  TH1D *HitsExpected_electron_signal = new TH1D("HitsExpected_electron_signal", "HitsExpected_electron_signal", 100, 0., 100.);
  TH1D *HitsExpected_electron_bg = new TH1D("HitsExpected_electron_bg", "HitsExpected_electron_bg", 100, 0., 100.);
  signalTree->Draw("kElectron1HitsExpected>>HitsExpected_electron_signal", phiCut && dsPlusMCut && DeltaECut && MBCCut);
  bgTree->Draw("kElectron1HitsExpected>>HitsExpected_electron_bg", phiCut && dsPlusMCut && DeltaECut && MBCCut);
  
  TH1D *HitsExpected_positron_signal = new TH1D("HitsExpected_positron_signal", "HitsExpected_positron_signal", 100, 0., 100.);
  TH1D *HitsExpected_positron_bg = new TH1D("HitsExpected_positron_bg", "HitsExpected_positron_bg", 100, 0., 100.);
  signalTree->Draw("kElectron2HitsExpected>>HitsExpected_positron_signal", phiCut && dsPlusMCut && DeltaECut && MBCCut);
  bgTree->Draw("kElectron2HitsExpected>>HitsExpected_positron_bg", phiCut && dsPlusMCut && DeltaECut && MBCCut);
    
  
  TCanvas *dsPlusM = new TCanvas("dsPlusM");
  dsPlusM->Divide(1,2);
  dsPlusM->cd(1);
  dsPlusM_signal->Draw();
  dsPlusM->cd(2);
  dsPlusM_bg->Draw();
  
  TCanvas *DeltaE = new TCanvas("DeltaE");
  DeltaE->Divide(1,2);
  DeltaE->cd(1);
  DeltaE_signal->Draw();
  DeltaE->cd(2);
  DeltaE_bg->Draw();
  
  TCanvas *MBC = new TCanvas("MBC");
  MBC->Divide(1,2);
  MBC->cd(1);
  MBC_signal->Draw();
  MBC->cd(2);
  MBC_bg->Draw();
  
  TCanvas *DeltaM = new TCanvas("DeltaM");
  DeltaM->Divide(1,2);
  DeltaM->cd(1);
  DeltaM_signal->Draw();
  DeltaM->cd(2);
  DeltaM_bg->Draw();
  
  TCanvas *d0_e = new TCanvas("d0_e");
  d0_e->Divide(1,2);
  d0_e->cd(1);
  d0_e_signal->Draw();
  d0_e->cd(2);
  d0_e_bg->Draw();
  
  TCanvas *d0_p = new TCanvas("d0_p");
  d0_p->Divide(1,2);
  d0_p->cd(1);
  d0_p_signal->Draw();
  d0_p->cd(2);
  d0_p_bg->Draw();
  
  TCanvas *RatioHitsToExpected_electron = new TCanvas("RatioHitsToExpected_electron");
  RatioHitsToExpected_electron->Divide(1,2);
  RatioHitsToExpected_electron->cd(1);
  RatioHitsToExpected_electron_signal->Draw();
  RatioHitsToExpected_electron->cd(2);
  RatioHitsToExpected_electron_bg->Draw();
  
  TCanvas *RatioHitsToExpected_positron = new TCanvas("RatioHitsToExpected_positron");
  RatioHitsToExpected_positron->Divide(1,2);
  RatioHitsToExpected_positron->cd(1);
  RatioHitsToExpected_positron_signal->Draw();
  RatioHitsToExpected_positron->cd(2);
  RatioHitsToExpected_positron_bg->Draw();
  
  TCanvas *HitsExpected_electron = new TCanvas("HitsExpected_electron");
  HitsExpected_electron->Divide(1,2);
  HitsExpected_electron->cd(1);
  HitsExpected_electron_signal->Draw();
  HitsExpected_electron->cd(2);
  HitsExpected_electron_bg->Draw();
  
  TCanvas *HitsExpected_positron = new TCanvas("HitsExpected_positron");
  HitsExpected_positron->Divide(1,2);
  HitsExpected_positron->cd(1);
  HitsExpected_positron_signal->Draw();
  HitsExpected_positron->cd(2);
  HitsExpected_positron_bg->Draw();
  */
  
}


