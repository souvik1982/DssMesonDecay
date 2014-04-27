void Illustrations()
{
  gROOT->SetStyle("Plain");
  Int_t linsav = gStyle->GetLineWidth();
  //gStyle->SetLineWidth(2);
  
  TLatex t;
  t.SetTextAlign(22);
  t.SetTextSize(0.1);
  
  TCanvas *FeynmanDiagram_Dsstee = new TCanvas("FeynmanDiagram_Dsstee", "FeynmanDiagram_Dsstee", 500, 500);
  FeynmanDiagram_Dsstee->Range(0, 0, 100, 100);
  
  TArrow *line;
  gStyle->SetLineWidth(3);
  line=new TArrow(5, 40, 40, 40, 0.03, "->-"); line->Draw();
  line=new TArrow(40, 40, 70, 5, 0.03, "->-"); line->Draw();
  gStyle->SetLineWidth(0);
  line=new TArrow(60, 60, 60, 95, 0.02, "->-"); line->Draw();
  line=new TArrow(60, 60, 95, 60, 0.02, "->-"); line->Draw();
  
  TCurlyLine *gamma;
  gamma=new TCurlyLine(40, 40, 60, 60); gamma->SetWavy(); gamma->Draw();
  
  t.DrawLatex(20, 50, "D_{s}^{*+}");
  t.DrawLatex(60, 30, "D_{s}^{+}");
  t.DrawLatex(45, 55, "#gamma");
  t.DrawLatex(65, 80, "e^{-}");
  t.DrawLatex(80, 67, "e^{+}");
  
  TCanvas *FeynmanDiagram_Dsstgamma = new TCanvas("FeynmanDiagram_Dsstgamma", "FeynmanDiagram_Dsstgamma", 500, 500);
  FeynmanDiagram_Dsstgamma->Range(0, 0, 100, 100);
  
  TArrow *line;
  gStyle->SetLineWidth(3);
  line=new TArrow(5, 50, 50, 50, 0.03, "->-"); line->Draw();
  line=new TArrow(50, 50, 90, 10, 0.03, "->-"); line->Draw();
  
  TCurlyLine *gamma;
  gStyle->SetLineWidth(0);
  gamma=new TCurlyLine(50, 50, 90, 90); gamma->SetWavy(); gamma->Draw();
  
  t.DrawLatex(20, 60, "D_{s}^{*+}");
  t.DrawLatex(70, 20, "D_{s}^{+}");
  t.DrawLatex(70, 80, "#gamma");
  
  // Diff d_0
  TCanvas *DeltaD0 = new TCanvas("DeltaD0", "DeltaD0", 500, 500);
  DeltaD0->Range(0, 0, 100, 100);
  
  TLine *axis;
  axis=new TLine(30, 0, 30, 100); axis->Draw();
  axis=new TLine(0, 30, 100, 30); axis->Draw();
  
  TEllipse *curv;
  // gStyle->SetLineColor(kBlue);
  curv=new TEllipse(60, 0, 42); curv->Draw();
  curv=new TEllipse(0, 60, 42); curv->Draw();
  // gStyle->SetLineColor(kRed);
  curv=new TEllipse(20, 100, 56); curv->Draw();
  curv=new TEllipse(100, 20, 56); curv->Draw();
  
  // gStyle->SetLineColor(134);
  axis=new TLine(30, 30, 20, 100); axis->Draw();
  axis=new TLine(30, 30, 100, 20); axis->Draw();
  line=new TArrow(30, 30, 28, 44, 0.03, "<>"); line->Draw();
  line=new TArrow(30, 30, 44, 28, 0.03, "<>"); line->Draw();
  
  TMarker *point;
  point=new TMarker(60, 60, 20); point->Draw();
  
  t.SetTextSize(0.05);
  // t.SetTextColor(kBlue);
  t.DrawLatex(10, 25, "signal e^{-}");
  t.DrawLatex(30, 10, "signal e^{+}");
  t.DrawLatex(35, 23, "#Deltad_{0}~0");
  // t.SetTextColor(kRed);
  t.DrawLatex(55, 80, "conversion e^{-}");
  t.DrawLatex(80, 60, "conversion e^{+}");
  t.DrawLatex(30, 40, "d_{0}^{e^{-}}");
  t.DrawLatex(40, 30, "d_{0}^{e^{+}}");
  
  // Diff phi_0
  TCanvas *DeltaPhi0 = new TCanvas("DeltaPhi0", "DeltaPhi0", 500, 500);
  DeltaPhi0->Range(0, 0, 100, 100);
  
  TLine *axis;
  axis=new TLine(30, 0, 30, 100); axis->Draw();
  axis=new TLine(0, 30, 100, 30); axis->Draw();
  
  TEllipse *curv;
  // gStyle->SetLineColor(kBlue);
  curv=new TEllipse(60, 0, 42); curv->Draw();
  curv=new TEllipse(0, 60, 42); curv->Draw();
  // gStyle->SetLineColor(kRed);
  curv=new TEllipse(20, 100, 56); curv->Draw();
  curv=new TEllipse(100, 20, 56); curv->Draw();
  
  // gStyle->SetLineColor(134);
  axis=new TLine(30, 30, 20, 100); axis->Draw();
  axis=new TLine(30, 30, 100, 20); axis->Draw();
  axis=new TLine(28, 44, 3.25, 40.47); axis->Draw();
  axis=new TLine(44, 28, 40.47, 3.25); axis->Draw();
  axis=new TLine(3.25, 41.18, 25, 41.18); axis->Draw();
  axis=new TLine(37, 5, 55, 5); axis->Draw();
  
  TArc *arc;
  arc=new TArc(8.20, 41.18, 10, 0, 8.13); arc->Draw();
  arc=new TArc(40.5, 5, 5, 0, 81.9); arc->Draw("only");
  
  TMarker *point;
  point=new TMarker(60, 60, 20); point->Draw();
  
  t.SetTextSize(0.05);
  // t.SetTextColor(kBlue);
  t.DrawLatex(10, 25, "signal e^{-}");
  t.DrawLatex(30, 10, "signal e^{+}");
  t.DrawLatex(35, 23, "#Delta#phi_{0}~0");
  // t.SetTextColor(kRed);
  t.DrawLatex(55, 80, "conversion e^{-}");
  t.DrawLatex(80, 60, "conversion e^{+}");
  t.DrawLatex(12, 37, "#phi_{0}^{e^{-}}");
  t.DrawLatex(49, 9, "#phi_{0}^{e^{+}}");
  
  // e+e- Collision
  
  TCanvas *Collision = new TCanvas("Collision", "Collision", 400, 300);
  Collision->Range(0, 0, 400, 300);
  
  gStyle->SetLineWidth(1);
  line=new TArrow(5, 200, 200, 200, 0.03, "->-"); line->Draw();
  line=new TArrow(395, 200, 200, 200, 0.03, "->-"); line->Draw();
  gStyle->SetLineWidth(3);
  line=new TArrow(200, 200, 150, 150, 0.03, "->-"); line->Draw();
  gStyle->SetLineWidth(2);
  line=new TArrow(200, 200, 275, 275, 0.03, "->-"); line->Draw();
  line=new TArrow(150, 150, 275, 25, 0.03, "->-"); line->Draw();
  gStyle->SetLineWidth(1);
  line=new TArrow(150, 150, 25, 25, 0.03, "->-"); line->Draw();
  line=new TArrow(150, 150, 25, 150, 0.03, "->-"); line->Draw();
  
  t.SetTextSize(0.05);
  t.SetTextColor(kBlack);
  t.DrawLatex(100, 220, "e^{+}");
  t.DrawLatex(300, 220, "e^{-}");
  t.DrawLatex(270, 250, "D_{s}^{-}");
  t.DrawLatex(200, 175, "D_{s}^{*+}");
  t.DrawLatex(240, 80, "D_{s}^{+}");
  t.DrawLatex(100, 165, "e^{-}");
  t.DrawLatex(120, 100, "e^{+}");
  
}
