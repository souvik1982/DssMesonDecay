#include <iostream>
#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"

double yield_mean_Dsstee;
double bg_Dsstee, bg_stat, bg_syst;
double eff_Dsstee;
double yield_Dsstgamma, yield_Dsstgamma_stat;
double eff_Dsstgamma;
double branchingFr_mode, branchingFr_mode_error;
double L=586;
double sigma=948;

void setDecayValues(int d)
{
  if (d==1) // "KKpi"
  {
    yield_mean_Dsstee=14.7;
    bg_Dsstee=1.048;
    bg_stat=0.37;
    bg_syst=0.79;
    eff_Dsstee=0.0729;
    yield_Dsstgamma=9570.19;
    yield_Dsstgamma_stat=157.75;
    eff_Dsstgamma=0.3597;
    branchingFr_mode=0.055;
    branchingFr_mode_error=0.0028;
  }
  else if (d==2) // "KsK"
  {
    yield_mean_Dsstee=3.87;
    bg_Dsstee=0.849;
    bg_stat=0.43;
    bg_syst=0.74;
    eff_Dsstee=0.0597;
    yield_Dsstgamma=1903.29;
    yield_Dsstgamma_stat=66.81;
    eff_Dsstgamma=0.2653;
    branchingFr_mode=0.0149;
    branchingFr_mode_error=0.0009;
  }
  else if (d==3) // "pieta"
  {
    yield_mean_Dsstee=3.21;
    bg_Dsstee=1.4;
    bg_stat=0.7;
    bg_syst=0.49;
    eff_Dsstee=0.0855;
    yield_Dsstgamma=1136.83;
    yield_Dsstgamma_stat=63.76;
    eff_Dsstgamma=0.3448;
    branchingFr_mode=0.0158*0.3931;
    branchingFr_mode_error=branchingFr_mode*pow(pow(0.0021/0.0158, 2)+pow(0.0020/0.3931, 2), 0.5);
  }
  else if (d==4) // "pietaprime"
  {
    yield_mean_Dsstee=1.2;
    bg_Dsstee=0.0;
    bg_stat=0.63;
    bg_syst=0;
    eff_Dsstee=0.053;
    yield_Dsstgamma=711.25;
    yield_Dsstgamma_stat=41.84;
    eff_Dsstgamma=0.221;
    branchingFr_mode=0.038*0.446*0.3931;
    branchingFr_mode_error=branchingFr_mode*pow(pow(0.004/0.038, 2)+pow(0.0014/0.446, 2)+pow(0.0020/0.3931, 2), 0.5);
  }
  else if (d==5) // "KKpipi0"
  {
    yield_mean_Dsstee=6.55;
    bg_Dsstee=1.703;
    bg_stat=0.47;
    bg_syst=0.56;
    eff_Dsstee=0.0255;
    yield_Dsstgamma=4590.13;
    yield_Dsstgamma_stat=204.61;
    eff_Dsstgamma=0.1475;
    branchingFr_mode=0.056;
    branchingFr_mode_error=0.005;
  }
  else if (d==6) // "pipipi"
  {
    yield_mean_Dsstee=5.32;
    bg_Dsstee=1.572;
    bg_stat=0.45;
    bg_syst=0.59;
    eff_Dsstee=0.0992;
    yield_Dsstgamma=2904.97;
    yield_Dsstgamma_stat=113.38;
    eff_Dsstgamma=0.4755;
    branchingFr_mode=0.0111;
    branchingFr_mode_error=0.0008;
  }
  else if (d==7) // "KsKmpipi"
  {
    yield_mean_Dsstee=3.57;
    bg_Dsstee=1.575;
    bg_stat=0.53;
    bg_syst=0.4;
    eff_Dsstee=0.0356;
    yield_Dsstgamma=1788.91;
    yield_Dsstgamma_stat=92.32;
    eff_Dsstgamma=0.2094;
    branchingFr_mode=0.0164;
    branchingFr_mode_error=0.0012;
  }
  else if (d==8) // "pipi0eta"
  {
    yield_mean_Dsstee=8.11;
    bg_Dsstee=2.621;
    bg_stat=0.59;
    bg_syst=0.23;
    eff_Dsstee=0.0316;
    yield_Dsstgamma=1463.87;
    yield_Dsstgamma_stat=103.26;
    eff_Dsstgamma=0.1107;
    branchingFr_mode=0.130*1.*0.3931;
    branchingFr_mode_error=branchingFr_mode*pow(pow(0.022/0.130, 2)+pow(0.0020/0.3931, 2), 0.5);
  }
  else if (d==9) // "pietaprimerho"
  {
    yield_mean_Dsstee=4.26;
    bg_Dsstee=1.835;
    bg_stat=0.49;
    bg_syst=0.25;
    eff_Dsstee=0.0638;
    yield_Dsstgamma=1535.41;
    yield_Dsstgamma_stat=112.98;
    eff_Dsstgamma=0.3348;
    branchingFr_mode=0.038*0.294;
    branchingFr_mode_error=branchingFr_mode*pow(pow(0.004/0.038, 2)+pow(0.009/0.294, 2), 0.5);
  }
  yield_mean_Dsstee*=1;
}


void ToyMC()
{

  TH1D *h_K = new TH1D("h_K", "K (linear addition);", 75, 0.0, 0.015);
  TH1D *h_K_weighted = new TH1D("h_K_weighted", "K (weighted addition);", 75, 0.0, 0.015);
  
  TRandom *rand=new TRandom(0);
  
  for (int i=1; i<=100000; ++i)
  {
    double yield_Dsstee_sum=0;
    double effB_Dsstee_sum=0;
    double yield_Dsstgamma_sum=0;
    double effB_Dsstgamma_sum=0;
    double K_weightedSum=0;
    double weightSum=0;
    for (int d=1; d<=9; ++d)
    //int d=1;
    {
      setDecayValues(d);
      
      double yield_Dsstee=rand->Poisson(yield_mean_Dsstee);
      
      double signal_Dsstee=yield_Dsstee-bg_Dsstee;
      //if (signal_Dsstee<0) signal_Dsstee=0;
      double signal_Dsstee_stat=pow(yield_Dsstee+bg_stat*bg_stat, 0.5);
      yield_Dsstee_sum+=signal_Dsstee;
      effB_Dsstee_sum+=eff_Dsstee*branchingFr_mode;
      yield_Dsstgamma_sum+=yield_Dsstgamma;
      effB_Dsstgamma_sum+=eff_Dsstgamma*branchingFr_mode;
      // double K_mode=(signal_Dsstee/eff_Dsstee)/(yield_Dsstgamma/eff_Dsstgamma);
      // double dK_mode=K_mode*pow(pow(signal_Dsstee_stat/signal_Dsstee, 2)+pow(yield_Dsstgamma_stat/yield_Dsstgamma, 2), 0.5);
      double K_mode=signal_Dsstee/(L*sigma*branchingFr_mode*eff_Dsstee*0.942);
      double dK_mode=K_mode*(signal_Dsstee_stat/pow(pow(signal_Dsstee, 2)+1e-11, 0.5));
      K_weightedSum+=K_mode/(dK_mode*dK_mode);
      weightSum+=1/(dK_mode*dK_mode);
      /*
      std::cout<<"signal_Dsstee = "<<signal_Dsstee<<std::endl;
      std::cout<<"eff_Dsstee = "<<eff_Dsstee<<std::endl;
      std::cout<<"yield_Dsstgamma = "<<yield_Dsstgamma<<std::endl;
      std::cout<<"eff_Dsstgamma = "<<eff_Dsstgamma<<std::endl;
      std::cout<<"mode "<<d<<" K = "<<K_mode<<"+-"<<dK_mode<<std::endl;
      */
    }
     
     // double K=(yield_Dsstee_sum/effB_Dsstee_sum)/(yield_Dsstgamma_sum/effB_Dsstgamma_sum);
     double K=(yield_Dsstee_sum/effB_Dsstee_sum)/(L*sigma*0.942);
     double K_weighted=K_weightedSum/weightSum;
     // std::cout<<"K = "<<K<<std::endl;
     // std::cout<<"K(weighted) = "<<K_weighted<<std::endl;
     h_K->Fill(K);
     h_K_weighted->Fill(K_weighted);
  }
  TLine *line;
  TCanvas *c_K = new TCanvas("c_K", "c_K", 400, 800);
  c_K->Divide(1,2);
  c_K->cd(1);
  h_K->Draw();
  line=new TLine(0.0065, 0, 0.0065, h_K->GetMaximum()); line->Draw();
  c_K->cd(2);
  h_K_weighted->Draw();
  line=new TLine(0.0065, 0, 0.0065, h_K_weighted->GetMaximum()); line->Draw();

}
