#include <iostream>

double pi=3.14159265358979;

double factorial(unsigned int f)
{
  double result=1;
  for (unsigned int i=1; i<=f; ++i) result*=i;
  return result;
}

double quad(double a, double b, double c=0, double d=0, double e=0, double f=0)
{
  return pow(a*a+b*b+c*c+d*d+e*e+f*f, 0.5);
}

double poisson(unsigned int i, double mu)
{
  return pow(mu,i)*exp(-mu)/factorial(i);
}

double upperError(unsigned int i)
{
  double max=log(poisson(i,double(i)));
  double step=pow(double(i), 0.5)/1000.0;
  double diff=0.0;
  double i_step=double(i);
  while(diff<0.5)
  {
    i_step+=step;
    diff=max-log(poisson(i, i_step));
  }
  double i_report=i_step-0.5*step;
  return i_report-i;
}

double lowerError(unsigned int i)
{
  double max=log(poisson(i,double(i)));
  double step=pow(double(i), 0.5)/1000.0;
  double diff=0.0;
  double i_step=double(i);
  while(diff<0.5)
  {
    i_step-=step;
    diff=max-log(poisson(i, i_step));
  }
  double i_report=i_step+0.5*step;
  return i-i_report;
}

double upperErrorZero()
{
  double cumprob=0.0;
  double step=0.0001;
  double i_step=0.0;
  while (cumprob<0.68)
  {
    i_step+=step;
    cumprob+=0.5*(poisson(0, i_step)+poisson(0, i_step-step))*step;
  }
  return i_step;
}

void AsymmetricErrors()
{
  std::cout<<"Error on 10 is + "<<upperError(10)<<" - "<<lowerError(10)<<std::endl;
  
  std::string decay="pietaprime";
  // options: KKpi, KsK, pieta, pietaprime, KKpipi0, pipipi, KsKmpipi, pipi0eta, pietaprimerho
  unsigned int yield=0, bg_N=0;
  double bg=0.0;
  double n_gamma=0.0, n_gamma_stat=0.0;
  double eff_ee=0.0, eff_ee_err=0.0, eff_gamma=0.0, eff_gamma_err=0.0;
  if (decay=="KKpi")
  {
    yield=14;
    bg=1.048;
    bg_N=8;
    n_gamma=9113.72; n_gamma_stat=110.027;
    eff_ee=0.0729; eff_ee_err=0.0019;
    eff_gamma=0.339; eff_gamma_err=0.0015;
  }
  else if (decay=="KsK")
  {
    yield=1;
    bg=0.849;
    bg_N=4;
    n_gamma=1902.08; n_gamma_stat=56.9293;
    eff_ee=0.0597; eff_ee_err=0.0017;
    eff_gamma=0.257268; eff_gamma_err=0.000426529;
  }
  else if (decay=="pieta")
  {
    yield=4;
    bg=1.4;
    bg_N=4;
    n_gamma=1036.66; n_gamma_stat=46.1277;
    eff_ee=0.0855; eff_ee_err=0.0021;
    eff_gamma=0.330947; eff_gamma_err=0.00154357;
  }
  else if (decay=="pietaprime")
  {
    yield=4;
    bg=0.0;
    bg_N=1;
    n_gamma=691.007; n_gamma_stat=33.504;
    eff_gamma=0.210096; eff_gamma_err=0.00128908;
    eff_ee=0.053; eff_ee_err=0.0016;
  }
  else if (decay=="KKpipi0")
  {
    yield=6;
    bg=1.703;
    bg_N=13;
    n_gamma=3592.41; n_gamma_stat=117.642;
    eff_gamma=0.122486; eff_gamma_err=0.000999057;
    eff_ee=0.0255; eff_ee_err=0.0011;
  }
  else if (decay=="pipipi")
  {
    yield=7;
    bg=1.572;
    bg_N=12;
    n_gamma=2744.55; n_gamma_stat=93.0426;
    eff_gamma=0.458331; eff_gamma_err=0.00175936;
    eff_ee=0.0992; eff_ee_err=0.0022;
  }
  else if (decay=="KsKmpipi")
  {
    yield=4;
    bg=1.575;
    bg_N=9;
    n_gamma=1569.61; n_gamma_stat=74.2455;
    eff_gamma=0.191342; eff_gamma_err=0.00115724;
    eff_ee=0.0356; eff_ee_err=0.0013;
  }
  else if (decay=="pipi0eta")
  {
    yield=7;
    bg=2.621;
    bg_N=20;
    n_gamma=3169.7; n_gamma_stat=161.093;
    eff_gamma=0.183935; eff_gamma_err=0.00128925;
    eff_ee=0.0316; eff_ee_err=0.0013;
  }
  else if (decay=="pietaprimerho")
  {
    yield=4;
    bg=1.835;
    bg_N=14;
    n_gamma=1531.29; n_gamma_stat=79.7781;
    eff_gamma=0.317065; eff_gamma_err=0.00153116;
    eff_ee=0.0638; eff_ee_err=0.0018;
  }
  else if (decay=="Sum")
  {
    yield=51;
    bg=12.603;
    bg_N=10;
  }
  else if (decay=="Others")
  {
    yield=37;
    bg=11;
    bg_N=-10;
  }
  
  double yield_upperError=upperError(yield);
  double yield_lowerError=lowerError(yield);
  
  double background_upperError=bg*upperError(bg_N)/double(bg_N);
  double background_lowerError=bg*lowerError(bg_N)/double(bg_N);
  if (decay=="pietaprime")
  {
    background_upperError=upperErrorZero();
    background_lowerError=0;
  }
  
  double signal=yield-bg;
  double signal_upperError=quad(yield_upperError, background_upperError);
  double signal_lowerError=quad(yield_lowerError, background_lowerError);
  
  double eff_ratio=eff_gamma/eff_ee;
  double eff_ratio_err=eff_ratio*quad(eff_gamma_err/eff_gamma, eff_ee_err/eff_ee);
  
  double R=(signal*eff_gamma)/(n_gamma*eff_ee);
  double R_upperError=R*quad(signal_upperError/signal, n_gamma_stat/n_gamma);
  double R_lowerError=R*quad(signal_lowerError/signal, n_gamma_stat/n_gamma);
  
  std::cout<<"Decay "<<decay<<std::endl;
  std::cout<<" ================ "<<std::endl;
  std::cout<<" yield = "<<yield<<" + "<<yield_upperError<<" - "<<yield_lowerError<<" (stat) "<<std::endl;
  std::cout<<" background = "<<bg<<" + "<<background_upperError<<" - "<<background_lowerError<<" (stat) "<<std::endl;
  std::cout<<" signal = "<<signal<<" + "<<signal_upperError<<" - "<<signal_lowerError<<" (stat) "<<std::endl;
  std::cout<<" effRatio = "<<eff_ratio<<" + "<<eff_ratio_err<<std::endl;
  std::cout<<" R = "<<R*100.<<" + "<<R_upperError*100.<<" - "<<R_lowerError*100.<<" (stat) "<<std::endl;
  std::cout<<" ============== "<<std::endl;
  
}
