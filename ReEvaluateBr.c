#include <iostream>

double quad(double a, double b, double c=0, double d=0, double e=0, double f=0)
{
  return pow(a*a+b*b+c*c+d*d+e*e+f*f, 0.5);
}

void ReEvaluateBr()
{

  double m1=0.062;
  double m1_stat=0.005;
  double m1_syst=0.006;
  
  double m2=0.007206620326048;
  // double m2_stat=0.001375398732972;
  double m2_stat=0.0013;
  double m2_syst=0.000950618987392;

  
  double a=1+m1+m2;
  
  double b0=1/a;
  double b0_stat=quad(m1_stat, m2_stat)/(a*a);
  double b0_syst=quad(m1_syst, m2_syst)/(a*a);
  
  double b1=m1/a;
  double db1dm1=1/a-m1/(a*a);
  double db1dm2=-m1/(a*a);
  double b1_stat=quad(db1dm1*m1_stat, db1dm2*m2_stat);
  double b1_syst=quad(db1dm1*m1_syst, db1dm2*m2_syst);
  
  double b2=m2/a;
  double db2dm1=-m2/(a*a);
  double db2dm2=1/a-m2/(a*a);
  double b2_stat=quad(db2dm1*m1_stat, db2dm2*m2_stat);
  double b2_syst=quad(db2dm1*m1_syst, db2dm2*m2_syst);
  
  std::cout<<"B(Ds* -> Ds gamma) = "<<b0<<" +- "<<b0_stat<<" +- "<<b0_syst<<std::endl;
  std::cout<<"B(Ds* -> Ds pi0) = "<<b1<<" +- "<<b1_stat<<" +- "<<b1_syst<<std::endl;
  std::cout<<"B(Ds* -> Ds e+e-) = "<<b2<<" +- "<<b2_stat<<" +- "<<b2_syst<<std::endl;
}
