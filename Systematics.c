#include <iostream>

double quad(double a, double b, double c=0, double d=0, double e=0, double f=0)
{
  return pow(a*a+b*b+c*c+d*d+e*e+f*f, 0.5);
}

void Systematics()
{
  double dalitz=1.174e-2;
  double dalitz_err=0.035e-2;
  double gamma=98.823e-2;
  double gamma_err=0.034e-2;
  
  double rat=2*dalitz/gamma;
  double rat_err=rat*quad((dalitz_err/dalitz), (gamma_err/gamma));
  
  double es=1069.0/29974.0;
  double es_err=es/pow(1069.0, 0.5);
  
  double es1=10.0/29974.0;
  double es1_err=es1/pow(10.0, 0.5);
  // es1_err=es1_err/1000;
  
  double ec=2.0/149888.0;
  double ec_err=ec/pow(2.0, 0.5);
  
  double ec1=54.0/149888.0;
  double ec1_err=ec1/pow(54.0, 0.5);
  
  double y=306.0;
  double y_err=pow(y, 0.5);
  double y1=141.0;
  double y1_err=pow(y1, 0.5);
  
  std::cout<<"rat = "<<rat<<" +- "<<rat_err<<std::endl;
  std::cout<<"es = "<<es<<" +- "<<es_err<<std::endl;
  std::cout<<"ec = "<<ec<<" +- "<<ec_err<<std::endl;
  std::cout<<"y = "<<y<<" +- "<<y_err<<std::endl;
  std::cout<<"y1 = "<<y1<<" +- "<<y1_err<<std::endl;
  
  double n1=y/(es+ec/rat);
  double eff=es+ec/rat;
  double eff_err=quad(es_err, (ec/rat)*quad(ec_err/ec, rat_err/rat));
  double n1_err=n1*quad(y_err/y, eff_err/eff);
  
  std::cout<<"n1 (method 1) = "<<n1<<" +- "<<n1_err<<std::endl;
  
  double egamma=25713.0/149888.0;
  double egamma_err=egamma/pow(25713.0, 0.5);
  double ygamma=58602.0;
  double ygamma_err=pow(ygamma, 0.5);
  
  double ngamma=ygamma/egamma;
  double ngamma_err=ngamma*quad(ygamma_err/ygamma, egamma_err/egamma);
  
  std::cout<<"egamma = "<<egamma<<" +- "<<egamma_err<<std::endl;
  std::cout<<"ngamma = "<<ngamma<<" +- "<<ngamma_err<<std::endl;
  
  double rat_m1=n1/ngamma;
  double rat_m1_err=rat_m1*quad(n1_err/n1, ngamma_err/ngamma);
  double b_m1=rat_m1*gamma/2;
  double b_m1_err=b_m1*quad(rat_m1_err/rat_m1, gamma_err/gamma);
  
  std::cout<<"rat_m1 = "<<rat_m1<<" +- "<<rat_m1_err<<std::endl;
  std::cout<<"b_m1 = "<<b_m1<<" +- "<<b_m1_err<<std::endl;
  
  std::cout<<" == method 2 =="<<std::endl;
  
  std::cout<<"es1 = "<<es1<<" +- "<<es1_err<<std::endl;
  std::cout<<"ec1 = "<<ec1<<" +- "<<ec1_err<<std::endl;
  
  double n2=(y*es1-y1*es)/(ec*es1-es*ec1);
  double n2_err=n2*quad(quad(es1*y_err, es1_err*y, es*y1_err, es_err*y1)/(y*es1-y1*es),
                        quad(ec*es1_err, es1*ec_err, ec1*es_err, es*ec1_err)/(ec*es1-ec1*es));
                        
  double n1_m3=(y*ec1-y1*ec)/(es*ec1-ec*es1);
  //double n1_m3_err=n1_m3*quad(quad(ec1*y_err, y*ec1_err, ec*y1_err, y1*ec_err)/(y*ec1-y1*ec),
  //                                quad(ec*es1_err, es1*ec_err, ec1*es_err, es*ec1_err)/(ec1*es-ec*es1));
  
  double n1_m3_err=quad(ec1*y_err, ec*y1_err, 
                        (y/(es*ec1-ec*es1)-(y*ec1-y1*ec)*es/pow(es*ec1-ec*es1, 2))*ec1_err,
                        (-y/(es*ec1-ec*es1)+(y*ec1-y1*ec)*es1/pow(es*ec1-ec*es1, 2))*ec_err,
                        (-(y*ec1-y1*ec)*ec1/pow(es*ec1-ec*es1, 2))*es_err,
                        ((y*ec1-y1*ec)*ec/pow(es*ec1-ec*es1, 2))*es1_err);
  
  std::cout<<"n2 (method 2) = "<<n2<<" +- "<<n2_err<<std::endl;
  
  double n1_m2=n2*rat;
  double n1_m2_err=n1_m2*quad(n2_err/n2, rat_err/rat);
  double rat_m2=n1_m2/ngamma;
  double rat_m2_err=rat_m2*quad(n1_m2_err/n1_m2, ngamma_err/ngamma);
  double b_m2=rat_m2*gamma/2;
  double b_m2_err=b_m2*quad(rat_m2_err/rat_m2, gamma_err/gamma);
  
  std::cout<<"n1 (method 2) = "<<n1_m2<<" +- "<<n1_m2_err<<std::endl;
  std::cout<<"rat (method 2) = "<<rat_m2<<" +- "<<rat_m2_err<<std::endl;
  std::cout<<"b_m2 = "<<b_m2<<" +- "<<b_m2_err<<std::endl;
  
  std::cout<<" == method 3 =="<<std::endl;
  
  double rat_m3=n1_m3/ngamma;
  double rat_m3_err=rat_m3*quad(n1_m3_err/n1_m3, ngamma_err/ngamma);
  double b_m3=rat_m3*gamma/2;
  double b_m3_err=b_m3*quad(rat_m3_err/rat_m3, gamma_err/gamma);
  std::cout<<"n1 (method 3) = "<<n1_m3<<" +- "<<n1_m3_err<<std::endl;
  std::cout<<"rat (method 3) = "<<rat_m3<<" +- "<<rat_m3_err<<std::endl;
  std::cout<<"b_m3 = "<<b_m3<<" +- "<<b_m3_err<<std::endl; 
  
  
}
