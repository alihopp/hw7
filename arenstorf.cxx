include <iostream>
#include <fstream>
#include <cmath>

//-------------------
using namespace
//-------------------

// equations: motion of the light

void MOL(double* k, double x0, double x1, double y0, double y1, double mu)
{
  double r = sqrt((x0+mu)*(x0+mu) + y0*y0,2);
  double s = sqrt(pow((x0-1+mu),2) + y0*y0,2);
  
  k[0] = x1;
  k[1] = y1;
  k[2] = x0 + 2.0*y1 - (1-mu)*(x0+mu)/(r*r*r) - mu*(x-1+mu)/(s*s*s);
  k[3] = y0 - 2*x0 - (1-mu)*y0/(r*r*r) - mu*y0/(s*s*s);

  /* k[4] = z2; k[5] = (1-mu)*z/(r*r*r) - mu*z/(s*s*s) 
     not relevant because we are circle in plane (x,y) */
}


//-------------------------
void k_th(double k1, double k2, double k3, double k4, double k5, double k6, double k7, double dt, double mu)
{
  double ; //kofihj

  MOL(k1, x0, x1, y0, y1, mu);

  MOL(k2, x0+dt*1/5*k1[0], x1+dt*1/5*k2[1], y0+dt*1/5*k3[2], y1+dt*1/5*k4[3],mu);

  MOL(k3, x0+dt*(3/40*k1[0]+9/40*k2[0]), x1+dt*(3/40*k1[1]+9/40*k2[1]), y0+dt*(3/40*k1[3]+9/40*k2[3]), x0+dt*(3/40*k1[4]+9/40*k2[4]),mu);

  MOL(k4, x0+dt*(44/45*k1[0]-56/15*k2[0]+32/9*k3[0], x1+dt*(44/45*k1[1]-56/15*k2[1]+32/9*k3[1], y0+dt*(44/45*k1[2]-56/15*k2[2]+32/9*k3[2], y1+dt*(44/45*k1[3]-56/15*k2[3]+32/9*k3[3],mu);

  MOL(k5, x0+dt*(19372/6561*k1[0]-25360/2187*k2[0]+64448/6561*k3[0]-212/729*k4[0]), x1+dt*(19372/6561*k1[1]-25360/2187*k2[1]+64448/6561*k3[1]-212/729*k4[1]), y0+dt*(19372/6561*k1[2]-25360/2187*k2[2]+64448/6561*k3[2]-212/729*k4[2]), y1+dt*(19372/6561*k1[3]-25360/2187*k2[3]+64448/6561*k3[3]-212/729*k4[3]));

  MOL(k6,x0+dt*(9017/3168*k1[0]-355/33*k2[0]+46732/5247*k3[0]+49/176*k4[0]-5103/18656*k5[0]), x1+dt*(9017/3168*k1[1]-355/33*k2[1]+46732/5247*k3[1]+49/176*k4[1]-5103/18656*k5[1]),y0+dt*(9017/3168*k1[2]-355/33*k2[2]+46732/5247*k3[2]+49/176*k4[2]-5103/18656*k5[2]),y1+dt*(9017/3168*k1[3]-355/33*k2[3]+46732/5247*k3[3]+49/176*k4[3]-5103/18656*k5[3]));

  MOL(k7,x0+dt*(35/384*k1[0]+0*k2[0]+500/1113*k3[0]+125/192*k4[0]-2187/6784*k5[0]+11/84*k6[0]), x1+dt*(35/384*k1[1]+0*k2[1]+500/1113*k3[1]+125/192*k4[1]-2187/6784*k5[1]+11/84*k6[1]),y0+dt*(35/384*k1[2]+0*k2[2]+500/1113*k3[2]+125/192*k4[2]-2187/6784*k5[2]+11/84*k6[2]), y1+dt*(35/384*k1[3]+0*k2[3]+500/1113*k3[3]+125/192*k4[3]-2187/6784*k5[3]+11/84*k6[3]));

}

//-------------------------- 
void fifth_k(double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7, double dt, double mu) 
{
 double RK5;

 kth(k1,k2,k3,k4,k5,k6,k7,dt,mu);
 
 for(int i; i<4 ;i++) 
   {
     RK5[i] += dt*(35/384*k1[i] + 500/1113*k3[i] + 125/192*k4[i] - 2187*6784*k5[i] + 11/84*k6[i]);
   }
}

//---------------
void fourth(double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7, double dt, double mu)
{
  double RK4;

  kth(k1,k2,k3,k4,k5,k6,k7,dt,mu);

  for(int i; i<4; i++)
    {
      RK4 += dt*(5149/57600*k1[i] + 7571/16695*k3[i] + 393/640*k4[i] - 92097/339200*k5[i] + 187/2100*k6[i] +1/40*k7[i];
    }

}


//------------------------

void main()
{
  const double mu = 0.012277471;
  double t = 0.0;
  double T = 17.065216560157;

  // Initial condition
  x[0] = 0.994;
  x[1] = 0.0;
  x[2] = 0.0;
  x[3] = -2.00158510637908;

  y[0] = 0.994;
  y[1] = 0.0;
  y[2] = 0.0;
  y[3] = -2.00158510637908;

  ofstream out ("file.txt");

  fifth_k(k1,k2,k3,k4,k5,k6,k7,dt,mu);
  fourth_k(k1,k2,k3,k4,k5,k6,k7,dt,mu);
  
  out.close();
  return 0;
}
