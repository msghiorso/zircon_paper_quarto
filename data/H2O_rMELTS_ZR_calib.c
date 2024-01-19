#ifdef DEBUG
#undef DEBUG
#endif

#ifdef ABS
#undef ABS
#endif
#define ABS(x)   ((x) < 0 ? -(x) : (x))
#ifdef MAX
#undef MAX
#endif
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((a) < (b) ? (a) : (b))

static const double tr = 298.15;
static const double pr = 1.0;
static const double trl = 1673.0;
static const double r = 8.3143;
static char *identifier = "Wed Jan  6 09:29:48 2021";

#define SQUARE(x)  ((x)*(x))
#define CUBICE(x)  ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))
#define QUINTIC(x) ((x)*(x)*(x)*(x)*(x))

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>

static double psat2(double t)
{
   double w, wsq, v, ff;
   int i;

   static const double a[9] = { 0.0,
      -7.8889166, 2.5514255, -6.716169, 33.239495, 
      -105.38479, 174.35319, -148.39348, 48.631602
   };
 
   if (t <= 314.0) { 
      return exp(6.3573118-8858.843/t + 607.56335/pow(t, (double) 0.6));
   } else { 
      v = t/647.25;
      w = ABS(1.0-v);
      wsq = sqrt(w);
      ff = 0.0;
      for(i=1;i<=8;i++) {
         ff = ff + a[i]*w;
         w = w*wsq;
      }
      return 220.93*exp(ff/v);
   }
}

static void kubik(double b, double c, double d, 
   double *x1, double *x2, double *x2i, double *x3)
{
   double q, p, r, phi3, ff;

   static const double pi = 3.14159263538979;

   *x2 = 0.0;
   *x2i = 0.0;
   *x3 = 0.0;

   if (c == 0.0 && d == 0.0) { 
      *x1 = -b;
      return; 
   }

   q = (2.0*CUBICE(b)/27.0 - b*c/3.0 + d)/2.0;
   p = (3.0*c - b*b)/9.0;
   ff = ABS(p);
   r = sqrt(ff);
   ff = r*q;
   if (ff < 0.0) r = -r;
   ff = q/CUBICE(r);
 
   if (p > 0.0) { 
      phi3 = log(ff + sqrt(ff*ff+1.0))/3.0;
      *x1 = -r*(exp(phi3) - exp(-phi3)) - b/3.0;
      *x2i = 1.0;
   } else { 
      if (q*q + p*p*p > 0.0) { 
         phi3 = log(ff + sqrt(ff*ff-1.0))/3.0;
         *x1 = -r*(exp(phi3) + exp(-phi3)) - b/3.0;
         *x2i = 1.0;
      } else { 
         phi3 = atan(sqrt(1.0-ff*ff)/ff)/3.0;
         *x1 = -2.0*r*cos(phi3) - b/3.0;
         *x2 = 2.0*r*cos(pi/3.0-phi3) - b/3.0;
         *x2i = 0.0;
         *x3 = 2.0*r*cos(pi/3.0+phi3) - b/3.0;
      } 
   } 
   return; 
}

static void whaar(double p, double t, double *gH2O, double *hH2O, double *sH2O,
  double *cpH2O, double *dcpdtH2O, double *vH2O, double *dvdtH2O, 
  double *dvdpH2O, double *d2vdt2H2O, double *d2vdtdpH2O, double *d2vdp2H2O)
{
   /* Calculates and returns the thermodynamic properties of Water and 
      using the Haar equation of State. This routine uses a Redlich-Kwong 
      equation to obtain a first guess for the density of water, 
      thus speeding things up                                                */ 
 
   /* gi are in (BAR CC / G)  =  10 * (j / G)                                */ 
 
   static const double gi[41] = { 0.0,
      -.53062968529023e4,  .22744901424408e5,  .78779333020687e4,
      -.69830527374994e3,  .17863832875422e6, -.39514731563338e6, 
       .33803884280753e6, -.13855050202703e6, -.25637436613260e7, 
       .48212575981415e7, -.34183016969660e7,  .12223156417448e7, 
       .11797433655832e8, -.21734810110373e8,  .10829952168620e8, 
      -.25441998064049e7, -.31377774947767e8,  .52911910757704e8, 
      -.13802577177877e8, -.25109914369001e7,  .46561826115608e8, 
      -.72752773275387e8,  .41774246148294e7,  .14016358244614e8, 
      -.31555231392127e8,  .47929666384584e8,  .40912664781209e7, 
      -.13626369388386e8,  .69625220862664e7, -.10834900096447e8, 
      -.22722827401688e7,  .38365486000660e7,  .68833257944332e5, 
       .21757245522644e6, -.26627944829770e5, -.70730418082074e6, 
      -.225e1,            -1.68e1,             .055e1, 
      -93.0e1 
   };
 
   static const int ki[41] = { 0,
      1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 
      6, 6, 6, 6, 7, 7, 7, 7, 9, 9, 9, 9, 3, 3, 1, 5, 2, 2, 2, 4
   };
      
   static const int li[41] = { 0,
      1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 
      1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 0, 3, 3, 3, 0, 2, 0, 0
   };
 
   static const double ci[19] = { 0.0,
      .19730271018e2,      .209662681977e2,   -.483429455355, 
      .605743189245e1,   22.56023885,        -9.87532442, 
     -.43135538513e1,      .458155781,        -.47754901883e-1, 
      .41238460633e-2,    -.27929052852e-3,    .14481695261e-4, 
     -.56473658748e-6,     .16200446e-7,      -.3303822796e-9, 
      .451916067368e-11,  -.370734122708e-13,  .137546068238e-15
   };
 
   static const double rhoi[41] = { 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.319, 0.310, 0.310, 1.55
   };
 
   static const double ttti[41] = { 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 640.0, 640.0, 641.6, 270.0
   };
      
   static const double alpi[41] = { 0.0, 
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,    0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,    0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,    0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 34.0, 40.0, 30.0, 1050.0
   };
 
   static const double beti[41] = { 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,    0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,    0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,    0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0e4, 2.0e4, 4.0e4, 25.0
   }; 

   static const double bi[6]  = {
      0.7478629,  -0.3540782,    0.0,      0.007159876, 
      0.0,        -0.003528426
   };
   
   static const double bbi[6] = { 
      1.1278334,  -0.5944001,   -5.010996, 0.0,
      0.63684256,  0.0 
   };

   static const double r     = 4.6152;
   static const double gref  = -54955.23970014; /* -54955.2356146121147 25 c and 1 bar */
   static const double href  = -34099.89230644;
   /* 
   static const double sref  =  69.94917790942; 
   */
   static const double t0    = 647.073;
   static const double rr    = 8.31441;
   static const double alpha = 11.0;
   static const double beta  = 133.0/3.0;
   static const double gamma = 7.0/2.0;
   static const double P0    = 1.01325;
   
   double ps = 220.55;
   double taui[7], ermi[10], ark, brk, oft, buk, cuk, duk, x1, x2, x2i, x3, 
     vol, rhn, dp, dr, rh, pr, dpr, q10, qm, x, tr, dtrdt, dpdrh, dpdt, drhdt,
     d2pdt2, d2pdtdrh, d2pdrh2, temp;
   double b, dbdt, d2bdt2, d3bdt3;
   double bb, dbbdt, d2bbdt2, d3bbdt3; 
   double y, dydt, dydrh, d2ydt2, d2ydtdrh, d2ydrh2, d3ydt3, d3ydt2drh, 
          d3ydtdrh2, d3ydrh3;
   double Z, dZdt, dZdrh, d2Zdt2, d2Zdtdrh, d2Zdrh2, d3Zdt3, d3Zdt2drh, 
          d3Zdtdrh2, d3Zdrh3;
   double Ab, dAbdt, dAbdrh, d2Abdt2, d2Abdtdrh, d2Abdrh2, d3Abdt3, d3Abdt2drh, 
          d3Abdtdrh2, d3Abdrh3;
   double Ar, dArdt, dArdrh, d2Ardt2, d2Ardtdrh, d2Ardrh2, d3Ardt3, d3Ardt2drh, 
          d3Ardtdrh2, d3Ardrh3;
   double Ai, dAidt, d2Aidt2, d3Aidt3;
   double A, dAdt, dAdrh, d2Adt2, d2Adtdrh, d2Adrh2, d3Adt3, d3Adt2drh, 
          d3Adtdrh2, d3Adrh3;
   int i, count;

 
   /* The values (T/T0)**i are stored in the array TAUI[i] */
   taui[0] = 1.0; taui[1] = t/t0; for(i=2;i<=6;i++) taui[i] = taui[i-1]*taui[1];

   b  =  bi[1]*log(taui[1]) 
      + bi[0]  + bi[3]/taui[3]  + bi[5]/taui[5];
   bb = bbi[0] + bbi[1]/taui[1] + bbi[2]/taui[2] + bbi[4]/taui[4];
 
   if (t <= 647.25) ps = psat2(t);

   /**************************************************************************
    set initial guess for rho using thb-fit to redlich-kwong
    **************************************************************************/ 

   ark = 1.279186e8 - 2.241415e4 * t;
   brk = 1.428062e1 + 6.092237e-4 * t;

   oft = ark/(p*sqrt(t));
   buk = -10.0*rr*t/p;
   cuk = oft - brk*brk + brk*buk;
   duk = - brk*oft;
   kubik(buk, cuk, duk, &x1, &x2, &x2i, &x3);

   if (x2i != 0.0) vol = x1;
   else            vol = (p < ps) ? MAX(x1, MAX(x2, x3)) : MIN(x1, MIN(x2, x3));

   rhn = (vol <= 0.0) ? 1.9 : (1.0/vol)*18.0152;

   /**************************************************************************/

   dp = DBL_MAX; dr = DBL_MAX;
   for (count=1; count<=100 && (dp > 10.0*DBL_EPSILON || 
                                dr > 10.0*DBL_EPSILON); count++) {
      rh = rhn;
      if (rh <= 0.0) rh = 1.e-8;
      if (rh >  1.9) rh = 1.9;
      y = rh*b/4.0;
      ermi[0] = 1.0; ermi[1] = 1.0-exp(-rh); 
        for (i=2;i<=9;i++) ermi[i] = ermi[i-1]*ermi[1];

      pr = 0.0; dpr = 0.0;

      for (i=1; i<=36; i++) {
        pr  = pr + gi[i]/taui[li[i]]*ermi[ki[i]-1];
        dpr = dpr + (2.0+rh*(ki[i]*exp(-rh)-1.0)/ermi[1])*gi[i]/taui[li[i]]
                     *ermi[ki[i]-1];
      } 
      for(i=37; i<=40; i++) {
        double del, tau, abc;
        del = rh/rhoi[i] - 1.0;
        tau = t/ttti[i]  - 1.0;
        abc = -alpi[i] * pow(del, (double) ki[i]) - beti[i] * tau*tau;
        q10 = (abc > - 100.00) ? gi[i] * pow(del, (double) li[i]) * exp(abc)
                               : 0.0;
        qm = li[i]/del - ki[i]*alpi[i]*pow(del, (double) (ki[i]-1));
        pr = pr + q10*qm*rh*rh/rhoi[i];
        dpr = dpr + (q10*qm*rh*rh/rhoi[i]) * (2.0/rh+qm/rhoi[i]) 
                  - rh*rh/(rhoi[i]*rhoi[i])*q10* 
                    (li[i]/del/del + ki[i]*(ki[i]-1)*alpi[i]
                    *pow(del, (double) (ki[i]-2)));
      }

      pr = rh*(rh*exp(-rh)*pr + r*t*((1.0 + alpha*y + beta*y*y)/CUBICE(1.0-y)
                              + 4.0*y*(bb/b - gamma)));
      dpr = rh*exp(-rh)*dpr 
          + r*t*( (1.0 + 2.0*alpha*y + 3.0*beta*y*y)/CUBICE(1.0-y) 
                 + 3.0*y*(1.0 + alpha*y + beta*y*y)/QUARTIC(1.0-y) 
                 + 2.0*4.0*y*(bb/b - gamma));

      if (dpr <= 0.0) rhn *= (p <= ps) ? 0.95 : 1.05;
      else { 
        if (dpr < 0.01) dpr = 0.01;
        x = (p - pr)/dpr;
        if (ABS(x) > 0.1) x = 0.1*x/ABS(x);
        rhn = rh + x;
      } 
      dp = ABS(1.0 - pr/p);
      dr = ABS(1.0 - rhn/rh);
   }
#ifdef DEBUG
   printf("count, rh, dp, dr: %d %14.6g %14.6g %14.6g\n", count, rh, dp, dr);
#endif
   rh = rhn;

   dbdt    = bi[1]/t - 3.0*bi[3]/(taui[3]*t)  - 5.0*bi[5]/(taui[5]*t);
   d2bdt2  = -bi[1]/(t*t) + 4.0*3.0*bi[3]/(taui[3]*t*t)  
           + 6.0*5.0*bi[5]/(taui[5]*t*t);
   d3bdt3  = 2.0*bi[1]/(t*t*t) - 5.0*4.0*3.0*bi[3]/(taui[3]*t*t*t)  
           - 7.0*6.0*5.0*bi[5]/(taui[5]*t*t*t);

   dbbdt   = - bbi[1]/(taui[1]*t) - 2.0*bbi[2]/(taui[2]*t) 
           - 4.0*bbi[4]/(taui[4]*t);
   d2bbdt2 = 2.0*bbi[1]/(taui[1]*t*t) + 3.0*2.0*bbi[2]/(taui[2]*t*t) 
           + 5.0*4.0*bbi[4]/(taui[4]*t*t);
   d3bbdt3 = - 3.0*2.0*bbi[1]/(taui[1]*t*t*t) 
           - 4.0*3.0*2.0*bbi[2]/(taui[2]*t*t*t) 
           - 6.0*5.0*4.0*bbi[4]/(taui[4]*t*t*t);

   y         = rh*b/4.0;
   dydt      = rh*dbdt/4.0;
   dydrh     = b/4.0;
   d2ydt2    = rh*d2bdt2/4.0;
   d2ydtdrh  = dbdt/4.0;
   d2ydrh2   = 0.0;
   d3ydt3    = rh*d3bdt3/4.0;
   d3ydt2drh = d2bdt2/4.0;
   d3ydtdrh2 = 0.0;
   d3ydrh3   = 0.0;

   ermi[0] = 1.0; ermi[1] = 1.0-exp(-rh);
     for (i=2; i<=9; i++) ermi[i] = ermi[i-1]*ermi[1];

   /* calculate base function                                                */ 
 
   Z         = - log(1.0-y) - (beta-1.0)/(1.0-y) 
             + (alpha+beta+1.0)/(2.0*(1.0-y)*(1.0-y)) + 4*y*(bb/b - gamma) 
             - (alpha-beta+3.0)/2.0 + log(rh*r*t/P0);
   dZdt      = 1.0/t + 4.0*(b*dbbdt-bb*dbdt)*y/(b*b) + 4.0*(bb/b-gamma)*dydt
             + dydt/(1.0-y) + (alpha+beta+1.0)*dydt/CUBICE(1.0-y)
             - (beta-1.0)*dydt/SQUARE(1.0-y);
   dZdrh     = 1.0/rh + 4.0*(bb/b-gamma)*dydrh + dydrh/(1.0-y)
             + (alpha+beta+1.0)*dydrh/CUBICE(1.0-y) 
             - (beta-1.0)*dydrh/SQUARE(1.0-y);
   d2Zdt2    = -1.0/(t*t) + 4.0*(((b*d2bbdt2-bb*d2bdt2)*(b*b)
                 -2.0*(b*dbbdt-bb*dbdt)*b*dbdt)*y)/QUARTIC(b)
             + SQUARE(dydt/(1.0-y)) 
             + 3.0*(alpha+beta+1.0)*SQUARE(dydt)/QUARTIC(1.0-y)
             - 2.0*(beta-1.0)*SQUARE(dydt)/CUBICE(1.0-y)
             + 8.0*(b*dbbdt-bb*dbdt)*dydt/(b*b) + 4.0*(bb/b-gamma)*d2ydt2
             + d2ydt2/(1.0-y) + (alpha+beta+1.0)*d2ydt2/CUBICE(1.0-y)
             - (beta-1.0)*d2ydt2/SQUARE(1.0-y);
   d2Zdtdrh  = 4.0*(bb/b-gamma)*d2ydtdrh + d2ydtdrh/(1.0-y)
             + (alpha+beta+1.0)*d2ydtdrh/CUBICE(1.0-y)
             - (beta-1.0)*d2ydtdrh/SQUARE(1.0-y)
             + 4.0*(b*dbbdt-bb*dbdt)*dydrh/(b*b) + dydt*dydrh/SQUARE(1.0-y)
             + 3.0*(alpha+beta+1.0)*dydt*dydrh/QUARTIC(1.0-y)
             - 2.0*(beta-1.0)*dydt*dydrh/CUBICE(1.0-y);
   d2Zdrh2   = -1.0/(rh*rh) + SQUARE(dydrh/(1.0-y)) 
             + 3.0*(alpha+beta+1.0)*SQUARE(dydrh)/QUARTIC(1.0-y)
             - 2.0*(beta-1.0)*SQUARE(dydrh)/CUBICE(1.0-y)
             + 4.0*(bb/b-gamma)*d2ydrh2 + d2ydrh2/(1.0-y)
             + (alpha+beta+1.0)*d2ydrh2/CUBICE(1.0-y)
             - (beta-1.0)*d2ydrh2/SQUARE(1.0-y);
   d3Zdt3    = 2.0/CUBICE(t) + 4.0*(((-2.0*(b*dbbdt-bb*dbdt)*SQUARE(dbdt)
                 + (b*d3bbdt3+d2bbdt2*dbdt-dbbdt*d2bdt2-bb*d3bdt3)*SQUARE(b)
                 - 2.0*(b*dbbdt-bb*dbdt)*b*d2bdt2)*y 
                 + ((b*d2bbdt2-bb*d2bdt2)*SQUARE(b) 
                 - 2.0*(b*dbbdt-bb*dbdt)*b*dbdt)*dydt)*QUARTIC(b)
                 - 4.0*((b*d2bbdt2-bb*d2bdt2)*SQUARE(b)
                 -2.0*(b*dbbdt-bb*dbdt)*b*dbdt)*CUBICE(b)*y*dbdt)
                 /(QUARTIC(b)*QUARTIC(b))
             + 2.0*CUBICE(dydt/(1.0-y)) 
             + 12.0*(alpha+beta+1.0)*CUBICE(dydt)/QUINTIC(1.0-y)
             - 6.0*(beta-1.0)*CUBICE(dydt)/QUARTIC(1.0-y)
             + 8.0*((b*d2bbdt2-bb*d2bdt2)*SQUARE(b)
                 -2.0*(b*dbbdt-bb*dbdt)*b*dbdt)*dydt/QUARTIC(b)
             + 12.0*(b*dbbdt-bb*dbdt)*d2ydt2/SQUARE(b)
             + 3.0*dydt*d2ydt2/SQUARE(1.0-y)
             + 9.0*(alpha+beta+1.0)*dydt*d2ydt2/QUARTIC(1.0-y)
             - 6.0*(beta-1.0)*dydt*d2ydt2/CUBICE(1.0-y)
             + 4.0*(bb/b-gamma)*d3ydt3 + d3ydt3/(1.0-y)
             + (alpha+beta+1.0)*d3ydt3/CUBICE(1.0-y)
             - (beta-1.0)*d3ydt3/SQUARE(1.0-y);
   d3Zdt2drh = 4.0*(((b*d2bbdt2-bb*d2bdt2)*SQUARE(b)
                 -2.0*(b*dbbdt-bb*dbdt)*b*dbdt)*dydrh)/QUARTIC(b)
             + 4.0*(bb/b-gamma)*d3ydt2drh + d3ydt2drh/(1.0-y)
             + (alpha+beta+1.0)*d3ydt2drh/CUBICE(1.0-y)
             - (beta-1.0)*d3ydt2drh/SQUARE(1.0-y)
             + 8.0*(b*dbbdt-bb*dbdt)*d2ydtdrh/SQUARE(b)
             + 2.0*dydt*d2ydtdrh/SQUARE(1.0-y)
             + 6.0*(alpha+beta+1.0)*dydt*d2ydtdrh/QUARTIC(1.0-y)
             - 4.0*(beta-1.0)*dydt*d2ydtdrh/CUBICE(1.0-y)
             + 2.0*SQUARE(dydt)*dydrh/CUBICE(1.0-y)
             + 12.0*(alpha+beta+1.0)*SQUARE(dydt)*dydrh/QUINTIC(1.0-y)
             - 6.0*(beta-1.0)*SQUARE(dydt)*dydrh/QUARTIC(1.0-y)
             + d2ydt2*dydrh/SQUARE(1.0-y)
             + 3.0*(alpha+beta+1.0)*d2ydt2*dydrh/QUARTIC(1.0-y)
             - 2.0*(beta-1.0)*d2ydt2*dydrh/CUBICE(1.0-y);
   d3Zdtdrh2 = 2.0*dydt*SQUARE(dydrh)/CUBICE(1.0-y)
             + 12.0*(alpha+beta+1.0)*dydt*SQUARE(dydrh)/QUINTIC(1.0-y)
             - 6.0*(beta-1.0)*dydt*SQUARE(dydrh)/QUARTIC(1.0-y)
             + 4.0*(bb/b-gamma)*d3ydtdrh2 + d3ydtdrh2/(1.0-y)
             + (alpha+beta+1.0)*d3ydtdrh2/CUBICE(1.0-y)
             - (beta-1.0)*d3ydtdrh2/SQUARE(1.0-y)
             + 2.0*d2ydtdrh*dydrh/SQUARE(1.0-y)
             + 6.0*(alpha+beta+1.0)*d2ydtdrh*dydrh/QUARTIC(1.0-y)
             - 4.0*(beta-1.0)*d2ydtdrh*dydrh/CUBICE(1.0-y)
             + 4.0*(b*dbbdt-bb*dbdt)*d2ydrh2/SQUARE(b)
             + dydt*d2ydrh2/SQUARE(1.0-y)
             + 3.0*(alpha+beta+1.0)*dydt*d2ydrh2/QUARTIC(1.0-y)
             - 2.0*(beta-1.0)*dydt*d2ydrh2/CUBICE(1.0-y);
   d3Zdrh3   = 2.0/CUBICE(rh) + 2.0*CUBICE(dydrh/(1.0-y))
             + 12.0*(alpha+beta+1.0)*CUBICE(dydrh)/QUINTIC(1.0-y)
             - 6.0*(beta-1.0)*CUBICE(dydrh)/QUARTIC(1.0-y)
             + 3.0*dydrh*d2ydrh2/SQUARE(1.0-y)
             + 9.0*(alpha+beta+1.0)*dydrh*d2ydrh2/QUARTIC(1.0-y)
             - 6.0*(beta-1.0)*dydrh*d2ydrh2/CUBICE(1.0-y)
             + 4.0*(bb/b-gamma)*d3ydrh3 + d3ydrh3/(1.0-y)
             + (alpha+beta+1.0)*d3ydrh3/CUBICE(1.0-y)
             - (beta-1.0)*d3ydrh3/SQUARE(1.0-y);

   Ab         = r*t*Z;
   dAbdt      = r*Z + r*t*dZdt;
   dAbdrh     = r*t*dZdrh;
   d2Abdt2    = 2.0*r*dZdt + r*t*d2Zdt2;
   d2Abdtdrh  = r*dZdrh + r*t*d2Zdtdrh;
   d2Abdrh2   = r*t*d2Zdrh2;
   d3Abdt3    = 3.0*r*d2Zdt2 + r*t*d3Zdt3;
   d3Abdt2drh = 2.0*r*d2Zdtdrh + r*t*d3Zdt2drh;
   d3Abdtdrh2 = r*d2Zdrh2 + r*t*d3Zdtdrh2;
   d3Abdrh3   = r*t*d3Zdrh3;

   /* calculate residual function                                            */
 
   for (i=1, Ar=0.0, dArdt=0.0, dArdrh=0.0, d2Ardt2=0.0, d2Ardtdrh=0.0, 
        d2Ardrh2=0.0, d3Ardt3=0.0, d3Ardt2drh=0.0, d3Ardtdrh2=0.0, 
        d3Ardrh3=0.0; i<=36; i++) {
     Ar         += gi[i]/ki[i]/taui[li[i]]*ermi[ki[i]];
     dArdt      += -li[i]*gi[i]/ki[i]/(taui[li[i]]*t)*ermi[ki[i]];
     dArdrh     += gi[i]/taui[li[i]]*ermi[ki[i]-1]*exp(-rh);
     d2Ardt2    += (li[i]+1.0)*li[i]*gi[i]/ki[i]/(taui[li[i]]*t*t)*ermi[ki[i]];
     d2Ardtdrh  += -li[i]*gi[i]/(taui[li[i]]*t)*ermi[ki[i]-1]*exp(-rh);
     d2Ardrh2   += (ki[i] > 1) ? gi[i]/taui[li[i]]
                   *((ki[i]-1.0)*ermi[ki[i]-2]*exp(-rh)-ermi[ki[i]-1])*exp(-rh)
                     : -gi[i]/taui[li[i]]*exp(-rh);
     d3Ardt3    += -(li[i]+2.0)*(li[i]+1.0)*li[i]*gi[i]/ki[i]
                     /(taui[li[i]]*t*t*t)*ermi[ki[i]];
     d3Ardt2drh += (li[i]+1.0)*li[i]*gi[i]/(taui[li[i]]*t*t)
                   *ermi[ki[i]-1]*exp(-rh);
     d3Ardtdrh2 += (ki[i] > 1) ? -li[i]*gi[i]/(taui[li[i]]*t)
                   *((ki[i]-1.0)*ermi[ki[i]-2]*exp(-rh)-ermi[ki[i]-1])*exp(-rh)
                     : li[i]*gi[i]/(taui[li[i]]*t)*exp(-rh);
     if (ki[i] > 2) d3Ardrh3 += gi[i]/taui[li[i]]
                   *(((ki[i]-2.0)*ermi[ki[i]-3]*exp(-rh) - 3.0*ermi[ki[i]-2])
                   *(ki[i]-1.0)*exp(-rh)+ermi[ki[i]-1])*exp(-rh);
     else d3Ardrh3 += (ki[i] > 1) ? -gi[i]/taui[li[i]]*
                                     (4.0*exp(-rh)-1.0)*exp(-rh)
                                  : gi[i]/taui[li[i]]*exp(-rh);
   }
   for(i=37; i<=40; i++) {
     double del = rh/rhoi[i] - 1.0;
     double tau = t/ttti[i] - 1.0;
     double Q         = -alpi[i]*pow(del, (double) ki[i]) - beti[i]*tau*tau;
     double dQdt      = -beti[i]*2.0*tau/ttti[i];
     double dQdrh     = (ki[i] == 0) ? 0.0
                      : -alpi[i]*ki[i]*pow(del, (double) (ki[i]-1))/rhoi[i]; 
     double d2Qdt2    = -beti[i]*2.0/SQUARE(ttti[i]);
     double d2Qdtdrh  = 0.0;
     double d2Qdrh2   = (ki[i] == 0 || ki[i] == 1) ? 0.0
                      : -alpi[i]*ki[i]*(ki[i]-1.0)*pow(del, (double) (ki[i]-2))
                        /SQUARE(rhoi[i]);
     double d3Qdt3    = 0.0;
     double d3Qdt2drh = 0.0;
     double d3Qdtdrh2 = 0.0;
     double d3Qdrh3   = (ki[i] == 0 || ki[i] == 1 || ki[i] == 2) ? 0.0
                      : -alpi[i]*ki[i]*(ki[i]-1.0)*(ki[i]-2.0)
                        *pow(del, (double) (ki[i]-3))/CUBICE(rhoi[i]);
     double expQ      = (Q > -100.0) ? exp(Q) : 0.0;
     Ar         += gi[i]*pow(del, (double) li[i])*expQ;
     dArdt      += gi[i]*pow(del, (double) li[i])*expQ*dQdt;
     dArdrh     += (li[i] == 0) ? gi[i]*expQ*dQdrh
                 : gi[i]*li[i]*pow(del, (double) (li[i]-1))*expQ/rhoi[i]
                   + gi[i]*pow(del, (double) li[i])*expQ*dQdrh;
     d2Ardt2    += gi[i]*pow(del, (double) li[i])*expQ*(SQUARE(dQdt)+d2Qdt2);
     d2Ardtdrh  += (li[i] == 0) ? gi[i]*expQ*dQdt*dQdrh + gi[i]*expQ*d2Qdtdrh
                 : gi[i]*li[i]*pow(del, (double) (li[i]-1))*expQ*dQdt/rhoi[i]
                   + gi[i]*pow(del, (double) li[i])*expQ*(dQdt*dQdrh+d2Qdtdrh);
     if (li[i] == 0) d2Ardrh2 += gi[i]*expQ*SQUARE(dQdrh) + gi[i]*expQ*d2Qdrh2;
     else if (li[i] == 1) d2Ardrh2 += 2.0*gi[i]*expQ*dQdrh/rhoi[i]
         + gi[i]*del*expQ*SQUARE(dQdrh) + gi[i]*del*expQ*d2Qdrh2;
     else d2Ardrh2 += gi[i]*li[i]*(li[i]-1.0)*pow(del, (double) (li[i]-2))
           *expQ/SQUARE(rhoi[i])
         + 2.0*gi[i]*li[i]*pow(del, (double) (li[i]-1))*expQ*dQdrh/rhoi[i]
         + gi[i]*pow(del, (double) li[i])*(expQ*SQUARE(dQdrh)+expQ*d2Qdrh2);
     d3Ardt3    += gi[i]*pow(del, (double) li[i])
                     *expQ*(CUBICE(dQdt)+3.0*dQdt*d2Qdt2+d3Qdt3);
     d3Ardt2drh += (li[i] == 0) ? 
         gi[i]*(expQ*SQUARE(dQdt)*dQdrh + expQ*(d2Qdt2*dQdrh + dQdt*d2Qdtdrh)) 
       + gi[i]*(expQ*dQdt*d2Qdtdrh + expQ*d3Qdt2drh)
       : gi[i]*li[i]*pow(del, (double) (li[i]-1))
          *(expQ*SQUARE(dQdt)+expQ*d2Qdt2)/rhoi[i]
       + gi[i]*pow(del, (double) li[i])*(expQ*dQdt*(dQdt*dQdrh+d2Qdtdrh)
         + expQ*(d2Qdt2*dQdrh+dQdt*d2Qdtdrh+d3Qdt2drh));
     d3Ardtdrh2 += (li[i] == 1) ? 2.0*gi[i]*li[i]*(expQ*dQdt*dQdrh 
           + expQ*d2Qdtdrh)/rhoi[i]
         + gi[i]*del*(expQ*dQdt*SQUARE(dQdrh) + expQ*2.0*dQdrh*d2Qdtdrh 
           + expQ*dQdt*d2Qdrh2 + expQ*d3Qdtdrh2)
       : gi[i]*li[i]*(li[i]-1.0)*pow(del, (double) (li[i]-2))
           *expQ*dQdt/SQUARE(rhoi[i])
         + 2.0*gi[i]*li[i]*pow(del, (double) (li[i]-1))
           *(expQ*dQdt*dQdrh + expQ*d2Qdtdrh)/rhoi[i]
         + gi[i]*pow(del, (double) li[i])*(expQ*dQdt*SQUARE(dQdrh)
           + expQ*2.0*dQdrh*d2Qdtdrh + expQ*dQdt*d2Qdrh2 + expQ*d3Qdtdrh2);
     if (li[i] == 0) d3Ardrh3 += gi[i]*(expQ*CUBICE(dQdrh) 
           + expQ*2.0*dQdrh*d2Qdrh2 + expQ*dQdrh*d2Qdrh2 + expQ*d3Qdrh3);
     else if (li[i] == 1) d3Ardrh3 += 3.0*gi[i]*(expQ*SQUARE(dQdrh) 
           + expQ*d2Qdrh2)/rhoi[i]
         + gi[i]*del*(expQ*CUBICE(dQdrh) + expQ*2.0*dQdrh*d2Qdrh2 
           + expQ*dQdrh*d2Qdrh2 + expQ*d3Qdrh3);
     else if (li[i] == 2) d3Ardrh3 += 3.0*gi[i]*2.0*expQ*dQdrh/SQUARE(rhoi[i])
         + 3.0*gi[i]*2.0*del*(expQ*SQUARE(dQdrh) + expQ*d2Qdrh2)/rhoi[i]
         + gi[i]*SQUARE(del)*(expQ*CUBICE(dQdrh) + expQ*2.0*dQdrh*d2Qdrh2 
           + expQ*dQdrh*d2Qdrh2 + expQ*d3Qdrh3);
     else d3Ardrh3 += gi[i]*li[i]*(li[i]-1.0)*(li[i]-2.0)
           *pow(del, (double) (li[i]-3))*expQ/CUBICE(rhoi[i])
         + 3.0*gi[i]*li[i]*(li[i]-1.0)*pow(del, (double) (li[i]-2))
           *expQ*dQdrh/SQUARE(rhoi[i])
         + 3.0*gi[i]*li[i]*pow(del, (double) (li[i]-1))
           *(expQ*SQUARE(dQdrh) + expQ*d2Qdrh2)/rhoi[i]
         + gi[i]*pow(del, (double) li[i])*(expQ*CUBICE(dQdrh) 
           + expQ*2.0*dQdrh*d2Qdrh2 + expQ*dQdrh*d2Qdrh2 + expQ*d3Qdrh3);
   } 
 
   /* calculate ideal gas function                                           */ 
 
   tr        = t/1.0e2;
   dtrdt     = 1.0/100.0;

   Z         = 1.0 + (ci[1]/tr + ci[2])*log(tr);
   dZdt      = (-ci[1]*dtrdt/SQUARE(tr))*log(tr) + (ci[1]/tr + ci[2])*dtrdt/tr;
   d2Zdt2    = (2.0*ci[1]*SQUARE(dtrdt)/CUBICE(tr))*log(tr) 
             + (-ci[1]*dtrdt/SQUARE(tr))*dtrdt/tr 
             + (-ci[1]*dtrdt/SQUARE(tr))*dtrdt/tr
             - (ci[1]/tr + ci[2])*SQUARE(dtrdt/tr);
   d3Zdt3    = (-3.0*2.0*ci[1]*CUBICE(dtrdt)/QUARTIC(tr))*log(tr)
             + (2.0*ci[1]*SQUARE(dtrdt)/CUBICE(tr))*dtrdt/tr
              + (2.0*ci[1]*SQUARE(dtrdt)/CUBICE(tr))*dtrdt/tr 
             - (-ci[1]*dtrdt/SQUARE(tr))*SQUARE(dtrdt/tr) 
             + (2.0*ci[1]*SQUARE(dtrdt)/CUBICE(tr))*dtrdt/tr
             - (-ci[1]*dtrdt/SQUARE(tr))*SQUARE(dtrdt/tr)
             - (-ci[1]*dtrdt/SQUARE(tr))*SQUARE(dtrdt/tr)
             + 2.0*(ci[1]/tr + ci[2])*CUBICE(dtrdt/tr);

   for (i=3; i<=18; i++) { 
     Z      += ci[i]*pow(tr, (double) (i-6)); 
     dZdt   += (i-6.0)*ci[i]*pow(tr, (double) (i-7))*dtrdt;
     d2Zdt2 += (i-6.0)*(i-7.0)*ci[i]*pow(tr, (double) (i-8))*SQUARE(dtrdt);
     d3Zdt3 += (i-6.0)*(i-7.0)*(i-8.0)*ci[i]*pow(tr, (double) (i-9))
               *CUBICE(dtrdt);
   }

   Ai         = -r*t*Z;
   dAidt      = - r*Z - r*t*dZdt;
   d2Aidt2    = - 2.0*r*dZdt - r*t*d2Zdt2;
   d3Aidt3    = - 3.0*r*d2Zdt2 - r*t*d3Zdt3;
 
   /* Calculate the sum                                                      */

   A         = Ab         + Ar         + Ai;
   dAdt      = dAbdt      + dArdt      + dAidt;
   dAdrh     = dAbdrh     + dArdrh;
   d2Adt2    = d2Abdt2    + d2Ardt2    + d2Aidt2;
   d2Adtdrh  = d2Abdtdrh  + d2Ardtdrh;
   d2Adrh2   = d2Abdrh2   + d2Ardrh2;
   d3Adt3    = d3Abdt3    + d3Ardt3    + d3Aidt3;
   d3Adt2drh = d3Abdt2drh + d3Ardt2drh;
   d3Adtdrh2 = d3Abdtdrh2 + d3Ardtdrh2;
   d3Adrh3   = d3Abdrh3   + d3Ardrh3;

   /* calculate g = A + p/rh  and  v = 1/rh */

   p           = rh*rh*dAdrh;
   dpdrh       = 2.0*rh*dAdrh + rh*rh*d2Adrh2;
   dpdt        = rh*rh*d2Adtdrh;
   drhdt       = -dpdt/dpdrh;
   d2pdt2      = rh*rh*d3Adt2drh;
   d2pdtdrh    = 2.0*rh*d2Adtdrh + rh*rh*d3Adtdrh2;
   d2pdrh2     = 2.0*dAdrh + 4.0*rh*d2Adrh2 + rh*rh*d3Adrh3;
 
   *gH2O       = A + p/rh;
   *hH2O       = A + p/rh - t*dAdt; /* g + ts */
   *sH2O       = - dAdt;
   *cpH2O      = -t*d2Adt2 + (t/(rh*rh))*SQUARE(dpdt)/(dpdrh);
   *vH2O       = 1.0/rh;
   *dvdtH2O    = -(1.0/SQUARE(rh))*drhdt;
   *dvdpH2O    = -(1.0/SQUARE(rh))/dpdrh;
          temp = (-d2pdt2 - 2.0*d2pdtdrh*drhdt - d2pdrh2*SQUARE(drhdt))/dpdrh;
   *d2vdt2H2O  = 2.0*SQUARE(drhdt)/CUBICE(rh) - temp/SQUARE(rh);
          temp = (-d2pdtdrh/dpdrh - d2pdrh2*drhdt/dpdrh)/dpdrh;
   *d2vdtdpH2O = -2.0*dpdt/(CUBICE(rh)*SQUARE(dpdrh)) - temp/SQUARE(rh);
   *d2vdp2H2O  = (1.0/SQUARE(rh))*d2pdrh2/CUBICE(dpdrh) 
               + (2.0/CUBICE(rh))/SQUARE(dpdrh);
          temp = - d2Adt2 - t*d3Adt3 + SQUARE(dpdt/rh)/dpdrh
               + t*(2.0*dpdrh*dpdt*d2pdt2 - SQUARE(dpdt)*d2pdtdrh)
                 /SQUARE(rh*dpdrh);
   *dcpdtH2O   = temp + t*(*d2vdt2H2O)*dpdt;

   *gH2O       *= 1.80152;
   *hH2O       *= 1.80152;
   *sH2O       *= 1.80152;
   *cpH2O      *= 1.80152;
   *dcpdtH2O   *= 1.80152;
   *vH2O       *= 18.0152;
   *dvdtH2O    *= 18.0152;
   *dvdpH2O    *= 18.0152;
   *d2vdt2H2O  *= 18.0152;
   *d2vdtdpH2O *= 18.0152;
   *d2vdp2H2O  *= 18.0152;

   *gH2O       += -285829.96 - (298.15*69.9146) - gref;
   *hH2O       += -285829.96                    - href;

#ifdef DEBUG
   printf("t, p: %g, %g\n", t, p);
   printf("   g:       %g\n", *gH2O);
   printf("   h:       %g\n", *hH2O);
   printf("   s:       %g\n", *sH2O);
   printf("   cp:      %g\n", *cpH2O);
   printf("   dcpdt:   %g\n", *dcpdtH2O);
   printf("   v:       %g\n", *vH2O);
   printf("   dvdt:    %g\n", *dvdtH2O);
   printf("   dvdp:    %g\n", *dvdpH2O);
   printf("   d2vdt2:  %g\n", *d2vdt2H2O);
   printf("   d2vdtdp: %g\n", *d2vdtdpH2O);
   printf("   d2vdp2:  %g\n", *d2vdp2H2O);
#endif
}

/*

    double a = -33676.0 + phase->h/r, b = 18.3527 - phase->s/r;
    double phiP = (0.110/t + 4.432e-5 + 1.405e-7*t - 2.394e-11*t*t)*p
        + (7.337e-8/t - 1.170e-8 - 9.502e-13*t)*p*p
        + (1.876e-10/t + 4.586e-13)*CUBE(p) - 1.191e-14*QUARTIC(p)/t;
    double dphiPdt = (-0.110/SQUARE(t) + 1.405e-7 - 2.0*2.394e-11*t)*p
        + (-7.337e-8/SQUARE(t) - 9.502e-13)*p*p
        - 1.876e-10*CUBE(p)/SQUARE(t) + 1.191e-14*QUARTIC(p)/SQUARE(t);
    double d2phiPdt2 = (2.0*0.110/CUBE(t) - 2.0*2.394e-11)*p
        + 2.0*7.337e-8*p*p/CUBE(t) + 2.0*1.876e-10*CUBE(p)/CUBE(t)
        - 2.0*1.191e-14*QUARTIC(p)/CUBE(t);
    double d3phiPdt3 = -6.0*0.110*p/QUARTIC(t)
        - 6.0*7.337e-8*p*p/QUARTIC(t) - 6.0*1.876e-10*CUBE(p)/QUARTIC(t)
        + 6.0*1.191e-14*QUARTIC(p)/QUARTIC(t);
    double gRobie = r*t*(2.9147*log(t) - 9.6863e-4*t + 6.8593e-8*t*t
        + 77.8899/sqrt(t) - 28954.8/t - 2263.27/SQUARE(t) - 15.8997);
    double dgRobiedt = r*(2.9147*log(t) - 9.6863e-4*t + 6.8593e-8*t*t
        + 77.8899/sqrt(t) - 28954.8/t - 2263.27/SQUARE(t) - 15.8997) 
        + r*t*(2.9147/t - 9.6863e-4 + 2.0*6.8593e-8*t - 0.5*77.8899/pow(t,
            (double) 1.5) + 28954.8/SQUARE(t) + 2.0*2263.27/CUBE(t));
    double d2gRobiedt2 = 2.0*r*(2.9147/t - 9.6863e-4 + 2.0*6.8593e-8*t
        - 0.5*77.8899/pow(t, (double) 1.5) + 28954.8/SQUARE(t)
        + 2.0*2263.27/CUBE(t)) + r*t*(-2.9147/SQUARE(t) + 2.0*6.8593e-8
        + 1.5*0.5*77.8899/pow(t, (double) 2.5) - 2.0*28954.8/CUBE(t)
        - 6.0*2263.27/QUARTIC(t));
    double d3gRobiedt3 = 3.0*r*(-2.9147/SQUARE(t) + 2.0*6.8593e-8
        + 1.5*0.5*77.8899/pow(t, (double) 2.5) - 2.0*28954.8/CUBE(t)
        - 6.0*2263.27/QUARTIC(t)) + r*t*(2.0*2.9147/CUBE(t)
        - 2.5*1.5*0.5*77.8899/pow(t, (double) 3.5) + 6.0*28954.8/QUARTIC(t)
        + 24.0*2263.27/QUINTIC(t));
    double gH2O, hH2O, sH2O, cpH2O, dcpdtH2O, vH2O, dvdtH2O, dvdpH2O,
    d2vdt2H2O, d2vdtdpH2O, d2vdp2H2O, dgdt, d2gdt2, d3gdt3;
            
    whaar(1.0, t, &gH2O, &hH2O, &sH2O, &cpH2O, &dcpdtH2O, &vH2O, &dvdtH2O,
                  &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O);
            
    gl     = r*t*(a/t + b + phiP) + gH2O - gRobie;
    dgdt   = r*(a/t + b + phiP) + r*t*(-a/SQUARE(t) + dphiPdt) - dgRobiedt;
    sl     = sH2O - dgdt;
    hl     = gl + t*sl;
    d2gdt2 = 2.0*r*(-a/SQUARE(t) + dphiPdt) + r*t*(2.0*a/CUBE(t) + d2phiPdt2) - d2gRobiedt2;
    d3gdt3 = 3.0*r*(2.0*a/CUBE(t) + d2phiPdt2) + r*t*(-6.0*a/QUARTIC(t) + d3phiPdt3) - d3gRobiedt3;
            
    cpl    = cpH2O - t*d2gdt2;
    dcpldt = dcpdtH2O - d2gdt2 - t*d3gdt3;
            
    vl       = r*(0.110 + 4.432e-5*t + 1.405e-7*t*t - 2.394e-11*CUBE(t))
        + 2.0*r*(7.337e-8 - 1.170e-8*t - 9.502e-13*t*t)*p
        + 3.0*r*(1.876e-10 + 4.586e-13*t)*p*p
        - 4.0*r*1.191e-14*CUBE(p);
    dvldt    = r*(4.432e-5 + 2.0*1.405e-7*t - 3.0*2.394e-11*t*t)
        - 2.0*r*(1.170e-8 + 2.0*9.502e-13*t)*p + 3.0*r*4.586e-13*p*p;
    dvldp    = 2.0*r*(7.337e-8 - 1.170e-8*t - 9.502e-13*t*t)
        + 6.0*r*(1.876e-10 + 4.586e-13*t)*p - 12.0*r*1.191e-14*p*p;
    d2vldt2  = r*(2.0*1.405e-7 - 6.0*2.394e-11*t) - 4.0*r*9.502e-13*p;
    d2vldp2  = 6.0*r*(1.876e-10 + 4.586e-13*t) - 24.0*r*1.191e-14*p;
    d2vldtdp = - 2.0*r*(1.170e-8 + 2.0*9.502e-13*t) + 6.0*r*4.586e-13*p;

*/

static double rMELTS_ZR_g(double t, double p) {
    double result = 0.0;
    double a = -33676.0, b = 18.3527;
    double phiP = (0.110/t + 4.432e-5 + 1.405e-7*t - 2.394e-11*t*t)*p
        + (7.337e-8/t - 1.170e-8 - 9.502e-13*t)*p*p
        + (1.876e-10/t + 4.586e-13)*CUBICE(p) - 1.191e-14*QUARTIC(p)/t;
    double gRobie = r*t*(2.9147*log(t) - 9.6863e-4*t + 6.8593e-8*t*t
        + 77.8899/sqrt(t) - 28954.8/t - 2263.27/SQUARE(t) - 15.8997);
    double gH2O, hH2O, sH2O, cpH2O, dcpdtH2O, vH2O, dvdtH2O, dvdpH2O, 
        d2vdt2H2O, d2vdtdpH2O, d2vdp2H2O;
    whaar(1.0, t, &gH2O, &hH2O, &sH2O, &cpH2O, &dcpdtH2O, &vH2O, &dvdtH2O,
                  &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O);
    result += r*t*(a/t + b + phiP) + gH2O - gRobie;
    return result;
}

static double rMELTS_ZR_dgdt(double t, double p) {
    double result = 0.0;
    double a = -33676.0, b = 18.3527;
    double phiP = (0.110/t + 4.432e-5 + 1.405e-7*t - 2.394e-11*t*t)*p
        + (7.337e-8/t - 1.170e-8 - 9.502e-13*t)*p*p
        + (1.876e-10/t + 4.586e-13)*CUBICE(p) - 1.191e-14*QUARTIC(p)/t;
    double dphiPdt = (-0.110/SQUARE(t) + 1.405e-7 - 2.0*2.394e-11*t)*p
        + (-7.337e-8/SQUARE(t) - 9.502e-13)*p*p
        - 1.876e-10*CUBICE(p)/SQUARE(t) + 1.191e-14*QUARTIC(p)/SQUARE(t);
    double dgRobiedt = r*(2.9147*log(t) - 9.6863e-4*t + 6.8593e-8*t*t
        + 77.8899/sqrt(t) - 28954.8/t - 2263.27/SQUARE(t) - 15.8997) 
        + r*t*(2.9147/t - 9.6863e-4 + 2.0*6.8593e-8*t - 0.5*77.8899/pow(t,
            (double) 1.5) + 28954.8/SQUARE(t) + 2.0*2263.27/CUBICE(t));
    double gH2O, hH2O, sH2O, cpH2O, dcpdtH2O, vH2O, dvdtH2O, dvdpH2O, 
        d2vdt2H2O, d2vdtdpH2O, d2vdp2H2O;
    whaar(1.0, t, &gH2O, &hH2O, &sH2O, &cpH2O, &dcpdtH2O, &vH2O, &dvdtH2O,
                  &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O);
    result += r*(a/t + b + phiP) + r*t*(-a/SQUARE(t) + dphiPdt) - sH2O - dgRobiedt;
    return result;
}

static double rMELTS_ZR_dgdp(double t, double p) {
    double result = r*(0.110 + 4.432e-5*t + 1.405e-7*t*t - 2.394e-11*CUBICE(t))
        + 2.0*r*(7.337e-8 - 1.170e-8*t - 9.502e-13*t*t)*p
        + 3.0*r*(1.876e-10 + 4.586e-13*t)*p*p
        - 4.0*r*1.191e-14*CUBICE(p);
    return result;
}

static double rMELTS_ZR_d2gdt2(double t, double p) {
    double result = 0.0;
    double a = -33676.0;
    double dphiPdt = (-0.110/SQUARE(t) + 1.405e-7 - 2.0*2.394e-11*t)*p
        + (-7.337e-8/SQUARE(t) - 9.502e-13)*p*p
        - 1.876e-10*CUBICE(p)/SQUARE(t) + 1.191e-14*QUARTIC(p)/SQUARE(t);
    double d2phiPdt2 = (2.0*0.110/CUBICE(t) - 2.0*2.394e-11)*p
        + 2.0*7.337e-8*p*p/CUBICE(t) + 2.0*1.876e-10*CUBICE(p)/CUBICE(t)
        - 2.0*1.191e-14*QUARTIC(p)/CUBICE(t);
    double d2gRobiedt2 = 2.0*r*(2.9147/t - 9.6863e-4 + 2.0*6.8593e-8*t
        - 0.5*77.8899/pow(t, (double) 1.5) + 28954.8/SQUARE(t)
        + 2.0*2263.27/CUBICE(t)) + r*t*(-2.9147/SQUARE(t) + 2.0*6.8593e-8
        + 1.5*0.5*77.8899/pow(t, (double) 2.5) - 2.0*28954.8/CUBICE(t)
        - 6.0*2263.27/QUARTIC(t));
    double gH2O, hH2O, sH2O, cpH2O, dcpdtH2O, vH2O, dvdtH2O, dvdpH2O, 
        d2vdt2H2O, d2vdtdpH2O, d2vdp2H2O;
    whaar(1.0, t, &gH2O, &hH2O, &sH2O, &cpH2O, &dcpdtH2O, &vH2O, &dvdtH2O,
                  &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O);
    result += 2.0*r*(-a/SQUARE(t) + dphiPdt) + r*t*(2.0*a/CUBICE(t) + d2phiPdt2) 
            - cpH2O/t - d2gRobiedt2;
    return result;
}

static double rMELTS_ZR_d2gdtdp(double t, double p) {
    double result = r*(4.432e-5 + 2.0*1.405e-7*t - 3.0*2.394e-11*t*t)
        - 2.0*r*(1.170e-8 + 2.0*9.502e-13*t)*p + 3.0*r*4.586e-13*p*p;
    return result;
}

static double rMELTS_ZR_d2gdp2(double t, double p) {
    double result = 2.0*r*(7.337e-8 - 1.170e-8*t - 9.502e-13*t*t)
        + 6.0*r*(1.876e-10 + 4.586e-13*t)*p - 12.0*r*1.191e-14*p*p;
    return result;
}

static double rMELTS_ZR_d3gdt3(double t, double p) {
    double result = 0.0;
    double a = -33676.0;
    double d2phiPdt2 = (2.0*0.110/CUBICE(t) - 2.0*2.394e-11)*p
        + 2.0*7.337e-8*p*p/CUBICE(t) + 2.0*1.876e-10*CUBICE(p)/CUBICE(t)
        - 2.0*1.191e-14*QUARTIC(p)/CUBICE(t);
    double d3phiPdt3 = -6.0*0.110*p/QUARTIC(t)
        - 6.0*7.337e-8*p*p/QUARTIC(t) - 6.0*1.876e-10*CUBICE(p)/QUARTIC(t)
        + 6.0*1.191e-14*QUARTIC(p)/QUARTIC(t);
    double d3gRobiedt3 = 3.0*r*(-2.9147/SQUARE(t) + 2.0*6.8593e-8
        + 1.5*0.5*77.8899/pow(t, (double) 2.5) - 2.0*28954.8/CUBICE(t)
        - 6.0*2263.27/QUARTIC(t)) + r*t*(2.0*2.9147/CUBICE(t)
        - 2.5*1.5*0.5*77.8899/pow(t, (double) 3.5) + 6.0*28954.8/QUARTIC(t)
        + 24.0*2263.27/QUINTIC(t));
    double gH2O, hH2O, sH2O, cpH2O, dcpdtH2O, vH2O, dvdtH2O, dvdpH2O, 
        d2vdt2H2O, d2vdtdpH2O, d2vdp2H2O;
    whaar(1.0, t, &gH2O, &hH2O, &sH2O, &cpH2O, &dcpdtH2O, &vH2O, &dvdtH2O,
                  &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O);
    result += 3.0*r*(2.0*a/CUBICE(t) + d2phiPdt2) + r*t*(-6.0*a/QUARTIC(t) + d3phiPdt3) 
            + cpH2O/t/t - dcpdtH2O/t - d3gRobiedt3;
    return result;
}

static double rMELTS_ZR_d3gdt2dp(double t, double p) {
    double result = r*(2.0*1.405e-7 - 6.0*2.394e-11*t) - 4.0*r*9.502e-13*p;
    return result;
}

static double rMELTS_ZR_d3gdtdp2(double t, double p) {
    double result = - 2.0*r*(1.170e-8 + 2.0*9.502e-13*t) + 6.0*r*4.586e-13*p;
    return result;
}

static double rMELTS_ZR_d3gdp3(double t, double p) {
    double result = 6.0*r*(1.876e-10 + 4.586e-13*t) - 24.0*r*1.191e-14*p;
    return result;
}


static double rMELTS_ZR_s(double T, double P) {
    double result = -rMELTS_ZR_dgdt(T, P);
    return result;
}

static double rMELTS_ZR_v(double T, double P) {
    double result = rMELTS_ZR_dgdp(T, P);
    return result;
}

static double rMELTS_ZR_cv(double T, double P) {
    double result = -T*rMELTS_ZR_d2gdt2(T, P);
    double dvdt = rMELTS_ZR_d2gdtdp(T, P);
    double dvdp = rMELTS_ZR_d2gdp2(T, P);
    result += T*dvdt*dvdt/dvdp;
    return result;
}

static double rMELTS_ZR_cp(double T, double P) {
    double result = -T*rMELTS_ZR_d2gdt2(T, P);
    return result;
}

static double rMELTS_ZR_dcpdt(double T, double P) {
    double result = -T*rMELTS_ZR_d3gdt3(T, P) - rMELTS_ZR_d2gdt2(T, P);
    return result;
}

static double rMELTS_ZR_alpha(double T, double P) {
    double result = rMELTS_ZR_d2gdtdp(T, P)/rMELTS_ZR_dgdp(T, P);
    return result;
}

static double rMELTS_ZR_beta(double T, double P) {
    double result = -rMELTS_ZR_d2gdp2(T, P)/rMELTS_ZR_dgdp(T, P);
    return result;
}

static double rMELTS_ZR_K(double T, double P) {
    double result = -rMELTS_ZR_dgdp(T, P)/rMELTS_ZR_d2gdp2(T, P);
    return result;
}

static double rMELTS_ZR_Kp(double T, double P) {
    double result = rMELTS_ZR_dgdp(T, P);
    result *= rMELTS_ZR_d3gdp3(T, P);
    result /= pow(rMELTS_ZR_d2gdp2(T, P), 2.0);
    return result - 1.0;
}

static double rMELTS_ZR_dparam_g(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        result += 0.0;
        break;
    }
    return result;
}

static double rMELTS_ZR_dparam_dgdt(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        result += 0.0;
        break;
    }
    return result;
}

static double rMELTS_ZR_dparam_dgdp(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        result += 0.0;
        break;
    }
    return result;
}

static double rMELTS_ZR_dparam_d2gdt2(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        result += 0.0;
        break;
    }
    return result;
}

static double rMELTS_ZR_dparam_d2gdtdp(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        result += 0.0;
        break;
    }
    return result;
}

static double rMELTS_ZR_dparam_d2gdp2(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        result += 0.0;
        break;
    }
    return result;
}

static double rMELTS_ZR_dparam_d3gdt3(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        result += 0.0;
        break;
    }
    return result;
}

static double rMELTS_ZR_dparam_d3gdt2dp(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        result += 0.0;
        break;
    }
    return result;
}

static double rMELTS_ZR_dparam_d3gdtdp2(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        result += 0.0;
        break;
    }
    return result;
}

static double rMELTS_ZR_dparam_d3gdp3(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    case 0: /* T_r */ 
        result += 0.0;
        break;
    case 1: /* P_r */ 
        result += 0.0;
        break;
    }
    return result;
}

static int rMELTS_ZR_get_param_number(void) {
    return 2;
}

static const char *paramNames[2] = { "T_r", "P_r" };

static const char *paramUnits[2] = { "K", "bar" };

static const char **rMELTS_ZR_get_param_names(void) {
    return paramNames;
}

static const char **rMELTS_ZR_get_param_units(void) {
    return paramUnits;
}

static void rMELTS_ZR_get_param_values(double **values) {
    (*values)[0] = tr;
    (*values)[1] = pr;
}

static int rMELTS_ZR_set_param_values(double *values) {
    return 1;
}

static double rMELTS_ZR_get_param_value(int index) {
    double result = 0.0;
    switch (index) {
    case 0:
        result = tr;
        break;
    case 1:
        result = pr;
        break;
    default:
         break;
    }
    return result;
}

static int rMELTS_ZR_set_param_value(int index, double value) {
    int result = 1;
    switch (index) {
    case 0:
        break;
    case 1:
        break;
    default:
         break;
    }
    return result;
}

const char *H2O_rMELTS_ZR_calib_identifier(void) {
    return identifier;
}

const char *H2O_rMELTS_ZR_calib_name(void) {
    return "H2O";
}

const char *H2O_rMELTS_ZR_calib_formula(void) {
    return "H2O";
}

const double H2O_rMELTS_ZR_calib_mw(void) {
    return 18.0152;
}

static const double elmformula[106] = {
        0.0,2.0,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0
    };

const double *H2O_rMELTS_ZR_calib_elements(void) {
    return elmformula;
}

double H2O_rMELTS_ZR_calib_g(double T, double P) {
    return rMELTS_ZR_g(T, P);
}

double H2O_rMELTS_ZR_calib_dgdt(double T, double P) {
    return rMELTS_ZR_dgdt(T, P);
}

double H2O_rMELTS_ZR_calib_dgdp(double T, double P) {
    return rMELTS_ZR_dgdp(T, P);
}

double H2O_rMELTS_ZR_calib_d2gdt2(double T, double P) {
    return rMELTS_ZR_d2gdt2(T, P);
}

double H2O_rMELTS_ZR_calib_d2gdtdp(double T, double P) {
    return rMELTS_ZR_d2gdtdp(T, P);
}

double H2O_rMELTS_ZR_calib_d2gdp2(double T, double P) {
    return rMELTS_ZR_d2gdp2(T, P);
}

double H2O_rMELTS_ZR_calib_d3gdt3(double T, double P) {
    return rMELTS_ZR_d3gdt3(T, P);
}

double H2O_rMELTS_ZR_calib_d3gdt2dp(double T, double P) {
    return rMELTS_ZR_d3gdt2dp(T, P);
}

double H2O_rMELTS_ZR_calib_d3gdtdp2(double T, double P) {
    return rMELTS_ZR_d3gdtdp2(T, P);
}

double H2O_rMELTS_ZR_calib_d3gdp3(double T, double P) {
    return rMELTS_ZR_d3gdp3(T, P);
}

double H2O_rMELTS_ZR_calib_s(double T, double P) {
    return rMELTS_ZR_s(T, P);
}

double H2O_rMELTS_ZR_calib_v(double T, double P) {
    return rMELTS_ZR_v(T, P);
}

double H2O_rMELTS_ZR_calib_cv(double T, double P) {
    return rMELTS_ZR_cv(T, P);
}

double H2O_rMELTS_ZR_calib_cp(double T, double P) {
    return rMELTS_ZR_cp(T, P);
}

double H2O_rMELTS_ZR_calib_dcpdt(double T, double P) {
    return rMELTS_ZR_dcpdt(T, P);
}

double H2O_rMELTS_ZR_calib_alpha(double T, double P) {
    return rMELTS_ZR_alpha(T, P);
}

double H2O_rMELTS_ZR_calib_beta(double T, double P) {
    return rMELTS_ZR_beta(T, P);
}

double H2O_rMELTS_ZR_calib_K(double T, double P) {
    return rMELTS_ZR_K(T, P);
}

double H2O_rMELTS_ZR_calib_Kp(double T, double P) {
    return rMELTS_ZR_Kp(T, P);
}

int H2O_rMELTS_ZR_get_param_number(void) {
    return rMELTS_ZR_get_param_number();
}

const char **H2O_rMELTS_ZR_get_param_names(void) {
    return rMELTS_ZR_get_param_names();
}

const char **H2O_rMELTS_ZR_get_param_units(void) {
    return rMELTS_ZR_get_param_units();
}

void H2O_rMELTS_ZR_get_param_values(double **values) {
    rMELTS_ZR_get_param_values(values);
}

int H2O_rMELTS_ZR_set_param_values(double *values) {
    return rMELTS_ZR_set_param_values(values);
}

double H2O_rMELTS_ZR_get_param_value(int index) {
    return rMELTS_ZR_get_param_value(index);
}

int H2O_rMELTS_ZR_set_param_value(int index, double value) {
    return rMELTS_ZR_set_param_value(index, value);
}

double H2O_rMELTS_ZR_dparam_g(double T, double P, int index) {
    return rMELTS_ZR_dparam_g(T, P, index);
}

double H2O_rMELTS_ZR_dparam_dgdt(double T, double P, int index) {
    return rMELTS_ZR_dparam_dgdt(T, P, index);
}

double H2O_rMELTS_ZR_dparam_dgdp(double T, double P, int index) {
    return rMELTS_ZR_dparam_dgdp(T, P, index);
}

double H2O_rMELTS_ZR_dparam_d2gdt2(double T, double P, int index) {
    return rMELTS_ZR_dparam_d2gdt2(T, P, index);
}

double H2O_rMELTS_ZR_dparam_d2gdtdp(double T, double P, int index) {
    return rMELTS_ZR_dparam_d2gdtdp(T, P, index);
}

double H2O_rMELTS_ZR_dparam_d2gdp2(double T, double P, int index) {
    return rMELTS_ZR_dparam_d2gdp2(T, P, index);
}

double H2O_rMELTS_ZR_dparam_d3gdt3(double T, double P, int index) {
    return rMELTS_ZR_dparam_d3gdt3(T, P, index);
}

double H2O_rMELTS_ZR_dparam_d3gdt2dp(double T, double P, int index) {
    return rMELTS_ZR_dparam_d3gdt2dp(T, P, index);
}

double H2O_rMELTS_ZR_dparam_d3gdtdp2(double T, double P, int index) {
    return rMELTS_ZR_dparam_d3gdtdp2(T, P, index);
}

double H2O_rMELTS_ZR_dparam_d3gdp3(double T, double P, int index) {
    return rMELTS_ZR_dparam_d3gdp3(T, P, index);
}
