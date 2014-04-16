/******************************************************************
Author: Abhilash Raj
Project: B.Tech Project
Session: 2013-14
Guide: Prof T. K. Kundu
Topic: Simulation of pellet induration cycle for iron ore pellets
***************************************************************/
#include <stdio.h>
#include <math.h>

#define PI 3.141

int main() {
  // Furnace Parameters
  double grate_length, udd_length, ddd_length, rc_length, cz1_length,
    cz2_length, grate_speed, z;
  // Process conditions
  double hgs, phi, vg, pg, Cpg, Cps, eb, clf, Tdry,  as,
    ds, fg, fs, Rw, Hw, Rc, Hc, Twc, Hwnb, Twb,
    Xch2o, pc, Do2n2, Do2c, de, dc, Cco2e, Cco2, Dsco2,
    Dco2, EL, ug, Cu, dl, u0, rl, ps, el, nl, Tr2, Tr1, m,
    K1, K2, K3, KL, Kl, Kc, kmc, kpl, kpc, Sh, Sc, Sc2, Re,
    val, Rl, Hl;
  // Simulation Variables
  double zstep, tstep;
  int n, time, udd_time, ddd_time, rc_time, cz1_time, cz2_time;
  // Required Arrays


  Twc = 647.0;
  Hwnb = 2260000;
  Twb = 99.98 + 273.15;
  Xch2o = 0.05;
  grate_length = 105.0;
  udd_length = 9.0;
  ddd_length = 22.5;
  rc_length = 78.0;
  cz1_length = 90.0;
  grate_speed = 0.035;   // 2.1 m/min
  z = 0.55;
  zstep = 0.01;
  tstep = 1.0;
  clf = 0.8;
  phi = 1.0;
  vg = 0.3833;
  eb = 0.41;
  dc = 7.5e-5;
  ds = 0.012;
  dl = 7.5e-5;
  de = 0.012;
  u0 = 1.72e-5;
  Cu = 113.0;
  Do2n2 = 0.25;
  rl = dl/2;
  Hc = 28e6;
  Cco2 = 200.0;
  pc = 625.0;
  ps = 2.7e3;
  EL = 10.0;
  el = 0.3;
  nl = 1000.0;
  fs = ps;
  fg = 1300.0;
  as = 6*(1-eb)/(ds*phi);
  n =  z/zstep;
  time = (int) grate_length/grate_speed;
  udd_time = (int) udd_length/grate_speed;

  double Ts[n][time], Xo2[n], Tg[n], Yh2o[n][udd_time], Yc[n][udd_time], P[n], Xh2o[n];
  // Assumed inlet gas temperature for UDD zone
  Tg[0] = 189 + 273.15;
  Xo2[0] = 21.0;
  Xh2o[0] = 0.056;
  P[0] = 0;
  Yh2o[0][0] = 0.131;
  // All solids enter the furnace at room temp i.e. 298K
  int iz, it;
  for(iz=0;iz<n;iz++) {
      for(it=0;it<udd_time;it++) {
        Ts[iz][it] = 298.15;
        Yc[iz][it] = 1.4;
      }
  }

  printf("UDD_TIME is %d \n", udd_time);
  printf("Number of elements are %d \n", n);
  // Loop though all the elements of UDD zone
  for(it=0;it<(udd_time-1);it++) {
    //Loop through all the elements in the height for a length of 1s
    for(iz=0;iz<(n-1);iz++) {
      // Compute all variables;
      Cpg = 881 + 0.31*Tg[iz] - 7.98e-5*pow(Tg[iz], 2);

      ug = u0*pow((Tg[iz]/273.15), 1.5)*((Cu+273.15)/(Cu + Tg[iz]));

      val = Yh2o[iz][it]/(1-Yh2o[iz][it]);
      pg = 219.38*(1+val)/(0.622 + val);

      Do2c = 0.435e-5*pow((Tg[iz]/298.15), 1.5)/(pow((0.29), -9.41));
      Dco2 = 7.166666e-10*pow(Tg[iz], 1.75);
      Dsco2 = Dco2*(el/pow(el, -0.41));

      Re = de*vg*pg/ug;
      Sc = ug/(pg*Do2n2);
      Sc2 = ug/(pg*Dsco2);
      Sh = 2 + 0.6*pow(Re/eb, 0.5)*pow(Sc2, 0.333);

      if (Ts[iz][it] > 1050.15) {
        Cps = 999 + 0.0416*Ts[iz][it];
      }
      else if( Ts[iz][it] > 950.15 && Ts[iz][it] < 1050.15) {
        Cps = 1111;
      }
      else {
        Cps = 341.6 + 1.324*Ts[iz][it] - 4.032e-4*pow(Ts[iz][it], 2);
      }

      kpc = 595.6*Ts[iz][it]*exp(-17940.0/Ts[iz][it]);
      kpl = 9.12e3*exp(-40000/(1.987*Ts[iz][it]));
      kmc = Do2c*(2 + 0.6*pow(Sc,(1/3))*pow((Re/eb),(1/2)))/de;
      Kc = kpc*kmc/(kpc+kmc);
      Rc = 6*Yc[iz][it]*ps*Kc*Xo2[iz]/(100*pc*dc);

      K1 = ds*(ds-dl)/(dl*Dsco2);
      K2 = ds/(Sh* Dco2);
      K3 = pow((ds/dl),2)/kpl;

      KL = 1/(K1 + K2 + K3);

      Kl = exp(7.35 - 5211/Ts[iz][it]);
      Cco2e = 1000*Kl/(82.057*Ts[iz][it]);
      hgs = (phi*vg*pg*Cpg)/(6*(1-eb)*clf);


      Hl = 95.15909*(4.5*Ts[iz][it] - EL);
      if (Ts[iz][it] < 1053) {
        Rl = 0;
      }
      else {
        Rl = 44*4*PI*pow(rl,2)*nl*KL*(Cco2e - Cco2);
      }

      if (Tg[iz] > Twc) {
        Hw  = 3.1563e6 - 2396.6*Tg[iz];
      }
      else {
        Tr1 = Twb/Twc;
        Tr2 = Tg[iz]/Twc;
        m = (1-Tr2)/(1-Tr1);
        Hw = Hwnb*pow(m, 0.375);
      }

      Tdry = (Tg[iz]+Ts[iz][it])/2;
      if (Yh2o[iz][it] > Xch2o)
        Rw = as*hgs*(Tg[iz]-Tdry)/Hw;
      else
        Rw = as*hgs*(Tg[iz]-Tdry)*Xh2o[iz]/(Hw*Xch2o);

      Xo2[iz+1] = Xo2[iz] - Rc*zstep/fg;

      Xh2o[iz+1] = Xh2o[iz] + Rw*zstep/fg;

      Tg[iz+1] = Tg[iz] + (Rw*Hw - hgs*as*(Tg[iz] -Ts[iz][it]))*zstep/(fg*Cpg);

      P[iz+1] = P[iz] - (150*pow(1-eb, 2)*ug*vg/(pow(de,2)*pow(eb, 3)*pow(phi, 2)) +
                         1.75*pg*(1-eb)*pow(vg, 2)/(de*pow(eb, 3)*phi))*zstep;
      Ts[iz+1][it+1] = Ts[iz][it] + (hgs*as*(Tg[iz] - Ts[iz][it]) -Rw*Hw + Rc*Hc -Rl*Hl)*tstep/(fs*Cps);
      fg = fg + (Rw + Rc)*zstep;
      Yh2o[iz+1][it+1] = Yh2o[it][iz] - Rw*tstep/fs;
      Yc[iz+1][it+1] = Yc[it][iz] - Rc*tstep/fs;

    }
  }

  for(iz=0;iz<n;iz++){
    for(it=udd_time-1;it>=0;it--)
      printf("%lf ", iz, it, Ts[iz][it]);
    printf("\n");
  }

  return 0;
}
