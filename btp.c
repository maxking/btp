/******************************************************************
Author: Abhilash Raj
Project: B.Tech Project
Session: 2013-14
Guide: Prof T. K. Kundu
Topic: Simulation of pellet induration cycle for iron ore pellets
***************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define PI 3.141

typedef struct Zones {
  double len;
  double inlet_gas_temp;
  double outlet_gas_temp;
  double pressure_drop;
  double area;
  bool flow_dir;
  int time;
} Zone;



int main() {

  Zone uddz, dddz, fz, cz, cz2, current_z;
  int iz, it;
  // Furnace Parameters
  double grate_length, grate_speed, z;

  // Process conditions
  double hgs, phi, vg, pg, Cpg, Cps, eb, clf, Tdry,  as,
    ds, fg, fs, Rw, Hw, Rc, Hc, Twc, Hwnb, Twb, time_spent,
    Xch2o, pc, Do2n2, Do2c, de, dc, Cco2e, Cco2, Dsco2,
    Dco2, EL, ug, Cu, dl, u0, rl, ps, el, nl, Tr2, Tr1, m,
    K1, K2, K3, KL, Kl, Kc, kmc, kpl, kpc, Sh, Sc, Sc2, Re,
    val, Rl, Hl, furnace_area, dp, Pdiff;

  // Simulation Variables
  double zstep, tstep;
  int n, time,  start, stop, incr;

  // Initialize the constants
  furnace_area = 400;
  Twc = 647.0;
  Hwnb = 2256000;
  Twb = 99.98 + 273.15;
  Xch2o = 0.05;
  grate_length = 105.0;
  // Furnace Properties Start
  uddz.len = 9.0;
  uddz.pressure_drop = 4000;
  uddz.area = uddz.len*furnace_area/grate_length;
  dddz.len = 22.5;
  dddz.pressure_drop = 4000;
  dddz.area = dddz.len*furnace_area/grate_length;
  fz.len = 78.0;
  fz.pressure_drop = 2900;
  fz.area = fz.len*furnace_area/grate_length;
  cz.len = 90.0;
  cz.pressure_drop = 4800;
  cz.area = cz.len*furnace_area/grate_length;
  grate_speed = 0.035;   // 2.1 m/min
  z = 0.55;
  // Furnace Properties End
  zstep = 0.01;
  tstep = 5.0;
  clf = 10;
  phi = 0.9;
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
  pg = 1.205;
  as = 6*(1-eb)/(ds*phi);
  n =  z/zstep;
  time = (int) grate_length/grate_speed;
  uddz.time = (int) uddz.len/grate_speed;
  dddz.time = (int) dddz.len/grate_speed;
  fz.time = (int) fz.len/grate_speed;
  cz.time = (int) grate_length/grate_speed;

  // Initialize the parameters that change during the course
  // of the process.
  double Ts[n][time], Xo2[n], Tg[n][time], Yh2o[n][time],
    Yc[n][time], P[n], Xh2o[n];

  // Set initial values for various variables
  Xo2[0] = 21.0;
  Xh2o[0] = 0.056;
  for(iz=0;iz<n;iz++)
	P[iz] = 0;
  // Initialize assumed inlet gas temperatures
  for(it=0; it<uddz.time; it++)
    Tg[0][it] = 189 + 273.15;
  for(it=uddz.time; it<dddz.time; it++)
    Tg[n-1][it] = 213 + 273.15;
  for(it=dddz.time; it<fz.time; it++)
    Tg[n-1][it] = 907 + 273.15;

  printf("Initialized Tg\n");
  // All solids enter the furnace at room temp i.e. 298K
  for(iz=0;iz<n;iz++) {
	for(it=0;it<time;it++) {
	  Ts[iz][it] = 298.15;
	  Yc[iz][it] = 1.4;
	  Yh2o[iz][it] = 0.131;
	}
  }
  printf("UDD_TIME is %d\n", uddz.time);
  printf("DDD_TIME is %d\n", dddz.time);
  printf("FZ_TIME is  %d\n", fz.time);
  printf("CZ_TIME is %d\n", cz.time);
  printf("Number of elements are %d\n", n);

  for(it=0; it<fz.time; it++) {

	if(it == fz.time)
	  vg = 0.3833;
    if(it<uddz.time){
	  current_z = uddz;
      dp = uddz.pressure_drop;
      start = 0;
      incr = 1;
      P[0] = 0;
      fg = vg*uddz.area*pg/eb;
    }
    else if(it<dddz.time && it>=uddz.time) {
	  current_z = dddz;
      dp = dddz.pressure_drop;
      start = n-1;
      incr = -1;
      P[n-1] = 0;
      fg = vg*dddz.area*pg/eb;
    }
    else if(it<fz.time && it>=dddz.time) {
	  current_z = fz;
      dp = fz.pressure_drop;
      start = n-1;
      incr = -1;
      P[n-1] = 0;
      fg = vg*fz.area*pg/eb;
    }
    else {
      dp = cz.pressure_drop;
      start = 0;
      incr = 1;
      P[0] = 0;
      fg = vg*cz.area*pg/eb;
    }
	/* while(abs(dp - abs(P[n-1] - P[0])) > 100) { */
	/* iz = start; */
	/* P[0] = 0; */
	/* fg = vg*current_z.area*pg/eb; */
	Pdiff = 0;
	while(abs(Pdiff - dp) > 100) {
	  iz = start;
	  while(iz>=0 && iz<n) {
      // Compute all variables;
	  Cpg = 881 + 0.31*Tg[iz][it] - 7.98e-5*pow(Tg[iz][it], 2);

	  //printf("Cpg = %lf Tg = %lf\n", Cpg, Tg[iz][it]);
	  ug = u0*pow((Tg[iz][it]/273.15), 1.5)*((Cu+273.15)/(Cu + Tg[iz][it]));

	  //printf("ug = %lf Tg = %lf \n", ug, Tg[iz][it]);

	  val = Yh2o[iz][it]/(1-Yh2o[iz][it]);
	  pg = 219.38*(1+val)/((0.622 + val)*Tg[iz][it]);

	  //printf("Tg = %lf", Tg[iz][it]);

	  Do2c = 0.435e-5*pow((Tg[iz][it]/298.15), 1.5)/(pow((0.29), -9.41));
	  Dco2 = 7.166666e-10*pow(Tg[iz][it], 1.75);
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


	  hgs = (phi*vg*pg*Cpg)/(6*(1-eb)*clf);

	  printf("hgs = %lf phi = %lf vg = %lf pg = %lf Cpg = %lf eb= %lf clf  = %f\n",
	         hgs, phi, vg, pg, Cpg, eb, clf);

	  Hl = 95.15909*(4.5*Ts[iz][it] - EL);

	  if (Ts[iz][it] < 1053) {
		Rl = 0;
	  }
	  else {
		K1 = ds*(ds-dl)/(dl*Dsco2);
		K2 = ds/(Sh* Dco2);
		K3 = pow((ds/dl),2)/kpl;

		KL = 1/(K1 + K2 + K3);

		Kl = exp(7.35 - 5211/Ts[iz][it]);
		Cco2e = 1000*Kl/(82.057*Ts[iz][it]);
		Rl = 44*4*PI*pow(rl,2)*nl*KL*(Cco2e - Cco2);
	  }

	  if (Tg[iz][it] > Twc) {
		Hw  = 3.1563e6 - 2396.6*Tg[iz][it];
	  }
	  else {
		Tr1 = Twb/Twc;
		Tr2 = Tg[iz][it]/Twc;
		m = (1-Tr2)/(1-Tr1);
		Hw = Hwnb*pow(m, 0.375);
	  }

	  Tdry = (Tg[iz][it]+Ts[iz][it])/2;
	  if (Yh2o[iz][it] < Xch2o)
		Rw = as*hgs*(Tg[iz][it]-Tdry)/Hw;
	  else
		Rw = as*hgs*(Tg[iz][it]-Tdry)*Xh2o[iz]/(Hw*Xch2o);

	  Xo2[iz+incr] = Xo2[iz] - Rc*zstep/fg;

	  Xh2o[iz+incr] = Xh2o[iz] + Rw*zstep/fg;

	  Tg[iz+incr][it] = Tg[iz][it] + (Rw*Hw - hgs*as*(Tg[iz][it] -
													Ts[iz][it]))*zstep/(fg*Cpg);

	  //printf("fg = %lf, pg = %lf, vg = %lf\n", fg, pg, vg);
	  //printf("%lf = %lf + (%lf*%lf - %lf*%lf(%lf-%lf))*%lf/(%lf*%lf)\n", Tg[iz+1][it]
			 //       , Tg[iz][it], Rw, Hw, hgs, as, Tg[iz][it], Ts[iz][it], zstep, fg, Cpg);

	  /* P[iz+incr] = P[iz] - ((150.0*pow(1-eb, 2)*ug*vg/(pow(de,2)*pow(eb, 3)*pow(phi, 2)) + */
	  /* 						 1.75*pg*(1-eb)*pow(vg, 2)/(de*pow(eb, 3)*phi))*zstep); */

	  //printf("P[%d] = %lf\n", iz+incr, P[iz]);
	  //printf(" DP = %lf \n", ((150.0*pow(1-eb, 2)*ug*vg/(pow(de,2)*pow(eb, 3)*pow(phi, 2)) +
	  //1.75*pg*(1-eb)*pow(vg, 2)/(de*pow(eb, 3)*phi)))*zstep);
	  //printf("vg = %lf de = %lf pg = %lf ug = %lf\n", vg, de, pg, ug);

	  Ts[iz][it+1] = Ts[iz][it] + (hgs*as*(Tg[iz][it] - Ts[iz][it])
								   - Rw*Hw + Rc*Hc - Rl*Hl )*tstep/(fs*Cps);

	  /* printf("Ts[%d][%d] = %lf, Tg = %lf, hgs = %lf, as = %lf, Rw = %lf, Hw = %lf, fs = %lf, Cps = %lf\n", */
	  /* iz, it, Ts[iz][it], Tg[iz][it], hgs, as, Rw, Hw, fs, Cps) */;
	  //printf("Tg = %lf \t Ts = %lf \n", Tg[iz][it], Ts[iz][it]);

	  //fg = fg + (Rw + Rc)*zstep;
	  //printf("fg = %lf \n", fg);

	  Yc[iz][it+1] = Yc[iz][it] - Rc*tstep/fs;
	  Yh2o[iz][it+1] = Yh2o[iz][it] - Rw*tstep/fs;
	  iz += incr;
	}
	/* Pdiff = (150.0*pow(1-eb, 2)*ug*vg/(pow(de,2)*pow(eb, 3)*pow(phi, 2)) + */
	Pdiff = 1.75*pg*(1-eb)*pow(vg, 2)/(de*pow(eb, 3)*phi)*z;
	vg = sqrt(pow(vg, 2)/Pdiff*dp);
	/* printf("%lf\t%lf\n", Pdiff, vg); */
  }

  }
  printf("\n");
  for(iz=n-1;iz>0;iz--){
    for(it=0;it<fz.time;it+=50)
      printf("%.0lf ",Ts[iz][it] - 273.15);
    printf("\n");
  }
  return 0;
}
