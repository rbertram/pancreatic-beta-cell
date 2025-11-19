% BJ_04a.m

# This file contains the program for the Dual Oscillator Model (DOM)
# from "Calcium and Glycolysis Mediate Multiple Bursting Modes in Pancreatic
# Beta-Cells", by R. Bertram, L. Satin, M. Zhang, P. Smolen, and A. Sherman in
# Biophysical Journal, vol. 87, pp. 3074-3087, 2004. 

% Variables:
%    v -- voltage
%    n -- activation of delayed rectifier
%    c -- free cytosolic calcium concentration
%    adp -- cytosolic ADP concentration
%    cer -- concentration of free calcium in the endoplasmic reticulum
%    g6p -- glucose 6-phosphate concentration
%    fbp -- fructose 1,6-bisphosphate concentration

% Compound bursting: rgk=0.2, gkatp=25000, gkca=600, kg=10
% Glycolytic bursting: rgk=0.2, gkatp=27000, gkca=100, kg=10
% Simple bursting: rgk=0.4, gkatp=25000, gkca=600, kg=10
% Subthreshold oscillations: rgk=0.2, gkatp=30000, gkca=100, kg=10
% Accordion bursting: rgk=0.2, gkatp=23000, gkca=600, kg=10


function [t,x] = call_DOM();
   tmax = 600000;
   tspan= [0 tmax];
   x0 = [-60 0 0.1 808 130 269 0];
%   [t,x] = ode113(@DOM,tspan,x0);
   [t,x] = ode15s(@DOM,tspan,x0);

%  change defaults for plotting 
   set(0,'DefaultAxesFontSize',18)

   tsec=t/1000;
   tmin=tsec/60;
   tmax_min=tmax/60000;

%  set up variable vectors
   v=x(:,1); n=x(:,2); c=x(:,3); adp=x(:,4); cer=x(:,5); g6p=x(:,6); fbp=x(:,7);

%  Time course figure
   figure;

%  Top left panel
   subplot(2,2,1);
   h=plot(tmin,v,'black');
   set(h,'LineWidth',2);
   axis([0 tmax_min -70 -10]);
   box off;
   ylabel('V (mV)');

%  Top right panel
   subplot(2,2,2);
   h=plot(tmin,c,'black');
   set(h,'LineWidth',2);
   axis([0 tmax_min 0 0.2]);
   box off;
   ylabel('Ca (\muM)');

%  Bottom left panel
   subplot(2,2,3);
   h=plot(tmin,adp,'black');
   set(h,'LineWidth',2);
   axis([0 tmax_min 750 850]);
   box off;
   xlabel('t (m)');
   ylabel('ADP (\muM)');

%  Bottom right panel
   subplot(2,2,4);
   h=plot(tmin,fbp,'black');
   set(h,'LineWidth',2);
   axis([0 tmax_min 0 40]);
   box off;
   xlabel('t (m)');
   ylabel('FBP (\muM)');

%  Phase plane figure
   figure;

   subplot(2,1,1);
   h=plot(c,v,'black');
   set(h,'LineWidth',2);
   axis([0.0 0.15 -70 -10]);
   xlabel('Ca (\muM)');
   ylabel('V (mV)');
   
   subplot(2,1,2);
   h=plot(fbp,v,'black');
   set(h,'LineWidth',2);
   axis([0 40 -70 -10]);
   xlabel('FBP (\muM)');
   ylabel('V (mV)');

%  Print data set into text file
%   set up a matrix of output. Must transpose the column vectors.
   out = [tmin.'; v.'; n.'; c.'; adp.'; cer.'; g6p.'; fbp.'];
   fid = fopen('data.dat','w');
   fprintf(fid,'%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n',out);
   fclose(fid);

   function dxdt = DOM(t,x);

%       set up variable vectors
        v=x(1); n=x(2); c=x(3); adp=x(4); cer=x(5); g6p=x(6); fbp=x(7);

%       Parameters most likely to vary:
        Rgk=0.2; 
        gkatpbar=25000;
        gkca=600;
        kg=10;

%       Other parameters:

%       Membrane and Ca parameters  (sigmav=cyt volume/ER volume)
        sigmav=31; pleak=0.0002; kserca=0.4; lambdaer=1;
        epser=1; taun=20;

%       conversion parameter for glycolytic subsystem
        lambda=0.005;

        vk=-75; vca=25;
        gk=2700; gca=1000;
        cm=5300; kd=0.5;

%       calcium handling
        alpha=4.50e-6; kpmca=0.2;
        fcyt=0.01; fer=0.01;

%       KATP channel
        kdd=17; ktd=26; ktt=1;
        r=1; vg=2.2; taua=300000; r1=0.35;

%      glycolytic parameters
%        Rgk--glucokinase rate
%        atot--total adenine nucleotide concentration (micromolar)
%        k1--Kd for AMP binding
%        k2--Kd for FBP binding
%        k3--Kd for F6P binding
%        k4--Kd for ATP binding
%        famp,etc--Kd amplification factors for heterotropic binding
%        Rgpdh--glyceraldehyde phosphate dehydrogenase rate

        atot=3000;
        kg=10; k1=30; k2=1; k3=50000; k4=1000;
        famp=0.02; fatp=20; ffbp=0.2; fbt=20; fmt=20;
        pfkbas=0.06; cat=2; katpase=0.0003;

%       Electrical and Calcium Functions:

%       currents 
        minf = 1/(1+exp(-(20+v)/12));
        ninf = 1/(1+exp(-(16+v)/5));
        ikca = gkca/(1+(kd/c)^2)*(v-vk);
        ica = gca*minf*(v-vca);
        ik = gk*n*(v-vk);

%       calcium fluxes
        Jmem = -(alpha*ica + kpmca*c);
        Jserca = kserca*c;
        Jleak = pleak*(cer - c);
        Jer = epser*(Jleak - Jserca)/lambdaer;

%       Glycolytic and  Keizer-Magnus components:
     
%     nucleotide concentrations
        rad = sqrt((adp-atot)^2-4*adp^2);
        atp = 0.5*(atot-adp+rad);
        amp = adp^2/atp;
        f6p = 0.3*g6p;
        Rgpdh = 0.2*sqrt(fbp);

%        Iterative calculation of PFK
%        alpha=1 -- AMP bound
%        beta=1 -- FBP bound
%        gamma=1 -- F6P bound
%        delta=1 -- ATP bound

%      (alpha,beta,gamma,delta)
%      (0,0,0,0)
       weight1=1;
       topa1=0;
       bottom1=1;
          
%      (0,0,0,1)
       weight2=atp^2/k4;
       topa2=topa1;
       bottom2=bottom1+weight2;

%      (0,0,1,0)
       weight3=f6p^2/k3;
       topa3=topa2+weight3;
       bottom3=bottom2+weight3;

%      (0,0,1,1)
       weight4=(f6p*atp)^2/(fatp*k3*k4);
       topa4=topa3+weight4;
       bottom4=bottom3+weight4;

%      (0,1,0,0)
       weight5=fbp/k2;
       topa5=topa4;
       bottom5=bottom4+weight5;

%      (0,1,0,1)
       weight6=(fbp*atp^2)/(k2*k4*fbt);
       topa6=topa5;
       bottom6=bottom5+weight6;

%      (0,1,1,0)
       weight7=(fbp*f6p^2)/(k2*k3*ffbp);
       topa7=topa6+weight7;
       bottom7=bottom6+weight7;

%      (0,1,1,1)
       weight8=(fbp*f6p^2*atp^2)/(k2*k3*k4*ffbp*fbt*fatp);
       topa8=topa7+weight8;
       bottom8=bottom7+weight8;

%      (1,0,0,0)
       weight9=amp/k1;
       topa9=topa8;
       bottom9=bottom8+weight9;

%      (1,0,0,1)
       weight10=(amp*atp^2)/(k1*k4*fmt);
       topa10=topa9;
       bottom10=bottom9+weight10;

%     (1,0,1,0)
      weight11=(amp*f6p^2)/(k1*k3*famp);
      topa11=topa10+weight11;
      bottom11=bottom10+weight11;

%     (1,0,1,1)
      weight12=(amp*f6p^2*atp^2)/(k1*k3*k4*famp*fmt*fatp);
      topa12=topa11+weight12;
      bottom12=bottom11+weight12;

%     (1,1,0,0)
      weight13=(amp*fbp)/(k1*k2);
      topa13=topa12;
      bottom13=bottom12+weight13;

%     (1,1,0,1)
      weight14=(amp*fbp*atp^2)/(k1*k2*k4*fbt*fmt);
      topa14=topa13;
      bottom14=bottom13+weight14;

%     (1,1,1,0) -- the most active state of the enzyme
      weight15=(amp*fbp*f6p^2)/(k1*k2*k3*ffbp*famp);
      topa15=topa14;
      topb=weight15;
      bottom15=bottom14+weight15;

%     (1,1,1,1)
      weight16=(amp*fbp*f6p^2*atp^2)/(k1*k2*k3*k4*ffbp*famp*fbt*fmt*fatp);
      topa16=topa15+weight16;
      bottom16=bottom15+weight16;

%     Phosphofructokinase rate
      pfk=(pfkbas*cat*topa16 + cat*topb)/bottom16;

%     nucleotide concentrations
      rad = sqrt((adp-atot)^2-4*adp^2);
      atp = 0.5*(atot-adp+rad);
      amp = adp^2/atp;
      ratio = atp/adp;

%     KATP channel open probability
      mgadp = 0.165*adp;
      adp3m = 0.135*adp;
      atp4m = 0.05*atp;
      topo = 0.08*(1+2*mgadp/kdd) + 0.89*(mgadp/kdd)^2;
      bottomo = (1+mgadp/kdd)^2 * (1+adp3m/ktd+atp4m/ktt);
      katpo = topo/bottomo;
      ikatp = gkatpbar*katpo*(v-vk);

%     glycolytic input to adp
      y = vg*(Rgpdh/(kg+Rgpdh));
      fback = r+y;

      vdot = -(ik + ica + ikca + ikatp)/cm;
      ndot = (ninf - n)/taun;
      cdot = fcyt*(Jmem + Jer);
      adpdot = (atp-adp*exp(fback*(1-c/r1)))/taua;
      cerdot = -fer*sigmav*Jer;
      g6pdot = lambda*(Rgk - pfk);
      fbpdot = lambda*(pfk - 0.5*Rgpdh);
      dxdt=[vdot;ndot;cdot;adpdot;cerdot;g6pdot;fbpdot];
    end
end
