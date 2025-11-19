% BMB_04.m

% This MATLAB m.file contains the program for pancreatic beta-cells, 
% published by Bertram and Sherman, Bull. Math. Biol., 66:1313-1344, 2004.

% Variables:
%    v -- voltage
%    n -- activation variable for a delayed rectifier
%    c -- free cytosolic calcium concentration
%    cer -- concentration of free calcium in the endoplasmic reticulum
%    a -- fraction of activated KATP channels

function [t,x] = call_phantom();
   tspan= [0 120000];
   x0 = [-61 0 0.1 0.46 97];
%   [t,x] = ode113(@phantom,tspan,x0);
   [t,x] = ode15s(@phantom,tspan,x0);

%  change defaults for plotting
   set(0,'DefaultAxesFontSize',18)

   tsec=t/1000;

%  set up variable vectors
   v=x(:,1); n=x(:,2); c=x(:,3); a=x(:,4); cer=x(:,5);

%  Time course figure
   figure;

%  Top left panel
   subplot(2,2,1);
   h=plot(tsec,v,'black');
   set(h,'LineWidth',2);
   axis([0 120 -70 -10]);
   box off;
   ylabel('V (mV)');

%  Top right panel
   subplot(2,2,2);
   h=plot(tsec,c,'black');
   set(h,'LineWidth',2);
   axis([0 120 0 0.3]);
   box off;
   ylabel('c (\muM)');

%  Bottom left panel
   subplot(2,2,3);
   h=plot(tsec,cer,'black');
   set(h,'LineWidth',2);
   axis([0 120 70 130]);
   box off;
   xlabel('t (s)');
   ylabel('Ca_{er} (\muM)');

%  Bottom right panel
   subplot(2,2,4);
   h=plot(tsec,a,'black');
   set(h,'LineWidth',2);
   axis([0 120 0.4 0.5]);
   box off;
   xlabel('t (s)');
   ylabel('a');

%  Phase plane figure
   figure;

   subplot(2,1,1);
   h=plot(c,v,'black');
   set(h,'LineWidth',2);
   axis([0.1 0.2 -70 -10]);
   xlabel('Ca (\muM)');
   ylabel('V (mV)');
   
   subplot(2,1,2);
   h=plot(cer,v,'black');
   set(h,'LineWidth',2);
   axis([90 110 -70 -10]);
   xlabel('Ca_{er} (\muM)');
   ylabel('V (mV)');

%  Print data set into text file
%   set up a matrix of output. Must transpose the column vectors.
   out = [tsec.'; v.'; n.'; c.'; a.'; cer.'];
   fid = fopen('data.dat','w');
   fprintf(fid,'%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n',out);
   fclose(fid);

   function dxdt = phantom(t,x);
%     Parameters most likely to vary:
      gkatpbar=500;
      gkca=700;
      r=0.14;
      kserca=0.4;

%     Other parameters:

%     Ca parameters  (sigmav=cyt volume/ER volume)
      sigmav=5;
      kc=0.2;
      ip3=0;

%     ACh current
      gach=0;
      vach=0;

%     Ik
      vn=-16;
      vk=-75;
      taun=20;
      gk=3000;
      sn=5;

%     Ica
      vca=25;
      gca=1200;
      vm=-20;
      sm=12;

%     Ikca
      kd=0.3;

%     Miscellaneous
      iapp=0;
      lambda=1.25;
      cm=5300;

%     Calcium Handling
      f=0.01;
      fer=0.01;
      alpha=4.50e-6;

%     ER Calcium Handling
      perl=0.0005;
      dact=0.35;
      dip3=0.5;
      dinh=0.4;
      sa=0.1;
      taua=300000;

%     set up variable vectors
      v=x(1);
      n=x(2);
      c=x(3);
      a=x(4);
      cer=x(5);

%     Activation functions
      ninf = 1/(1+exp((vn-v)/sn));
      minf = 1/(1+exp((vm-v)/sm));
      ainf = 1/(1+exp((r-c)/sa));
      omega = 1/(1+(kd/c)^5);

%     Ionic currents
      ica = gca*minf*(v-vca);
      ikca = gkca*omega*(v-vk);
      ikatp = gkatpbar*a*(v-vk);
      ik = gk*n*(v-vk);
      iach = gach*(v-vach);

%     ER functions
      ainf_ip3 = 1/(1 + dact/c);
      hinfer = 1/(1 + c/dinh);
      jerp = kserca*c;
      binf = ip3/(ip3 + dip3);
      o = ainf_ip3^3*binf^3*hinfer^3;

%     Ca fluxes
      jmemtot = -(alpha*ica + kc*c);
      jerleak = perl*(cer - c);
      jerip3 = o*(cer - c);
      jertot = jerleak + jerip3 - jerp;

      vdot = -(ica + ik + ikatp + ikca + iach - iapp)/cm;
      ndot = lambda*(ninf - n)/taun;
      cdot = f*(jertot + jmemtot);
      adot = (ainf-a)/taua;
      cerdot = -fer*sigmav*jertot;
      dxdt=[vdot;ndot;cdot;adot;cerdot];
    end
end
     
