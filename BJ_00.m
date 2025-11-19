% BJ_00.m

% The original phantom burster model, from "The Phantom Burster Model for Pancreatic
% Beta-Cells", by R. Bertram, J. Previte, A. Sherman, T. Kinard, and L. Satin, 
% Biophysical Journal, 79:2880-2892, 2001.

% For Fig. 2, set gs1=20.
% For Fig. 3, set gs1=7, set v(0)=0.6, set total integration time = 120000.
% For Fig. 4, set gs1=3, set v(0)=0.6, set total integration time = 300000.

% Units: V = mV; t = ms; g = pS; I = fA


function [t,x] = call_phantom()
   tspan= [0 90000];
   x0 = [-50 0 0.0 0.6];
   [t,x] = ode113(@phantom,tspan,x0);

%  change defaults for plotting 
   set(0,'DefaultAxesFontSize',18)

   tsec=t/1000;

%  set up variable vectors
   v=x(:,1); n=x(:,2); s1=x(:,3); s2=x(:,4);

%  Top panel
   subplot(3,1,1)
   h=plot(tsec,v,'black');
   set(h,'LineWidth',2);
   axis([0 90 -70 -10]);
   box off;
   title('Phantom model');
   ylabel('V (mV)');

%  Middle panel
   subplot(3,1,2);
   h=plot(tsec,s1,'black',tsec,s2,'red');
   set(h,'LineWidth',2);
   axis([0 90 0 1]);
   box off;
   xlabel('t (s)');
   ylabel('s_1 and s_2');

%  Bottom panel
   subplot(3,1,3);
   h=plot(s1,v,'black');
   set(h,'LineWidth',2);
   axis([0 1 -70 -10]);
   box off;
   xlabel('s_1');
   ylabel('V (mV)');

%  Print data set into text file
%   set up a matrix of output. Must transpose the column vectors.
   out = [tsec.'; v.'; n.'; s1.'; s2.'];
   fid = fopen('data.dat','w');
   fprintf(fid,'%9.5f %9.5f %9.5f %9.5f %9.5f\n',out);
   fclose(fid);

   function dxdt = phantom(t,x);
%     Parameters most likely to vary
      gs1=5;
      gs2=32;
      autos1=1;
      autos2=1;

%     Other parameters
      taus1=1000; taus2=120000; tnbar=9.09;
      vs1=-40; vs2=-42; ss1=0.5; ss2=0.4;
      s1knot=1; s2knot=1;
      gl=25; vl=-40;
      gk=1300; vk=-80;
      gca=280; vca=100;
      lambda=1.1;
      cm=4524;
      vm=-22; sm=7.5;
      vn=-9; sn=10;
      ss1=0.5; ss2=0.4;
      
%     set up variable vectors
      v=x(1); n=x(2); s1=x(3); s2=x(4);
      
%     activation and time-constant functions
      minf = 1.0/(1.0+exp((vm-v)/sm));
      ninf = 1.0/(1.0+exp((vn-v)/sn));
      taun = tnbar/(1.0+exp((v-vn)/sn));
      s1inf = 1.0/(1.0+exp((vs1-v)/ss1));
      s2inf = 1.0/(1.0+exp((vs2-v)/ss2));

%     ionic currents
      ica = gca*minf*(v-vca);
      ik = gk*n*(v-vk);
      il = gl*(v-vl);
      is1 = gs1*s1*(v-vk);
      is2 = gs2*s2*(v-vk);

      vdot = -(ica + ik + il + is1 + is2)/cm;
      ndot = lambda*(ninf - n)/taun;
      s1dot = (s1inf - s1)/taus1;
      s2dot = (s2inf - s2)/taus2;
      dxdt=[vdot;ndot;s1dot;s2dot];
    end
end
