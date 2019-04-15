function pfit=JiDiapause(T,F,nrand)

% T should be 365 by at least 2
% F should be 365 by 1 (chl units)
% no nans

tstd=0.1*nanmean(T(:,1));
fstd=0.3*nanmean(F(:));

T=[T;T;T];
F=[F;F;F];

TX=repmat(T(:,1),1,nrand)+tstd.*rand(length(T),nrand);
FX=repmat(F,1,nrand)+fstd.*rand(length(F),nrand);
FX=max(FX,0);

%TXS=moving_average(TX,3,1);% 7-day (2*3+1) moving average for T
%FXS=moving_average(FX,1,1);% 3-day (2*1+1) moving average for F
%FXS(FXS<0)=0; 




% mortality as function of temperature, q10 function
q10=2.2; q10s = 0.8;
q10r = q10*ones(length(T),nrand)+q10s*ones(length(T),1)*randn(1,nrand);%randomize q10
q10r(q10r<1)=1; %limit q10.
q10r(q10r>4)=4; %limit q10

tbase=6;%base temperature
t_mort=q10r.^((TX-tbase)/10.);

mort_e_base = 0.40+0.08*ones(length(T),1)*randn(1,nrand);%compute base mortality
mort_n_base = 0.07+0.014*ones(length(T),1)*randn(1,nrand);
mort_c_base = 0.04+0.008*ones(length(T),1)*randn(1,nrand);
mort_a_base = 0.05+0.01*ones(length(T),1)*randn(1,nrand);

mort_e = mort_e_base.*t_mort;
mort_n = mort_n_base.*t_mort;
mort_c = mort_c_base.*t_mort;
mort_a = mort_a_base.*t_mort; 

% generate random number for egg production
mrand = randn(1,nrand);

dt = 1;% set temporal resolution to 1 day;

dcount = 0;

for dini=15:30:365    % initial time from 1 to 365   
   %disp(dini)
   dcount = dcount + 1;
   for k=1:nrand
       [h d1 d2 d3] = msd(dini,TX(:,k),FX(:,k),'cfin'); 
       hout(dcount,k)=h;
       d1out(dcount,k)=d1;
       d2out(dcount,k)=d2;
       d3out(dcount,k)=d3;
       d3tmp = d3;
       % c5-c6f, n_c5 normalized to 1
       dend1 = dini+d3;               % = first day of adulthood
       dur1 = [dini:dt:dend1];        % may need add additional days before c6f become reproductive
       dur1i = floor(dur1);   
       n_c6f=0.5*exp(-sum(mort_c(dur1i,k)*dt));   %nc6f, assume 50% become female

       %starting integration 
       c5a = 0;     
       c5a_ddd = 0;     
       for ddd = dend1:dt:500  %starting from the beginning of c6f; ended at day 500  (< 730)
           m = epr(floor(ddd),TX(:,k),FX(:,k),'cfin');
           ms=m;
           
           % c6fs die while producing
          dini2 = dend1;
          dend2 = ddd;
          dur2 = [dini2:dt:dend2];
          dur2i = floor(dur2);      
          [h d1 d2 d3] = msd(floor(dend2),TX(:,k),FX(:,k),'cfin'); 


          % eggs die while hatching
          dini3 = dend2;
          dend3 = ddd+h;
          dur3 = [dini3:dt:dend3];
          dur3i = floor(dur3);      
          [h d1 d2 d3] = msd(floor(dend3),TX(:,k),FX(:,k),'cfin'); 


          % Ns die while developing to C
          dini4 = dend3;
          dend4 = ddd+h+d1;
          dur4 = [dini4:dt:dend4];
          dur4i = floor(dur4);      
          [h d1 d2 d3] = msd(floor(dend4),TX(:,k),FX(:,k),'cfin'); 


          % Cs die while developing to C5
          dini5 = dend4;
          dend5 = ddd+h+d1+d2;
          dur5 = [dini5:dt:dend5];
          dur5i = floor(dur5);      

          tmp1 = exp(-sum(mort_a(dur2i,k)*dt));
          tmp2 = exp(-sum(mort_e(dur3i,k)*dt));
          tmp3 = exp(-sum(mort_n(dur4i,k)*dt));
          tmp4 = exp(-sum(mort_c(dur5i,k)*dt));

          c5atmp =  n_c6f.*tmp1.*ms.*tmp2.*tmp3.*tmp4;
          c5a_ddd_tmp = c5atmp.*(ddd-dend1);
          c5a=c5a + c5atmp;
          c5a_ddd = c5a_ddd + c5a_ddd_tmp;
       end
       c5a;
       tgen = c5a_ddd./c5a + d3tmp;    %average generation time
       p=[0:0.01:1];
       c5(k,1:length(p)) = c5a*p+(1-p)*(1-0.01).^tgen;   %consider diapause mortality over a generation time  
       tgen_out(dini,k) = tgen;
   end
   R = mean(log(c5));
   [C ID]=max(R);
   tmp2= (ID-1)/(length(p)-1);
   %pfit(dini) = tmp2;    %corresponding p value with max fitness R
   pfit(dcount) = tmp2;    %corresponding p value with max fitness R
end


   
   