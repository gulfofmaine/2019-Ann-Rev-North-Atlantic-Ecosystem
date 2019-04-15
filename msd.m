function [h d1 d2 d3]=msd(ds,T,F,sp)
%calculate stage durations based on starting day (ds),
%time series of temperature (T), food (F) and species
 
%h: hatching time
%d1: n1-c1   (from beginning of n1 to begining of c1)
%d2: c1-c5
%d3: c5-c6

if (strcmp(sp,'cfin'))

   a = [595,387,582,1387,759,716,841,966,1137,1428,2166,4083];
   alpha = 9.11;
   beta= -2.05;
%   kfood = 1.0;
%   ffood = 1.0-exp(-F./kfood);   %independent from T
   
   %development time depend on both T and F (Ohman&Hseieh, JPR2008)
   %msd=a*(T+alpha)^(beta)/(1-exp(-q*T^s*F))
   %notice q, s calculated by fitting Campbell's data, where F is in mgC/m3
   %converted ugChl/l to mgC by a coefficient cc
   cc = 50.0;
   q=[0 0 0 0.2231    0.2276    0.2104    0.2414    0.1885    0.1798    0.2166    0.1803    0.2112];
   s=[0 0 0 -0.9523   -0.9206   -0.9740   -0.8968   -1.0306   -1.0460   -0.9689   -1.0430   -0.9643];
   
end  
 
%for h
   dd = ds;
   stg = 0;
   while stg<1
     tdd = T(dd);
     msd_h = a(1).*(tdd+alpha).^beta;
     stg = stg+1./msd_h;
     dd = dd+1;
   end
   h = dd-ds;

%for d1 (n1-c1)
   dd = ds;
   stg = 0;
   while stg<6
     tdd = T(dd);
 %    fdd = ffood(dd);  %need exclude egg,n1 and n2
     fdd = 1-exp(-q(2:7).*tdd.^s(2:7).*F(dd).*cc);    %fdd start from N1
                                          %if tdd<0, complex number could occur
     if((tdd)<=0);fdd(:)=1;end
     fdd(1:2) = 1;    %exclude n1 and n2
     msd_n = a(2:7).*(tdd+alpha).^beta;
     stg = stg+1./msd_n(floor(stg)+1).*fdd(floor(stg)+1);
     dd = dd+1;
   end
   d1 = dd-ds;

%for d2 (c1-c5)
   dd = ds;
   stg = 0;
   while stg<4
     tdd = T(dd);
     %fdd = ffood(dd);
     fdd = 1-exp(-q(8:11).*tdd.^s(8:11).*F(dd).*cc);
     if((tdd)<=0);fdd(:)=1;end
     msd_c = a(8:11).*(tdd+alpha).^beta;
     stg = stg+1./msd_c(floor(stg)+1).*fdd(floor(stg)+1);
     dd = dd+1;
   end
   d2 = dd-ds;

%for d3 (c5-c6)
   dd = ds;
   stg = 0;
   while stg<1
     tdd = T(dd);
     %fdd = ffood(dd);
     fdd = 1-exp(-q(12).*tdd.^s(12).*F(dd).*cc);
     if((tdd)<=0);fdd(:)=1;end
     msd_c5 = a(12).*(tdd+alpha).^beta;
     stg = stg+1./msd_c5.*fdd;
     dd = dd+1;
   end
   d3 = dd-ds;

