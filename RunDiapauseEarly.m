%RunDiapause

%1. Load the data
CH=load('../../MonthlyChl/NAchlFill.mat');%Chlorophyll
SST=load('../../AVHRR_OI/NAsst.mat');%SST (in K)
load NApixels;%Ocean pixels

%2. Initialize the variables
[m,n,p]=size(CH.NAchl);
PDearly=nans(m,n,12);

d=linspace(0,365,13);
d=[d(1:12)-365,d(1:12),d(1:12)+365];%days for interpolation

%Get T and C for this pixel
for j=1:length(Roc);
    if(mod(j,50)==0)
        fprintf('%8d\t%4.1f%%\n',j,100*j/length(Roc));
        save PDearly PDearly
    end
    
    jj=Roc(j);%210;%
    kk=Coc(j);%400;%
    C=squeeze(CH.NAchl(jj,kk,:));
    T=squeeze(SST.NAsst(jj,kk,:));

    %early period
    Ce=nanmean(reshape(C(1:60),12,5)');%early period
    Te=nanmean(reshape(T(1:60),12,5)');%early period
    J=find(~isnan(Te));
    J2=find(~isnan(Ce));
    if(length(J)>=9 & length(J2)>=6);
        %enough data to try the model
        Ce3=[Ce,Ce,Ce];
        Te3=[Te,Te,Te];
        I=find(~isnan(Ce3));
        CeI=interp1(d(I),Ce3(I),1:365,'pchip');%interpolate to daily
        I=find(~isnan(Te3));
        TeI=interp1(d(I),Te3(I),1:365,'pchip');%interpolate to daily
        pfit=JiDiapause(TeI(:)-273.15,CeI(:),10);%pfit will be monthly prb{diapause}
        PDearly(jj,kk,:)=pfit;
    else
        PDearly(jj,kk,:)=0;
    end

end

