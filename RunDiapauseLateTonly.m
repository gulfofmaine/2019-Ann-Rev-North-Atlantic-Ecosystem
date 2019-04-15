%RunDiapause

%1. Load the data
CH=load('NAchlFill.mat');%Chlorophyll
SST=load('NAsst.mat');%SST (in K)
load NApixels;%Ocean pixels

%2. Initialize the variables
[m,n,p]=size(CH.NAchl);
PDearly=nans(m,n,12);
PDlate=nans(m,n,12);

d=linspace(0,365,13);
d=[d(1:12)-365,d(1:12),d(1:12)+365];%days for interpolation


dfile=dir('/Users/apershing/Dropbox/Documents/Papers/01_InProgress/NAtlantic_Review/Work/PDlateTonly*.mat');
Nfile=length(dfile);
fprintf('\n\t\tSaving to file %d\n\n',Nfile);
fname=['/Users/apershing/Dropbox/Documents/Papers/01_InProgress/NAtlantic_Review/Work/PDlateTonly',padstr0(Nfile,2),'.mat'];
save(fname,'PDlate');

%Get T and C for this pixel
Irnd=randperm(length(Roc));


for j=1:length(Roc);
    if(mod(j,50)==0)
        fprintf('%8d\t%4.1f%%\n',j,100*j/length(Roc));
        save(fname,'PDlate');
    end
    
    jj=Roc(Irnd(j));%210;%
    kk=Coc(Irnd(j));%400;%
    if(jj>=75);%stay in the north
        C=squeeze(CH.NAchl(jj,kk,:));
        T=squeeze(SST.NAsst(jj,kk,:));

        %late period
        Cl=nanmean(reshape(C(1:60),12,5)');%ealry period
        Tl=nanmean(reshape(T(169:228),12,5)');%late period
        J=find(~isnan(Tl));
        J2=find(~isnan(Cl));
        if(length(J)>=9 & length(J2)>=6);
            %enough data to try the model
            Cl3=[Cl,Cl,Cl];
            Tl3=[Tl,Tl,Tl];
            I=find(~isnan(Cl3));
            ClI=interp1(d(I),Cl3(I),1:365,'pchip');%interpolate to daily
            I=find(~isnan(Tl3));
            TlI=interp1(d(I),Tl3(I),1:365,'pchip');%interpolate to daily
            pfit=JiDiapause(TlI(:)-273.15,ClI(:),10);%pfit will be monthly prb{diapause}
            PDlate(jj,kk,:)=pfit;
        end
    end
end

