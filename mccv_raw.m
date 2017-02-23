function [tarcv, tarv, factor, param]=mccv_raw(cal, caltar, val, runtimes, trainnum)
%
% The function is used to construct the Monte Carlo cross-validation with Raw spectra
% meanwhile, the validation set can be predicted by ensemble MC model
% [tarcv, tarv, press, param]=mccv_raw(cal, caltar, val, runtimes, trainnum)
%
%Input:
%cal: calibration set;
%caltar: the concentrations of calibration set
% val: validation set;
%runtimes: the number of runs of random resampling, the default is 2.5*sample number
%trainnum: the number of samples in the training subset, the default is 0.6*sample number
%
%Output:
%tarcv: the cross-validation prediction results obtained by Monte Carlo
%tarv: the predicted values obtained by ensemble MC models
%factor: the number of PLS factor determined by Monte Carlo cross-validation
%param: the structured parameters of each member model of ensemble PLS
%            param.press: the press of the MCCV validation
%            param.samindex: the selected samples' indexes in each PLS model
%            param.coefs: the PLS regression coeffficients of each PLS member model
%            param.tarvs: the predicted value vectors obtained with different PLS model
%
% By Da Chen 15/08/2008
% Modified by Da Chen, Aug. 22, 2010
% Da Chen, Aug. 23, 2011

maxfactor=15;
[cala,calb]=size(cal);
if nargin<4
    if cala<50
       runtimes=5*cala;
    else
        runtimes=2.5*cala;
    end
end
if nargin<5
    if cala<30
        trainnum=floor(0.8*cala);
    else
        trainnum=floor(0.6*cala);
    end
end
valnum=cala-trainnum;

%To calculate the PRESS of multiple PLS models
samindex=zeros(runtimes,cala);
presses=zeros(runtimes,maxfactor);
coefs=zeros(runtimes,calb);
tarvs3=zeros(runtimes,maxfactor,valnum);

for i=1:runtimes
    list=randperm(cala);
    samindex(i,:)=list;
    newcal=cal(list(1:trainnum),:);
    newtar=caltar(list(1:trainnum));
    newval=cal(list((trainnum+1):cala),:);
    newvaltar=caltar(list((trainnum+1):cala));   
    [p,q,w,b] = pls1(newcal', newtar', maxfactor); 
    for f=1:maxfactor
          tar_v=plspred(newval', p, q, w, b,f);
          tarvs3(i,f,:)=tar_v';
          presses(i,f)=sum((tar_v'-newvaltar).^2);
    end
end

press=mean(presses);
ratio=press./min(press);
factor=min(find(ratio<1.25));
tarvs=tarvs3(:,factor,:);

calindex=samindex(:,1:trainnum);
valindex=samindex(:,trainnum+1:cala);
Ftarvs=zeros(runtimes,size(val,1));
parfor i=1:runtimes
     newcal=cal(calindex(i,:),:);
     newtar=caltar(calindex(i,:));
     coefs(i,:)=plsxt(newcal,newtar,factor);
     Ftarvs(i,:)=(val*coefs(i,:)')';
end    
parfor i=1:cala
    sam_index=find(valindex==i);
    tarcv(i)=mean(tarvs(sam_index));
    sam_index=0;
end    
tarcv=tarcv';
tarv=(mean(Ftarvs))';
param.tarvs=Ftarvs;
param.coefs=coefs;
param.press=press;
param.samindex=samindex;