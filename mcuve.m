function  [S]=mcuve(cal,caltar,factor,runtimes, noise, trainnum)
%Calculate the stability of each normalized variables using Monte Carlo resampling
% [S]=mcuve(cal,caltar,factor,runtimes, noise, trainnum)
% The calibration set (cal) was normalized, thereafter, its stability was
% ... calculated.
%
%Input:
%cal: the training set;
%caltar: the concentrations of the training set;
%factor: the PLS factors used;
%runtimes: the number of resampling times
%noise: the artificial noise variables
%          0: means the noise variables are not added
%          1: means the noise variables are added.
%trainnum: the number of samples contained in the training subset
%
%Output
%S: the stability of each variable
%By Da Chen, Feb 26, 2009
%version 1.0

[cala,calb]=size(cal);

if noise==1
    if calb<400
        z=rand(cala,400)*10.^(-10); %To add the random matrix with very small value
    else
        z=rand(cala,calb)*10.^(10);
    end
    cal=[cal z];
end
if nargin<6
    trainnum=floor(0.6*cala);
end

%To calculate the different PLS coefficients by randomly altering the
%training set
coefs=zeros(runtimes,size(cal,2));
parfor i=1:runtimes
    list=randperm(cala);
    rlist=list(1:trainnum);
    newcal=cal(rlist,:);
    newtar=caltar(rlist);
    [newcal]=zscore(newcal); %To normalize the vector
    [newtar]=zscore(newtar); %To normalize the vector
    [coef] = plsxt(newcal, newtar, factor);
    coefs(i,:)=coef;
end

S=mean(coefs)./std(coefs);