function [wavedata, wcoefs,WL]=dwt_rebu(data,wfilter,scale,cutlevel);
%The function is used to remove the approximation and detail wavelet 
%components defined by users as baseline and noise, respecgively
%syntax:[wavedata, wcoefs]=dwt_rebu(data,wfilter,scale,cutlevel)
%input
%data: the signal to be treated
%wfilter: the wavelet filter.
%scale: the scale used in the wavelet transform, and correaponding
%approximate component is deleted.
%cutlevel:the level of detail component to be deleted.
%output
%wavedata: the wavelet filterted data after elimination of approximation
%and detail components.
%wcoef: the wavelet coefficient after elimination of approximation and
%detail components at defined scales.

if cutlevel>scale
    fprintf('the cutlevel can not be larger than scale!!\n');
    return;
end
[wa,wb]=size(data);
%Decompose the spectrum matrix into wavelet coefficient matrix using the
%function "dwt_matrix"
[wcal,wt_components,WL]=dwt_matrix(data,wfilter,scale);
wavedata=[];
for i=1:scale+1
    windex(i)=sum(WL(1:i));
end
for i=1:wa
    if cutlevel>0
        wavedata(i,:)=data(i,:)-wt_components{i}(1,:)-sum(wt_components{i}(end:-1:(end-cutlevel+1),:));
    else
        wavedata(i,:)=data(i,:)-wt_components{i}(1,:);
    end
end
wcoefs=wcal;
if cutlevel>0
    wcoefs(:,[1:windex(1) windex(scale+1):-1:windex(scale+1-cutlevel)+1])=0*wcoefs(:,[1:windex(1) windex(scale+1):-1:windex(scale+1-cutlevel)+1]);
else
    wcoefs(:,[1:windex(1)])=0*wcoefs(:,[1:windex(1)]);
end

        
