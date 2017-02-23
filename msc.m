function [mscData, reference] = msc(signals)
%Performs Multiplicatve Signal Correction (MSC) on data
%
% [mscData, offsetParam, slopeParam] = msc(data, offset, slope)
%
%Input Variables:
%  signals:        Data matrix to be transformed
%
%Output Variables:
%  mscData:     MSC transformed data matrix.
%  reference:   The reference spectrum used for prediction
%
%  The mean spectrum is chosen as reference spectrum for the
%  multiplicative model term
%
%  Author Da Chen, January 2nd, 2009

cal=signals;
[cala,calb]=size(cal);

mcal=mean(cal);
x=ones(1,length(mcal));
reference=[x;mcal];
reference=reference';
cbelta=inv(reference'*reference)*reference'*cal';
for i=1:cala
   mscData(i,:)=(cal(i,:)-cbelta(1,i)*ones(1,length(mcal)))/cbelta(2,i);
end

