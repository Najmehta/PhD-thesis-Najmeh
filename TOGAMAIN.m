[wavedata, wcoefs,WL]=dwt_rebu(Spectra,'sym5',6,2)
[S]=mcuve(wcoefs,Con,7,500, 0,20);
data.cal=wcoefs;
data.caltar=Con;
param.factor=7;
param.numVariables=200;
param.numSelectedVars=10;
param.pretreatment=1;
for i=1:10
[var_indexi, errors,outinfo]=toga_mc(data, param, S);
var_index(i,:)=var_indexi;
end
for i=1:length(wcoefs)
fre(i)=length(find(var_index==i));
end
[a,b]=sort(-abs(fre));
var_index2=b(1:10);
newwc=wcoefs*0;
newwc(:,var_index2)=wcoefs(:,var_index2);
[wt_recon]=dwt_reconstructed(newwc, WL, 'sym5', 6);

