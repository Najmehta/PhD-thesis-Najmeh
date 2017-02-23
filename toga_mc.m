function [var_index, errors,outinfo]=toga_mc(data, param, stability);
% TOGA: Using GA to select optimal variables based on Monto Carlo validation, ...
% the initial numVariables variables are pre-selected by UVE method, and the
%  sampling possibility of each variable was also determined by its  stability
% Syntax: [var_index, errors,outinfo]=toga_mc(data, param, stability);
%
%Input
%data:                      A structured data set including calibration set and its concentrations
%                              data.cal:      the calibration samples
%                              data.caltar:  the concentration vectors of calibratio samples
%param: the structured parameters for Genetic Algorithm, including three
%           parameters:
           %          param.factor: A vector, the optimal factor corresponding to different models
           %          param.numVariables: A vector, the number of variables preselected by UVE, 
           %                                           e.g. a signal with 400 data points was 
           %                                           preselected by UVE, and only 50 variables with
           %                                           largest stability were selected.
           %           param.numSelectedVars:A vector, the number of selected variables
           %           param.pretreatment: 1 with raw spectra, 2 with centerized spectra, 3 with SNV spectra,
           %                                            default is to use all of three methods
%stability:     a matrix, the UVE stability of each variable contributed to different models
%                  default is equal stability of each varialbe e.g.1
%
%Output
%var_index: the index number of selected variables
%errors: the predition errors corresponding to the model with selected
%           variables
%outinfo: the types of preatment, 1: raw, 2: center; 3: SNV
%
%By Da Chen, Dec 08, 2008
%Version 1.1
%Modified by Da Chen Feb 26, 2009
%Version 1.2

[yrow,ycol]=size(data.caltar);
[xrow,xcol]=size(data.cal);
factors=param.factor;
fname=fieldnames(param);
if sum(strcmp(fname,'folds'))==0
    folds=xrow;
else
    folds=param.folds;
end
if sum(strcmp(fname,'pretreatment'))==0
    pretreatment=4;
else
    pretreatment=param.pretreatment;
end

numVariables=param.numVariables;
numSelectedVars=param.numSelectedVars;
if nargin<3
    stability=ones(xrow,xcol);
end
% To define the list of random sample indecies
if 2.5*xrow<200
    sampleNum=100; % The number of samples
else
    sampleNum=3*xrow;
end
trainnum=ceil(0.618*xrow);
parfor i=1:sampleNum
          list=randperm(xrow);
          sampleIndices(i,:)=list(1:trainnum); %The random sample list of each run of Genetic algorithm
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:ycol
         
          [staValue,staIndex]=sort(-abs(stability(ii,:)));
          varIndex=staIndex(1:numVariables);
          varStability=staValue(1:numVariables(ii));
          cal=data.cal(:,varIndex); %The selected variables based on UVE method
          
          weight=varStability/sum(varStability); % To calculate the sampling possibility of each variable
          possibility=zeros(1,1+length(weight));
          for i=1:length(weight)
                possibility(i+1)=sum(weight(1:i));
          end
          
          options=gaoptimset('CreationFcn', {@varselect_uvegacreate, numVariables(ii),possibility}, ...
                                          'CrossoverFcn', @varselect_crossoverscattered,...
                                         'MutationFcn', @varselect_mutationuniform,...
                                         'PopulationSize',30,...   %
                                         'PopInitRange', [1;numVariables(ii)],...
                                         'FitnessLimit', 0,...
                                         'StallGenLimit',30,...    %
                                         'StallTimeLimit',2500,...   %
                                         'TimeLimit',3000,...          %
                                         'Generations',100,...          %
                                         'CrossoverFraction',0.7,...
                                         'Display','iter',...
                                         'PlotFcn',@gaplotbestf,...
                                         'UseParallel','never');

          if pretreatment==1
                  FitnessFcn1= {@varselect_gafit_mc_raw, sampleIndices, cal, data.caltar(:,ii), factors(ii)};
                  [selectedVars1, errorRate1]=ga(FitnessFcn1, numSelectedVars(ii), options);
                  var_index(ii,:)=sort(varIndex(selectedVars1));
                  errors(ii,:)=errorRate1;
                  outinfo(ii)=1;
          elseif pretreatment==2
                   FitnessFcn2= {@varselect_gafit_mc_cen, sampleIndices, cal, data.caltar(:,ii), factors(ii)};
                   [selectedVars2, errorRate2]=ga(FitnessFcn2, numSelectedVars(ii), options);
                   var_index(ii,:)=sort(varIndex(selectedVars2));
                   errors(ii,:)=errorRate2;
                   outinfo(ii)=2;
          elseif  pretreatment==3
                   FitnessFcn3= {@varselect_gafit_mc_snv, sampleIndices, cal, data.caltar(:,ii), factors(ii)};
                   [selectedVars3, errorRate3]=ga(FitnessFcn3, numSelectedVars(ii), options); 
                   var_index(ii,:)=sort(varIndex(selectedVars3));
                   errors(ii,:)=errorRate3;
                   outinfo(ii)=3;      
          else 
                  FitnessFcn1= {@varselect_gafit_mc_raw, sampleIndices, cal, data.caltar(:,ii), factors(ii)};
                  [selectedVars1, errorRate1]=ga(FitnessFcn1, numSelectedVars(ii), options);
                  FitnessFcn2= {@varselect_gafit_mc_cen, sampleIndices, cal, data.caltar(:,ii), factors(ii)};
                  [selectedVars2, errorRate2]=ga(FitnessFcn2, numSelectedVars(ii), options);
                  FitnessFcn3= {@varselect_gafit_mc_snv, sampleIndices, cal, data.caltar(:,ii), factors(ii)};
                  [selectedVars3, errorRate3]=ga(FitnessFcn3, numSelectedVars(ii), options);
                  errorlist=[errorRate1, errorRate2, errorRate3];
                  selectedVars=[selectedVars1; selectedVars2; selectedVars3];
                  [minerror,bb]=min(errorlist);
                  var_index(ii,:)=sort(varIndex(selectedVars(bb,:)));
                  errors(ii,:)=minerror;
                  outinfo(ii)=bb;

           end
       
           possibility=[];
end
