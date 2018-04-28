%% Example script illustrating how to fit the CSS model
%% Add code to the MATLAB path
%%
clear; close all; setup;
%% specify the index of the voxel to fit
%%
ix = 774; % same voxel used for pRF_model_example_1_linear
ix = 1945; % TO1
%ix = 1996; % V3AB
%ix = 2036; % V3AB
%% Load data
%%
% load data from the first dataset
load('dataset01.mat','betamn','betase');
%% Load stimuli
%%
load('stimuli.mat','conimages');
%% Perform stimulus pre-processing
%%
% extract the stimuli we need and then concatenate along the third dimension
stimulus = conimages(1:69);
stimulus = cat(3,stimulus{:});

% resize the stimuli to 100 x 100 (to reduce computational time)
temp = zeros(100,100,size(stimulus,3));
for p=1:size(stimulus,3)
  temp(:,:,p) = imresize(stimulus(:,:,p),[100 100],'cubic');
end
stimulus = temp;

% ensure that all values are between 0 and 1
stimulus(stimulus < 0) = 0;
stimulus(stimulus > 1) = 1;

% inspect one of the stimuli
figure;
imagesc(stimulus(:,:,10));
axis image tight;
colormap(gray);
colorbar;
title('Stimulus');
%%
% reshape stimuli into a "flattened" format: 69 stimuli x 100*100 positions
stimulus = reshape(stimulus,100*100,69)';
%% Prepare for model fitting
%%
% to perform model fitting, we will be using fitnonlinearmodel.m.  this function
% is essentially a wrapper around MATLAB's lsqcurvefit.m function.  the benefit
% of fitnonlinearmodel.m is that it simplifies input and output issues, deals with
% resampling (cross-validation and bootstrapping), makes it easy to evaluate multiple
% initial seeds, and makes it easy to perform stepwise fitting of models.
%
% to prepare for the call to fitnonlinearmodel.m, we have to define various 
% input parameters.  this is what we will now do.

% define constants
res = 100;  % resolution of the pre-processed stimuli

% the parameters of the CSS model are [R C S G N] where
%   R is the row index of the center of the 2D Gaussian
%   C is the column index of the center of the 2D Gaussian
%   S is the standard deviation of the 2D Gaussian
%   G is a gain parameter
%   N is the exponent of the power-law nonlinearity

% define the initial seed for the model parameters
seed = [(1+res)/2 (1+res)/2 res 1 0.5];

% define bounds for the model parameters
bounds = [1-res+1 1-res+1 0   -Inf 0;
          2*res-1 2*res-1 Inf  Inf Inf];

% fitnonlinearmodel.m provides the capacity to perform stepwise fitting.
% here, we define a version of bounds where we insert a NaN in the first
% row in the spot that corresponds to the exponent parameter.  this 
% indicates to fix the exponent parameter and not optimize it.
boundsFIX = bounds;
boundsFIX(1,5) = NaN;

% issue a dummy call to makegaussian2d.m to pre-compute xx and yy.
% these variables are re-used to achieve faster computation.
[d,xx,yy] = makegaussian2d(res,2,2,2,2);

% we will now define a function that implements the CSS model.  this function
% accepts a set of parameters (pp, a vector of size 1 x 5) and a set of stimuli 
% (dd, a matrix of size A x 100*100) and outputs the predicted response to those 
% stimuli (as a vector of size A x 1).  for compactness, we implement this 
% function as an anonymous function where the parameters are given by pp
% and the stimuli are given by dd.
modelfunCSS       = @(pp,dd) pp(4)*((dd*vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),pp(3),xx,yy,0,0)/(2*pi*pp(3)^2))).^pp(5));

modelfunLinear    = @(pp,dd) pp(4)*(dd*vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),pp(3),xx,yy,0,0)/(2*pi*pp(3)^2)));


% notice that the overall structure of the model is 
%   RESP = GAIN*(STIM*GAU).^N
% where STIM*GAU represents the dot product between the stimulus and the 2D Gaussian.
% also, note that the division by (2*pi*pp(3)^2) makes it such that the integral
% of the Gaussian is equal to 1 (this aids the interpretation of model parameters).

% now that we have defined modelfun, we are ready to define the final model
% specification.  in the following, we specify a stepwise fitting scheme.
% in the first fit (the first row), we start at the seed and optimize all 
% parameters except the exponent parameter.  in the second fit (the second row),
% we start at the parameters estimated in the first fit and optimize all parameters.
% the purpose of the stepwise fitting is to help converge to a good solution 
% (i.e. avoid local minima).  (the anonymous functions in the second row accept
% a single input, ss, which refers to the parameters estimated in the first fit.)
modelCSS = {{seed       boundsFIX   modelfunCSS} ...
         {@(ss) ss   bounds      @(ss) modelfunCSS}};

modelLinear = {seed(1:4) bounds(:,1:4)  modelfunLinear};

% define the resampling scheme to use.  here, we use 0, which
% means to just fit the data (no cross-validation nor bootstrapping).
resampling = 0;

% define the metric that we will use to quantify goodness-of-fit.
% here, we use a version of the coefficient of determination (R^2)
% in which variance in the data is computed relative to 0.
% this is sensible since the data being fitted are beta weights that
% represent evoked BOLD responses relative to the baseline signal level
% (which corresponds to 0).
metric = @(a,b) calccod(a,b,[],[],0);


% finally, construct the options struct that will be passed to fitnonlinearmodel.m
optCSS = struct( ...
  'stimulus',    stimulus, ...
  'data',        betamn(ix,:)', ...
  'model',       {modelCSS}, ...
  'resampling',  resampling, ...
  'metric',      metric);

optLinear = struct( ...
  'stimulus',    stimulus, ...
  'data',        betamn(ix,:)', ...
  'model',       {modelLinear}, ...
  'resampling',  resampling, ...
  'metric',      metric);

% do a quick inspection of opt
optCSS
%% 
% 
%% Fit the model
%%
resultsCSS    = fitnonlinearmodel(optCSS);
resultsLinear = fitnonlinearmodel(optLinear);
%% 
% 
%% Inspect the results
%%
% these are the final parameter estimates
resultsCSS.params
resultsLinear.params
%%
% this is the R^2 between the model fit and the data
resultsCSS.trainperformance

% compare to the linear model fit
resultsLinear.trainperformance
%%
% visualize the parameter estimates
figure; hold on;
pp = resultsCSS.params;
  % draw a circle indicating the PRF location +/- 2 PRF sizes.
  % PRF size is defined as S/sqrt(N).
drawellipse(pp(2),pp(1),0,2*pp(3)/sqrt(pp(5)),2*pp(3)/sqrt(pp(5)),[],[],'k');
drawrectangle((1+res)/2,(1+res)/2,res,res,'k-');
plot([.5 res+.5], (res+1)/2*[1 1], 'k--', (res+1)/2*[1 1], [.5 res+.5] , 'k--')
axis([.5 res+.5 .5 res+.5]);
set(gca,'YDir','reverse');
axis square;
title('Estimated PRF location and size');

pp = resultsLinear.params;
  % draw a circle indicating the PRF location +/- 2 PRF sizes.
  % PRF size is defined as S/sqrt(N).
drawellipse(pp(2),pp(1),0,2*pp(3),2*pp(3),[],[],'r');
%%
% visualize the data and the model fit
figure; setfigurepos([100 100 450 250]); hold on;
bar(betamn(ix,:),1);
errorbar2(1:69,betamn(ix,:),betase(ix,:),'v','g-','LineWidth',2);
modelfit = modelfunCSS(resultsCSS.params,stimulus);
plot(1:69,modelfit,'k-','LineWidth',3);
ax = axis;
axis([0 70 ax(3:4)]);
xlabel('Stimulus number');
ylabel('BOLD signal (% change)');
title('Data and model fit');

modelfitLinear = modelfunLinear(resultsLinear.params,stimulus);
plot(1:69,modelfitLinear,'r-','LineWidth',2);
%% 
% 
%% Try a different resampling scheme: bootstrapping
%%
% define an options struct that specifies 25 bootstraps (i.e. draw 
% with replacement from the 69 data points5
optBOOT = optCSS;
optBOOT.resampling = {50 0};  % the 0 sets the random-number seed to 0
optBOOT.optimoptions = {'Display' 'off'};  % turn off reporting

% fit the model
resultsBOOTCSS = fitnonlinearmodel(optBOOT);

optBOOT = optLinear;
optBOOT.resampling = {50 0};  % the 0 sets the random-number seed to 0
optBOOT.optimoptions = {'Display' 'off'};  % turn off reporting

% fit the model
resultsBOOTLinear = fitnonlinearmodel(optBOOT);
%%
% visualize the parameter estimates
figure; hold on;
subplot(1,2,1)
for p=1:size(resultsBOOTCSS.params,1)
  pp = resultsBOOTCSS.params(p,:);
  h = drawellipse(pp(2),pp(1),0,2*pp(3)/sqrt(pp(5)),2*pp(3)/sqrt(pp(5)));
  set(h,'Color',rand(1,3));
end
plot([.5 res+.5], (res+1)/2*[1 1], 'k--', (res+1)/2*[1 1], [.5 res+.5] , 'k--')
drawrectangle((1+res)/2,(1+res)/2,res,res,'k-');
axis([.5 res+.5 .5 res+.5]);
set(gca,'YDir','reverse');
axis square;
title('Bootstrap results CSS');

% visualize the parameter estimates
subplot(1,2,2); hold on;
for p=1:size(resultsBOOTLinear.params,1)
  pp = resultsBOOTLinear.params(p,:);
  h = drawellipse(pp(2),pp(1),0,2*pp(3),2*pp(3));
  set(h,'Color',rand(1,3));
end
plot([.5 res+.5], (res+1)/2*[1 1], 'k--', (res+1)/2*[1 1], [.5 res+.5] , 'k--')
drawrectangle((1+res)/2,(1+res)/2,res,res,'k-');
axis([.5 res+.5 .5 res+.5]);
set(gca,'YDir','reverse');
axis square;
title('Bootstrap Linear');
%% 
% 
%% Example of how to simulate model responses
%%
% let's take the model fitted to the full dataset and compute
% the predicted response of the model to some new stimuli.

% let's compute the response of the model to a point stimulus
% that is positioned at different locations in the visual field.
resp = zeros(res,res);
for r=1:res
  for c=1:res
    stim0 = zeros(res,res);
    stim0(r,c) = 1;
    resp(r,c) = modelfunCSS(resultsCSS.params,flatten(stim0));
  end
end

% visualize the results
figure;
mx = max(resp(:));
imagesc(resp,[0 mx]);
axis image tight;
colormap(gray);
title('Predicted response to point stimuli');
pp = resultsCSS.params;
drawellipse(pp(2),pp(1),0,pp(3),pp(3), [], [], 'r-');
drawellipse(pp(2),pp(1),0,2*pp(3),2*pp(3), [], [], 'c--');
%% notice that the results are consistent with the definition of PRF size as S/sqrt(N)
%%
% visualize the results
figure;
mx = max(resp(:));
imagesc(resp,[0 mx]);
axis image tight;
colormap(gray);
title('Predicted response to point stimuli');
pp = resultsCSS.params;
drawellipse(pp(2),pp(1),0,pp(3)/sqrt(pp(5)),pp(3)/sqrt(pp(5)), [], [], 'r-');
drawellipse(pp(2),pp(1),0,2*pp(3)/sqrt(pp(5)),2*pp(3)/sqrt(pp(5)), [], [], 'c--');