%% Example script showing how to simulate responses of the SOC model

%% Add code to the MATLAB path

setup

%% Load some stimuli

load('stimuli.mat','images');
stimulus = images(69+[16 70]);  % pick two stimuli
stimulus = cellfun(@(x) x(:,:,1),stimulus,'UniformOutput',0);  % use just the first frame
stimulus = cat(3,stimulus{:});

%% Compute the response of the SOC model

% initialize
cache = [];

% here, we compute the response of the SOC model, specifying certain values for the 
% parameters of the model.  the values used for the <sd>, <n>, and <c> parameters 
% are matched to the typical values found in V3.
[response,cache] = socmodel(stimulus,150,[0 254 127],1.2,{37.5 -1 1 8 2 .01 2 0}, ...
                            1,.5,1.5,1/1.5,.12,.99,cache);
%%

% here, we compute the responses again, but this time we set the <c> parameter to 0
[response2,cache] = socmodel(stimulus,150,[0 254 127],1.2,{37.5 -1 1 8 2 .01 2 0}, ...
                             1,.5,1.5,1/1.5,.12,  0,cache);
%%

% visualize the results
figure; setfigurepos([100 100 700 300]);
for p=1:2
  subplot(2,3,(p-1)*3+1);
  imagesc(stimulus(:,:,p),[0 254]); axis image tight; colormap(gray); colorbar; title('Stimulus');
  subplot(2,3,(p-1)*3+2);
  imagesc(response(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('Response (c=0.99)');
  subplot(2,3,(p-1)*3+3);
  imagesc(response2(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('Response (c=0)');
end
%%
