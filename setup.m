% This script downloads the example dataset (if it has not already
% been downloaded) and adds analyzePRF to the MATLAB path.

% Download the dataset if it has not already been downloaded.
% (If you like, you may manually download the dataset from:
%    http://kendrickkay.net/analyzePRF/exampledataset.mat)
files = {'exampledataset.mat' 'socmodel.zip'};
source = {'http://kendrickkay.net/analyzePRF/exampledataset.mat' ...
    'http://kendrickkay.net/socmodel/socmodel.zip'};

for p=1:length(files)
  if ~exist(files{p},'file')
    fprintf('Downloading %s (please be patient).\n',files{p});
    urlwrite(source{p},files{p});
    fprintf('Downloading is done!\n');
    [~, ext] = fileparts(files{p});
    if strcmp(ext, 'zip'), unzip(files{p});  end
  end
end
clear p files;

% Add analyzePRF to the MATLAB path (in case the user has not already done so).
path0 = strrep(which('setup.m'),'/setup.m','');
addpath(genpath(path0));
clear path0;
