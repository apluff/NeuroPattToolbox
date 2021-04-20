function outputs = mainProcessingWithOutput(data, Fs, params, textboxHandle, onlyPatterns, suppressFigures)
% Function to filter data, calculate velocity vector fields and identify
% patterns in (row)x(column)x(time)x(repetitions) matrix DATA with sampling
% frequency FS Hz, using the parameter structure PARAMS. Outputs progress
% to the command line, also outputs to text box with handle TEXTBOXHANDLE
% if given as an input. If flag ONLYPATTERNS is true, filtered data and
% velocity fields will not be output. If flag SUPPRESSFIGURES is true, no
% new figures will be generated.

% Rory Townsend, Aug 2018
% rory.townsend@sydney.edu.au

if nargin < 4
    textboxHandle = [];
end
if nargin < 5
    onlyPatterns = false;
end
if nargin < 6
    suppressFigures = false;
end

%% Pre-processing
startTime = datetime;
outputProgress(sprintf('Beginning NeuroPatt pattern analysis at %s', ...
    datestr(startTime)), textboxHandle)
% Optionally subtract the baseline or z-score all channels
if params.zscoreChannels || params.subtractBaseline
    data = bsxfun(@minus, data, mean(data,3));
end
if params.zscoreChannels
    data = bsxfun(@rdivide, data, std(data,[],3));
end

% Find any invalid channels that must be interpolated over
% Any channels with any NaN values
nanChans = any(isnan(data(:,:,:)),3);
% Any channels that never change over time
zeroChans = all(data(:,:,:)==0, 3);
badChannels = find(nanChans | zeroChans);

%% Band-pass filter data
outputProgress('Filtering waveforms...', textboxHandle)
tic
% Time should always be in the third dimension, but keep this as a
% parameter for possible future changes
timeDim = 3;
ntimesteps = size(data, timeDim) - 1;

% Use band-pass filter then Hilbert transform or Morlet wavelet transform
if params.useHilbert
    wvcfs = filterSignal(data,params.hilbFreqLow, params.hilbFreqHigh, ...
        Fs, timeDim, true);

else
    % Morlet wavelets
    wvcfs = squeeze(morletWaveletTransform(data, Fs, params.morletCfreq, ...
        params.morletParam, timeDim));
end

toc

%% Calculate velocity vector fields using optical flow estimation
outputProgress('Calculating velocity vector fields...', textboxHandle)
tic

% Pre-allocate velocity field variables
vfs = zeros(size(wvcfs));
vfs = vfs(:,:,1:end-1,:);
meanCSteps = zeros(size(wvcfs,4), 1);

% Calculate velocity fields for every trial, and same average number of
% steps to converge
for itrial = 1:size(wvcfs,4)
    [vx, vy, csteps] = opticalFlow2(wvcfs(:,:,:,itrial), badChannels, ...
        params.opAlpha, params.opBeta, ~params.useAmplitude);
    % alternatively, you can use Matlab built-in function 'opticalflowHS', 
    % which does not provide the penalty term Beta but is optimised for speed.
    % alpha = 1 ;
    % opticFlow = opticalFlowHS('Smoothness',alpha,'MaxIteration',100) ;
    % Vx = zeros(size(wvcfs,1),size(wvcfs,2),size(wvcfs,3)-1) ;
    % Vy = zeros(size(wvcfs,1),size(wvcfs,2),size(wvcfs,3)-1) ;
    %  for iTime = 1:downR:size(sigIn,3)
    %        flow = estimateFlow(opticFlow,wvcfs(:,:,iTime,itrial));
    %        Vx(:,:,iTime) = flow.Vx ;
    %        Vy(:,:,iTime) = flow.Vy ;

    %  end

    vfs(:,:,:,itrial) = vx + 1i*vy;

    meanCSteps(itrial) = mean(csteps);
    fprintf('Processed trial %i\n', itrial)
end

toc
outputProgress(...
    sprintf('Optical flow took %0.1f steps on average to converge.\n', ...
    mean(meanCSteps)), textboxHandle)

%% Perform singular value decomposition of velocity fields
if params.performSVD && ~suppressFigures
    outputProgress('Performing SVD of velocity vector fields...', textboxHandle)
    tic
    % Open new figure and plot SVD modes
    figure('Name', 'Dominant SVD modes')
    plotTime = (1:ntimesteps)/Fs;
    plotcsvd(vfs, params.nSVDmodes, plotTime, params.useComplexSVD);
    toc
end

%% Find all patterns present
outputProgress('Identifying all patterns in velocity fields...', textboxHandle)
tic
% Set up pattern structures
allPatts = cell(1, size(wvcfs,4));
allLocs = allPatts;
% Loop over all repetitions to find patterns
for itrial = 1:size(wvcfs,4)
    thisvf = vfs(:,:,:,itrial);
    [patterns, patternTypes, patternColNames, pattLocs] = ...
        findAllPatterns(real(thisvf), imag(thisvf), params, ...
        angle(wvcfs(:,:,:,itrial)));
    allPatts{itrial} = patterns;
    allLocs{itrial} = pattLocs;
end

toc

%% Analyse evolution between patterns
% Number of time steps before and after a pattern ends to search for other
% patterns
nafter = round(0.05*Fs);
nbefore = round(0.01*Fs);

pattTypeStr = '';
for itype = 1:length(patternTypes)
    pattTypeStr = sprintf('%s%i.%s ', pattTypeStr,itype,patternTypes{itype});
end

[nobs, nexp] = pattEvolution(allPatts, ntimesteps, nafter, nbefore);
rateDiff = (nobs - nexp) / ntimesteps * Fs;
%disp('Observed minus expected pattern transitions/sec')
disp(pattTypeStr)
%disp(nanmean(rateDiff,3))
%disp(median(nobs,3) - median(nexp,3))

disp('Fractional change between observed and expected')
disp(nanmean((nobs-nexp)./nexp, 3));

% Test differences between observed and expected if multiple trials are
% present
if size(vfs, 4) > 1

    pvals = zeros(size(nobs,1));
    for initPatt = 1:size(nobs,1)
        for nextPatt = 1:size(nobs,2)
            thisObs = nobs(initPatt, nextPatt, :);
            thisExp = nexp(initPatt, nextPatt, :);
            [h, p] = ttest(thisObs(:),  thisExp(:));
            pvals(initPatt, nextPatt) = p;
        end
    end
    disp('Paired t-test p-values')
    fprintf('Bonferroni correction factor = %i\n', numel(pvals))
    disp(pvals)
end

% Set all outputs
if ~onlyPatterns
    outputs.filteredSignal = wvcfs;
    outputs.velocityFields = vfs;
end
outputs.badChannels = badChannels;
outputs.nTimeSteps = ntimesteps;
outputs.patternTypes = patternTypes;
outputs.patternResultColumns = patternColNames;
outputs.patterns = allPatts;
outputs.patternLocs = allLocs;
outputs.params = params;
outputs.pattTransitionsObs = nobs;
outputs.pattTransitionsExp = nexp;
outputs.Fs = Fs;
outputs.processTime = datetime - startTime;

function outputProgress(outputStr, textboxHandle)
% Outputs the current progress given by TEXTSTR to the terminal and also
% appends text to new line in textbox given by TEXTBOXHANDLE
disp(outputStr)
if ishandle(textboxHandle)
    currentStr = get(textboxHandle, 'String');
    set(textboxHandle, 'String', sprintf('%s\n%s', currentStr{1}, outputStr));
    drawnow
end