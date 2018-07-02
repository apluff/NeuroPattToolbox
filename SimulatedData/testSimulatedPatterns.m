% Script to test the accuracy of pattern detection and tracking with
% multiple simulated critical point patterns

close all

%% Set parameters
gridSize = [16 16 10];
alphas = [0.5];
%alphas = linspace(0.1, 1.5, 8);
betas = [0.01, 10];
noiselevels = logspace(-2, 1, 7);
%noiselevels = 0.5;
nreps = 20;

% Set parameters that apply to both patterns
freq = 0.1;
wavelength = 5;
gaussWidth = [];
minAmplitude = 1;
maxAmplitude = 2;
maxVelocity = 0.1;
minWidth = 3;
maxWidth = 5;

% Set details of first pattern
ampA = [];
typeA = 'source';
fullNameA = 'Source';
locA = [];
velA = [];
% ampA = 1;
% fullNameA = 'unstableNode';
% locA = [0.4,0.3] .* gridSize(1:2);
% velA = [-0.1 0];

% Set details of second pattern
ampB = [];
typeB = 'sink';
fullNameB = 'Sink';
locB = [];
velB = [];
% ampB = 1.5;
% fullNameB = 'stableNode';
% locB = [0.8,0.6] .* gridSize(1:2);
% velB = [0 0.05];

% Check if data needs to be generated each iteration
if strcmp(typeA, 'random') || strcmp(typeB, 'random') || ...
        isempty(locA) || isempty(locB) || isempty(velA) || isempty(velB)
    staticData = false;
else
    staticData = true;
end

% Pattern detection parameters
params.maxTimeGap = 0;
params.minCritRadius = 2;
params.minEdgeDist = 2;
params.minDuration = 0;
params.combineNodeFocus = true;

%% Set up output variables
meanErrors = zeros(length(alphas), length(betas), ...
    length(noiselevels), nreps);
semErrors = meanErrors;
meanIterations = meanErrors;
missingPatterns = meanErrors;
extraPatterns = meanErrors;
timeStepsWithNoExtra = meanErrors;

%% Generate both patterns and combine
if staticData
    waveA = generatePattern(gridSize, typeA, ampA, freq, wavelength, ...
        locA, velA, gaussWidth);
    waveB = generatePattern(gridSize, typeB, ampB, freq, wavelength, ...
        locB, velB, gaussWidth);
    wave = waveA + waveB;
end

%% Loop over all parameter values and noise levels
tic
plotVfs = false;
for ialpha = 1:length(alphas)
    fprintf('Processing waves with alpha = %s\n', string(alphas(ialpha)))
    
    for ibeta = 1:length(betas)
        rng(1)
        for inoise = 1:length(noiselevels)
            for irep = 1:nreps
                %% Generate random patterns
                % Generate random amplitudes
                if isempty(ampA)
                    iampA = minAmplitude + rand * (maxAmplitude - minAmplitude);
                else
                    iampA = ampA;
                end
                if isempty(ampB)
                    iampB = minAmplitude + rand * (maxAmplitude - minAmplitude);
                else
                    iampB = ampB;
                end
                
                % Generate random center locations that are not too
                % close to the edge and not too close together
                if isempty(locA) || isempty(locB)
                    ilocA = params.minEdgeDist + rand(1,2) .* ...
                        (gridSize(1:2) - 2*params.minEdgeDist);
                    locisvalid = false;
                    for ii=1:1000
                        ilocB = params.minEdgeDist + rand(1,2) .* ...
                            (gridSize(1:2) - 2*params.minEdgeDist);
                        if sqrt(sum((ilocB-ilocA).^2)) > 2*params.minCritRadius
                            locisvalid = true;
                            break
                        end
                    end
                    if ~locisvalid
                        error('Could not find valid critical point location! Consider reducing minCritRadius or minEdgeDistance.')
                    end
                else
                    ilocA = locA;
                    ilocB = locB;
                end
                % Generate random velocities with magnitude less than
                % maxVelocity parameter
                if isempty(velA)
                    vrand = rand * maxVelocity;
                    ivelA = [vrand, rand * sqrt(maxVelocity^2 - vrand^2)];
                else
                    ivelA = velA;
                end
                if isempty(velB)
                    vrand = rand * maxVelocity;
                    ivelB = [vrand, rand * sqrt(maxVelocity^2 - vrand^2)];
                else
                    ivelB = velB;
                end
                % Generate random Gaussian width
                if isempty(gaussWidth)
                    igausswidth = minWidth + rand * (maxWidth - minWidth);
                else
                    igausswidth = gaussWidth;
                end
                
                if ~staticData
                    % Generate random wave
                    [waveA, itypeA] = generatePattern(gridSize, typeA, ...
                        iampA, freq, wavelength, ilocA, ivelA, igausswidth);
                    [waveB, itypeB] = generatePattern(gridSize, typeB, ...
                        iampB, freq, wavelength, ilocB, ivelB, igausswidth);
                    wave = waveA + waveB;
                end
                
                %% Process data
                % Add noise proportional to total amplitude
                waveStd = std(real(wave),[],3);
                % If no wave is present, just add random noise
                if all(isnan(waveStd(:))) || all(waveStd(:)==0)
                    noise = randn(size(wave));
                    thisWave = noise;
                else
                    noise = bsxfun(@times, randn(size(wave)), ...
                        waveStd * noiselevels(inoise));
                    thisWave = real(wave) + noise;
                end
                
                % Band-pass filter data
                thisWave = squeeze(morletWaveletTransform(...
                    thisWave, 1, freq, 6, 3));
                
                % Make velocity field
                [vx, vy, nit] = opticalFlow(thisWave, [], ...
                    alphas(ialpha), betas(ibeta), 1);
                
                % Find all detected patterns
                [patterns, pattTypes, colNames, allPatternLocs, params] ...
                    = findAllPatterns(vx, vy, params);
                allPatternLocs = cat(1, allPatternLocs{:});
                
                %% Plot results
                % Plot example phase maps if noise level doesn't
                % change and there aren't too many parameter values
                if all([ialpha,ibeta,inoise,irep]==1) && ...
                        length(noiselevels) == 1 && ...
                        length(alphas)*length(betas) <= 6
                    plotVfs = true;
                    figure
                    for iplot = 1:2
                        subplot(length(betas), length(alphas)+1, ...
                            (iplot-1) * (length(alphas) + 1) + 1)
                        imagesc(angle(thisWave(:,:,2*iplot)), [-pi pi])
                        colormap(pmkmp_new)
                        axis xy
                    end
%                     plotSnapshotsAndVfs(angle(thisWave), vx, vy, ...
%                         allPatternLocs, 1, 1, 5)
                end
                
                % Also plot velocity fields for each parameter value
                if plotVfs
                    % Make quiver plots of velocity field
                    subplot(length(betas), length(alphas)+1, ...
                        (ibeta-1) * (length(alphas) + 1) + ialpha + 1)
                    quiver(vx(:,:,2), vy(:,:,2))
                    xlim([-0.5, size(vx,2)+1.5])
                    ylim([-0.5, size(vx,1)+1.5])
                    xticks([])
                    yticks([])
                    % Add critical point locations
                    hold on
                    thisTime = allPatternLocs(:,3) == 2;
                    thisLoc = allPatternLocs(thisTime, 1:2);
                    scatter(thisLoc(:,2), thisLoc(:,1), 'filled')
                    hold off
                    % Add title
                    title(sprintf('Alpha=%s, beta=%s', ...
                        string(alphas(ialpha)), string(betas(ibeta))))
                    drawnow
                end
                
                %% Test accuracy of critical point centres
                typeInds = zeros(1,2);
                if staticData
                    typeInds(1) = find(strcmp(pattTypes, fullNameA));
                    typeInds(2) = find(strcmp(pattTypes, fullNameB));
                end
                pattDists = nan(gridSize(3)-1, 2);
                for jtime = 1:gridSize(3)-1
                    % Detected critical points at this time step
                    jpatts = allPatternLocs(allPatternLocs(:,3)==jtime, :);
                    % True positions at this time step
                    jlocA = ilocA + jtime * ivelA;
                    jlocB = ilocB + jtime * ivelB;
                    jlocs = [jlocA; jlocB];
                    % Check both critical point types
                    for itype = 1:2
                        if staticData
                            thisPatt = jpatts(jpatts(:,4)==typeInds(itype), :);
                        else
                            thisPatt = jpatts;
                        end
                        if ~isempty(thisPatt)
                            thisLoc = repmat(jlocs(itype,:), ...
                                size(thisPatt,1), 1);
                            % Calculate distance between detected and true
                            minDist = min(sqrt(sum((thisPatt(:,[2,1]) - ...
                                thisLoc).^2, 2)));
                            if minDist < 2
                                pattDists(jtime, itype) = minDist;
                            end
                        end
                    end
                end
                
                % Save results to output arrays
                meanErrors(ialpha,ibeta,inoise,irep) = ...
                    nanmean(pattDists(:));
                semErrors(ialpha,ibeta,inoise,irep) = ...
                    nanstd(pattDists(:)) / sum(~isnan(pattDists(:)));
                meanIterations(ialpha,ibeta,inoise,irep) = mean(nit);
                missingPatterns(ialpha,ibeta,inoise,irep) = ...
                    sum(isnan(pattDists(:)));
                nextra = max([0, size(allPatternLocs,1) - ...
                    sum(~isnan(pattDists(:)))]);
                extraPatterns(ialpha,ibeta,inoise,irep) = nextra;
                %timeStepsWithNoExtra(ialpha,ibeta,inoise,irep) = ...
                %    gridSize(3)-1-length(unique(allPatternLocs(:,3)));
            end
        end
    end
end
toc

%% Plot results
figure
plotOptions = {'-o', 'CapSize', 0, 'MarkerFaceColor', 'auto', ...
    'MarkerEdgeColor', 'auto', 'MarkerSize', 4};
if length(noiselevels) > 1
    % Plot noise level as x-axis
    xdata = repmat(log10(noiselevels), length(alphas), 1)';
    xlab = 'Log noise level';
    plotData = @(x, ibeta) errorbar(xdata, ...
        squeeze(nanmean(x(:,ibeta,:,:),4))', ...
        squeeze(nanstd(x(:,ibeta,:,:),[],4))' / sqrt(nreps), ...
        plotOptions{:});
else
    % Plot alpha parameter as x-axis
    xdata = repmat(alphas, length(noiselevels), 1)';
    xlab = 'Alpha parameter';
    plotData = @(x, ibeta) errorbar(xdata, ...
        squeeze(nanmean(x(:,ibeta,:,:),4)), ...
        squeeze(nanstd(x(:,ibeta,:,:),[],4)) / sqrt(nreps), ...
        plotOptions{:});
end
xlims = [min(xdata(:)), max(xdata(:))];
xlims = xlims + [-0.1, 0.1];

% Set up axis array for linking plots
allax = zeros(1,4);

for ibeta = 1:length(betas)
    % Plot mean error
    ax1 = subplot(1, 4, 1);
    plotData(meanErrors, ibeta)
    hold on
    xlabel(xlab)
    ylabel('Mean displacement')
    title('Mean error')
    if ibeta == length(betas)
        legend(string(betas))
    end
    xlim(xlims)
    
    % Plot number of missing patterns
    ax2 = subplot(1, 4, 2);
    plotData(missingPatterns * 100 / (2 * gridSize(3)-1), ibeta);
    %plotData(timeStepsWithNoExtra / (gridSize(3)-1), ibeta);
    hold on
    xlabel(xlab)
    ylabel('% of missing patterns')
    title('Patterns missed')
    %title('% Time steps with no extras')
    xlim(xlims)
    ylim([0 100])
    
    % Plot number of excess patterns
    ax3 = subplot(1, 4, 3);
    plotData(extraPatterns/gridSize(3), ibeta);
    hold on
    xlabel(xlab)
    ylabel('# of extra patterns / time step')
    title('Excess patterns detected')
    xlim(xlims)
    
    % Plot mean number of iterations
    ax4 = subplot(1, 4, 4);
    plotData(meanIterations, ibeta);
    hold on
    xlabel(xlab)
    ylabel('Iteration steps')
    title('Convergence time')
    xlim(xlims)
    
    allax(ibeta,:) = [ax1, ax2, ax3, ax4];
end

% for iax = 1:4
%     linkaxes(allax(:,iax))
% end

%% Plot ROC curve
% figure
% truePositiveRate = 1 - missingPatterns / (2 * gridSize(3)-1);
% falsePositiveRate = 
% plotData = @(x, ibeta) errorbar(xdata, ...
%         squeeze(nanmean(x(:,ibeta,:,:),4)), ...
%         squeeze(nanstd(x(:,ibeta,:,:),[],4)) / sqrt(nreps), ...
%         plotOptions{:});

%% Print results
fprintf('Pattern A: %s\n', typeA)
distsA = pattDists(:,1);
fprintf('Mean critical point location error = %0.4e\n', nanmean(distsA))
fprintf('SEM = %0.4e, missed %i out of %i occurences\n\n', ...
    nanstd(distsA)/sum(~isnan(distsA)), sum(isnan(distsA)), gridSize(3)-1)
fprintf('Pattern B: %s\n', typeB)
distsB = pattDists(:,2);
fprintf('Mean critical point location error = %0.4e\n', nanmean(distsB))
fprintf('SEM = %0.4e, missed %i out of %i occurences\n\n', ...
    nanstd(distsB)/sum(~isnan(distsB)), sum(isnan(distsB)), gridSize(3)-1)
fprintf('%i extra critical points found (%0.3f per time step)\n', ...
    nextra, nextra/gridSize(3))
