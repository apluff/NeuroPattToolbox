% Script to test the accuracy of pattern detection and tracking with
% multiple simulated critical point patterns

%% Set parameters
gridSize = [12 12 20];
alphas = [0.1 0.5 1.5];
%alphas = linspace(0.1, 1.5, 8);
betas = [0.01 10];
%noiselevels = logspace(-3, 0, 7);
noiselevels = 0.03;
nreps = 1;

% Set parameters that apply to both patterns
freq = 0.1;
wavelength = 5;
gausswidth = 2;

% Set details of first pattern
ampA = 1;
typeA = 'source';
fullNameA = 'unstableNode';
locA = [0.4,0.3] .* gridSize(1:2);
velA = [-0.1 0];

% Set details of second pattern
ampB = 1.5;
typeB = 'sink';
fullNameB = 'stableNode';
locB = [0.8,0.6] .* gridSize(1:2);
velB = [0 0.05];

% Pattern detection parameters
params.maxTimeGap = 0;
params.minCritRadius = 2;
params.minEdgeDist = params.minCritRadius;
params.minDuration = 0;
params.combineNodeFocus = true;

%% Set up output variables
meanErrors = zeros(length(alphas), length(betas), ...
    length(noiselevels), nreps);
semErrors = meanErrors;
meanIterations = meanErrors;
missingPatterns = meanErrors;
extraPatterns = meanErrors;

%% Generate both patterns and combine
waveA = generatePattern(gridSize, typeA, ampA, freq, wavelength, ...
    locA, velA, gausswidth);
waveB = generatePattern(gridSize, typeB, ampB, freq, wavelength, locB, ...
    velB, gausswidth);
wave = waveA + waveB;

%% Loop over all parameter values and noise levels
tic
plotVfs = false;
for ialpha = 1:length(alphas)
    fprintf('Processing waves with alpha = %s\n', string(alphas(ialpha)))
    
    for ibeta = 1:length(betas)
        for inoise = 1:length(noiselevels)
            for irep = 1:nreps
                % Add noise proportional to total amplitude
                noise = randn(size(wave)) .* abs(wave) * ...
                    noiselevels(inoise);
                thisWave = wave + noise;
                
                % Make velocity field
                [vx, vy, nit] = opticalFlow(thisWave, [], ...
                    alphas(ialpha), betas(ibeta), 1);
                
                % Find all detected patterns
                [patterns, pattTypes, colNames, allPatternLocs, params] ...
                    = findAllPatterns(vx, vy, params);
                allPatternLocs = cat(1, allPatternLocs{:});
                
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
                
                % Test accuracy of critical point centres
                typeInds = zeros(1,2);
                typeInds(1) = find(strcmp(pattTypes, fullNameA));
                typeInds(2) = find(strcmp(pattTypes, fullNameB));
                pattDists = nan(gridSize(3)-1, 2);
                for itime = 1:gridSize(3)-1
                    % Detected critical points at this time step
                    ipatts = allPatternLocs(allPatternLocs(:,3)==itime, :);
                    % True positions at this time step
                    ilocA = locA + itime * velA;
                    ilocB = locB + itime * velB;
                    ilocs = [ilocA; ilocB];
                    % Check both critical point types
                    for itype = 1:2
                        thisPatt = ipatts(ipatts(:,4)==typeInds(itype), :);
                        if ~isempty(thisPatt)
                            thisLoc = repmat(ilocs(itype,:), ...
                                size(thisPatt,1), 1);
                            % Calculate distance between detected and true
                            minDist = min(sqrt(sum((thisPatt(:,[2,1]) - ...
                                thisLoc).^2, 2)));
                            pattDists(itime, itype) = minDist;
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
                nextra = size(allPatternLocs,1)- sum(~isnan(pattDists(:)));
                extraPatterns(ialpha,ibeta,inoise,irep) = nextra;
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

% Set up axis array for linking plots
allax = zeros(length(betas),4);

for ibeta = 1:length(betas)
    % Plot mean error
    ax1 = subplot(length(betas), 4, 4*(ibeta-1) + 1);
    plotData(meanErrors, ibeta)
    xlabel(xlab)
    ylabel('Mean displacement')
    title(sprintf('Mean error, Beta=%0.2e', betas(ibeta)))
    if ibeta == 1
        legend(string(alphas))
    end
    xlim([-3 0])
    
    % Plot number of missing patterns
    ax2 = subplot(length(betas), 4, 4*(ibeta-1) + 2);
    plotData(missingPatterns * 100 / (2 * gridSize(3)), ibeta);
    xlabel(xlab)
    ylabel('% of missing patterns')
    title('Patterns missed')
    xlim([-3 0])
    
    % Plot number of excess patterns
    ax3 = subplot(length(betas), 4, 4*(ibeta-1) + 3);
    plotData(extraPatterns/gridSize(3), ibeta);
    xlabel(xlab)
    ylabel('# of extra patterns / time step')
    title('Excess patterns detected')
    xlim([-3 0])
    
    % Plot mean number of iterations
    ax4 = subplot(length(betas), 4, 4*(ibeta-1) + 4);
    plotData(meanIterations, ibeta);
    xlabel(xlab)
    ylabel('Iteration steps')
    title('Convergence time')
    xlim([-3 0])
    
    allax(ibeta,:) = [ax1, ax2, ax3, ax4];
end

for iax = 1:4
    linkaxes(allax(:,iax))
end

%% Print results
fprintf('Pattern A: %s\n', typeA)
distsA = pattDists(:,1);
fprintf('Mean critical point location error = %0.4e\n', nanmean(distsA))
fprintf('SEM = %0.4e, missed %i out of %i occurences\n\n', ...
    nanstd(distsA)/sum(~isnan(distsA)), sum(isnan(distsA)), gridSize(3))
fprintf('Pattern B: %s\n', typeB)
distsB = pattDists(:,2);
fprintf('Mean critical point location error = %0.4e\n', nanmean(distsB))
fprintf('SEM = %0.4e, missed %i out of %i occurences\n\n', ...
    nanstd(distsB)/sum(~isnan(distsB)), sum(isnan(distsB)), gridSize(3))
fprintf('%i extra critical points found (%0.3f per time step)\n', ...
    nextra, nextra/gridSize(3))
