% Script to test accuracy of pattern detection for a range of simulated
% pattern types

close all

gridSize = [10 10 20];
alphas = linspace(0.1, 1.9, 10);
alphas = 0.1;
betas = [0.01 10];
noiselevels = logspace(-3, 1, 9);
noiselevels = 0.1;
nreps = 10;

%% Plane wave direction
% Set plane wave properties
amp_pw = 1;
freq_pw = 0.05;
wavelength_pw = 6.5;
dir_pw = [-1 -4]; % Moves up-right
% Find true optical flow direction
true_dir = angle(-dir_pw(1) - 1i*dir_pw(2));
% Generate plane wave
wave = generatePattern(gridSize, 'plane', amp_pw, freq_pw, ...
    wavelength_pw, [], dir_pw);

% Set up output arrays
errorMean = zeros(length(alphas), length(betas), length(noiselevels), nreps);
errorSEM = errorMean;
meanIterations = errorMean;

% Loop over all parameter values and noise levels
for ialpha = 1:length(alphas)
    fprintf('Processing waves with alpha = %s\n', string(alphas(ialpha)))
    for ibeta = 1:length(betas)
        for inoise = 1:length(noiselevels)
            for irep = 1: nreps
                % Add noise to wave
                thisNoise = randn(size(wave)) * amp_pw * noiselevels(inoise);
                thisWave = wave + thisNoise;
                
                % Test direction of vectors in velocity field
                [vx, vy, nit] = opticalFlow(thisWave, [], alphas(ialpha), ...
                    betas(ibeta), 1);
                vf_dir = angle(vx + 1i*vy);
                vf_error = abs(anglesubtract(vf_dir, true_dir));
                
                % Save mean errors, SEM, and steps to convergence
                errorMean(ialpha,ibeta,inoise, irep) = mean(vf_error(:));
                errorSEM(ialpha,ibeta,inoise, irep) = std(vf_error(:))/numel(vf_error);
                meanIterations(ialpha,ibeta,inoise, irep) = mean(nit);
            end
        end
    end
end

%% Plot results
figure
if length(noiselevels) > 1
    % Plot noise level as x-axis
    xdata = repmat(log10(noiselevels), length(alphas), 1)';
    xlab = 'Log noise level';
    extractData = @(x, ibeta) squeeze(x(:,ibeta,:,:))';
else
    % Plot alpha parameter as x-axis
    xdata = repmat(alphas, length(noiselevels), 1)';
    xlab = 'Alpha parameter';
    extractData = @(x, ibeta) squeeze(x(:,ibeta,:,:));
end

for ibeta = 1:length(betas)
    % Plot mean angular error
    subplot(2, length(betas), ibeta)
    errorbar(xdata, extractData(mean(errorMean,4), ibeta), ...
       extractData(std(errorMean,[],4), ibeta), '-o')
    xlabel(xlab)
    ylabel('Mean angular error')
    title(sprintf('Beta=%0.2e', betas(ibeta)))
    if ibeta == 1
        legend(string(alphas))
    end
    
    % Plot mean steps to convergence
    subplot(2, length(betas), length(betas) + ibeta)
    errorbar(xdata, extractData(mean(meanIterations,4), ibeta), ...
       extractData(std(meanIterations,[],4), ibeta), '-o')
    xlabel(xlab)
    ylabel('Mean steps to converge')
    title('Convergence time')
end

%fprintf('Plane wave mean angular error = %0.4e ± %0.4e (SEM)\n\n', ...
%    mean(vf_error(:)), std(vf_error(:))/numel(vf_error));
%histogram(vf_error(:))
