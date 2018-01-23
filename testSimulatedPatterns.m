% Script to test accuracy of pattern detection for a range of simulated
% pattern types

gridSize = [12 15 20];
alpha = 0.5;
beta = 1;

%% Plane wave direction
% Make plane wave
amp = 1;
freq = 0.01;
wavelength = 8;
vel = [-1 -4]; % Moves up-right
wave = generatePattern(gridSize, 'plane', amp, freq, wavelength, [], vel);
% Test direction of vectors in velocity field
[vx, vy, nit] = opticalFlow(wave, [], alpha, beta, 1);
vf_dir = angle(vx + 1i*vy);
true_dir = angle(-vel(1) - 1i*vel(2));
vf_error = abs(anglesubtract(vf_dir, true_dir));
fprintf('Plane wave mean angular error = %0.4e ± %0.4e (SEM)\n\n', ...
    mean(vf_error(:)), std(vf_error(:))/numel(vf_error))
%histogram(vf_error(:))

    
%% Spiral wave

% wave = generatePattern(gridSize, pattType, amp, freq, wavelength, loc, ...
%     vel, gausswidth);


%% Moving source and sink patterns

% Make source pattern
ampA = 1;
freq = 0.05;
wavelength = 5;
typeA = 'source';
fullNameA = 'unstableNode';
locA = [0.4,0.3] .* gridSize(1:2);
velA = [-0.1 0];
gausswidth = 1;
waveA = generatePattern(gridSize, typeA, ampA, freq, wavelength, ...
    locA, velA, gausswidth);

% Make sink pattern
ampB = 1.5;
typeB = 'sink';
fullNameB = 'stableNode';
locB = [0.8,0.6] .* gridSize(1:2);
velB = [0 0.05];
waveB = generatePattern(gridSize, typeB, ampB, freq, wavelength, locB, ...
    velB, gausswidth);

% Make velocity field
wave = waveA + waveB;
[vx, vy, nit] = opticalFlow(wave, [], alpha, beta, 1);

% Find all patterns
params.maxTimeGap = 0;
params.minCritRadius = 2;
params.minEdgeDist = params.minCritRadius;
params.minDuration = 0;
params.combineNodeFocus = true;
[patterns, pattTypes, colNames, allPatternLocs, params]...
    = findAllPatterns(vx, vy, params);
allPatternLocs = cat(1, allPatternLocs{:});

% Plot results
figure
plotSnapshotsAndVfs(angle(wave), vx, vy, allPatternLocs, 1, 1, 5)

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
            thisLoc = repmat(ilocs(itype,:), size(thisPatt,1), 1);
            % Calculate distance between detected and true
            minDist = min(sqrt(sum((thisPatt(:,[2,1])-thisLoc).^2, 2)));
            pattDists(itime, itype) = minDist;
        end
    end
end
extraPatts = size(allPatternLocs,1) - sum(~isnan(pattDists(:)));

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
    extraPatts, extraPatts/gridSize(3))

