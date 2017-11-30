% Function to set up time-varying spatiotemporal patterns to test optical
% flow calculations and pattern detection methods

clearvars
%close all

nx = 10;
ny = 12;
nt = 20;
[x, y, t] = meshgrid(1:nx, 1:ny, 1:nt);

% Optical flow parameters
alpha = 0.5;
beta = 1;

% Pattern detection parameters
params.maxTimeGap = 0;
params.minEdgeDist = 0;
params.minCritRadius = 0;

%% Phase plane wave
dir = [-1, -2]; % Moves up-right
w = 0.5;
k = 0.5;
planeWave = exp(1i * (w*t + k*(dir(1)*x + dir(2)*y)));

%% Phase sink pattern
loc = [2.5, 2.5];
vel = [0 0];
w = 0.3;
k = 1;
sink = exp(1i * (w*t + k*sqrt((x-loc(1)-vel(1)*t).^2 + ...
    (y-loc(2)-vel(2)*t).^2)));

%% Phase source pattern
loc = [2.5, 2.5];
vel = [0, 0];
w = 0.3;
k = 1;
source = exp(1i * (-w*t + k*sqrt((x-loc(1)-vel(1)*t).^2 + ...
    (y-loc(2)-vel(2)*t).^2)));

%% Phase spiral pattern
loc = [2.5, 2.5];
vel = [0, 0];
w = 0.3;
k=1;
spiralWave = exp(1i*(-w*t + ...
    angle(x-loc(1)-vel(1)*t + 1i*(y-loc(2)-vel(2)*t)) - ...
    k*sqrt((x-loc(1)-vel(1)*t).^2 + (y-loc(2)-vel(2)*t).^2)...
    ));

%% Phase saddle pattern
loc = [2.5, 2.5];
vel = [0, 0];
w = 0.3;
k = 0.1;
saddle = exp(1i * (-w*t + k*abs(x-loc(1)-vel(1)*t) - ...
    k*abs(y-loc(2)-vel(2)*t)));

%% Noise
noise = randn(size(spiralWave));

%% Make a combination of all other patterns
Aplane = 1;
Asink = 0;%linspace(1,0.1,nt);
Asource = 0;
Aspiral = 0;
Asaddle = 0;
Anoise = 0;

addAmp = @(A, wave) bsxfun(@times, permute(A, [3 1 2]), wave);
wave = addAmp(Aplane, planeWave) + addAmp(Asink, sink) + ...
    addAmp(Asource, source) + addAmp(Aspiral, spiralWave) + ...
    addAmp(Asaddle, saddle) + addAmp(Anoise, noise);

%% Compute optical flow for wave phase maps
[vx, vy, nit] = opticalFlow(wave, [], alpha, beta, 1);

%% Also compute optical flow for amplitude maps and view?
% figure
% [vx, vy, nit] = opticalFlow(real(wave), [], alpha, beta, 0);
% plotSnapshotsAndVfs(real(wave), vx, vy, 0, 1, 3)

%% Find patterns in velocity fields
[patterns, pattTypes, colNames, allPatternLocs, params]...
    = findAllPatterns(vx, vy, params);
allPatternLocs = cat(1, allPatternLocs{:});

%% View snapshots with detected patterns
figure
plotSnapshotsAndVfs(angle(wave), vx, vy, allPatternLocs, 1, 1, 3)
%title(sprintf('%0.1f average steps', mean(nit)))

figure
[U, S, V] = plotcsvd(vx+1i*vy, 4, 1:nt-1, 0);


