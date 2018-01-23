% Function to set up time-varying spatiotemporal patterns to test optical
% flow calculations and pattern detection methods

clearvars
close all

nx = 20;
ny = 22;
nt = 20;
[x, y, t] = meshgrid(1:nx, 1:ny, 1:nt);

% Optical flow parameters
alpha = 0.5;
beta = 1;

% Pattern detection parameters
params.maxTimeGap = 0;
params.minEdgeDist = 0;
params.minCritRadius = 0;

% Choose amplitude of patterns
Aplane = 0;
Asink = 1;%linspace(1,0.1,nt);
Asource = 0;
Aspiral = 0;
Asaddle = 0;
Anoise = 0;

% Function to generate Gaussian spatial profile
gaussian = @(c, loc) exp(-1/(2*c^2) * (((x(:,:,1)-loc(1)).^2 + ...
    (y(:,:,1)-loc(2)).^2)));

%% Phase plane wave
dir = [-1, -2]; % Moves up-right
w = 0.5; % Angular frequency
k = 0.5; % Spatial wavenumber
planeWave = exp(1i * (w*t + k*(dir(1)*x + dir(2)*y)));
wave = addAmplitudeProfiles(planeWave, [], Aplane);

%% Phase sink pattern
loc = [4, 3];
vel = [0, 0]; % Velocity (per time step) of critical point
w = 0.3; % Angular frequency
k = 1; % Spatial wavenumber
c = 2; % Width parameter of Gaussian mask
sink = exp(1i * (w*t + k*sqrt((x-loc(1)-vel(1)*t).^2 + ...
    (y-loc(2)-vel(2)*t).^2)));
wave = wave + addAmplitudeProfiles(sink, gaussian(c,loc), Asink);

%% Phase source pattern
loc = [7, 6];
vel = [0, 0]; % Velocity (per time step) of critical point
w = 0.8; % Angular frequency
k = 1; % Spatial wavenumber
c = 2; % Width parameter of Gaussian mask
source = exp(1i * (-w*t + k*sqrt((x-loc(1)-vel(1)*t).^2 + ...
    (y-loc(2)-vel(2)*t).^2)));
wave = wave + addAmplitudeProfiles(source, gaussian(c,loc), Asource);

%% Phase spiral pattern
loc = [2.5, 2.5];
vel = [0, 0]; % Velocity (per time step) of critical point
w = 0.3; % Angular frequency
k = 1; % Spatial wavenumber
c = 2; % Width parameter of Gaussian mask
spiralWave = exp(1i*(-w*t + ...
    angle(x-loc(1)-vel(1)*t + 1i*(y-loc(2)-vel(2)*t)) - ...
    k*sqrt((x-loc(1)-vel(1)*t).^2 + (y-loc(2)-vel(2)*t).^2)...
    ));
wave = wave + addAmplitudeProfiles(spiralWave, gaussian(c,loc), Aspiral);

%% Phase saddle pattern
loc = [2.5, 2.5];
vel = [0, 0]; % Velocity (per time step) of critical point
w = 0.3; % Angular frequency
k = 0.1; % Spatial wavenumber
c = 2; % Width parameter of Gaussian mask
saddle = exp(1i * (-w*t + k*abs(x-loc(1)-vel(1)*t) - ...
    k*abs(y-loc(2)-vel(2)*t)));
wave = wave + addAmplitudeProfiles(saddle, gaussian(c,loc), Asaddle);

%% Noise
noise = randn(size(wave));
wave = wave + Anoise * noise;

%% Compute optical flow for wave phase maps
[vx, vy, nit] = opticalFlow(wave, [], alpha, beta, 1);

%% Also compute optical flow for amplitude maps and view?
figure
[vxa, vya, nita] = opticalFlow(real(wave), [], alpha, beta, 0);
plotSnapshotsAndVfs(real(wave), vxa, vya, [], 0, 1, 3)

%% Find patterns in velocity fields
[patterns, pattTypes, colNames, allPatternLocs, params]...
    = findAllPatterns(vx, vy, params);
allPatternLocs = cat(1, allPatternLocs{:});

%% View snapshots with detected patterns
figure
plotSnapshotsAndVfs(angle(wave), vx, vy, allPatternLocs, 1, 1, 3)
title(sprintf('%0.1f average steps', mean(nit)))

figure
[U, S, V] = plotcsvd(vx+1i*vy, 4, 1:nt-1, 0);

