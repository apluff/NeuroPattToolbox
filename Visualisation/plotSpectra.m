function hfig = plotSpectra(signal,fs,fmin,fmax,dim)
% PLOTSPECTRA plots the Fourier spectrum and and wavelet spectrogram of
% input signal SIGNAL with sampling frequency FS. Spectra are computed
% between minimum frequency FMIN (default 1 Hz) and maximum frequency FMAX
% (default 100 Hz). Signal is assumed to be in the from XxYxTIME, otherwise
% optional input DIM specifies the dimension of time.

% Process inputs
if nargin < 3 || ~isscalar(fmin) || ~isnumeric(fmin) || fmin<=0
    fmin = 1;
end
if nargin < 4 || ~isscalar(fmax) || ~isnumeric(fmax) || fmax<=fmin
    fmax = 50;
end
if nargin < 5 && ndims(signal)>=3
    dim = 3;
end

% Reshape data so that required dim is first
dimOrder = [dim, 1:dim-1, dim+1:ndims(signal)];
signal = permute(signal, dimOrder);
signal = signal(:,:);

% Plot Fourier spectrum
fftData = fftshift(fft(signal,[],1), 1);
sigLength = size(signal,1);
fftFreqs = fs*(0:(sigLength/2))/sigLength;
fftData = fftData(ceil(sigLength/2):end, :);
hfig = figure;
subplot(1,5,1:2)
plot(fftFreqs, mean(abs(fftData),2))
xlim([fmin, fmax])
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Fourier power spectrum')

% Plot morlet spectrogram
% TODO: Add these parameters as inputs instead of hard-coding them here
morletParam = 7;
nsteps = 50;
cfreqs = linspace(fmin, fmax, nsteps);
cfx = morletWaveletTransform(signal, fs, cfreqs, morletParam, 1);
% Generate time-frequency power plot
subplot(1,5,3:5)
time = (1:sigLength)/fs;
imagesc(time, cfreqs, mean(abs(cfx),3))
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Morlet wavelet power')
axis xy

end

