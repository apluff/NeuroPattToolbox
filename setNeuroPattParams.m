function params = setNeuroPattParams(params, name, value, fs)
% SETNEUROPATTPARAMS Sets all parameters in the NeuroPatt toolbox
%   PARAMS = SETPATTERNPARAMS(FS) sets all parameters in the PARAMS
%   structure using parameter values listed below using sampling frequency
%   FS (in Hz)
%
%   PARAMS = SETPATTERNPARAMS(PARAMS, NAME, VALUE, FS) updates existing
%   parameter structure PARAMS by changing the parameter called NAME to be
%   set to VALUE, using sampling frequency FS. If PARAMS input is empty,
%   parameters apart from NAME will be filled by the default values listed
%   below.
%

%% Parse inputs
if nargin==1 && isnumeric(params)
    % If only sampling frequency is defined, set all parameters based on
    % values below
    fs = params;
    setAllParams = true;
    name = '';
    params = [];
elseif nargin==4 && (isempty(params) || isstruct(params)) && ...
        ischar(name) && isnumeric(fs)
    if isempty(params)
        setAllParams = true;
    else
        setAllParams = false;
    end
else
    error('Invalid input format! Inputs must be setPatternParams(FS) or setPatternParams(PARAMS,NAME,VALUE,FS).')
end

if setAllParams
%% Preprocessing parameters
% Temporal downsampling scale after filtering for faster optical flow and
% pattern detection calculations (default 1, corresponding to no
% downsampling)
downsampleScale = 5;
% Flag to de-mean all channels before processing (default true)
subtractBaseline = true;
% Flag to z-score all channels before processing (default false)
zscoreChannels = false;
    
%% Filtering parameters
% Flag to filter data, otherwise raw unfiltered data will be used (default
% true)
filterData = true;
% Flag to use the Hilbert transform for filtering, otherwise Morlet
% wavelets will be used (default false)
useHilbert = false;
% Morlet transform parameters
% Centre frequency (in Hz)
morletCfreq = 6;
% Morlet wavelet scale parameter: higher values have better frequency
% resolution but lower temporal resolution (default 5)
morletParam = 5;
% Hilbert transform parameters
% Frequency limits of band-pass filter (in Hz)
hilbFreqLow = 1;
hilbFreqHigh = 4;

%% Optical flow parameters
% Smoothing parameter: higher values will give a smoother velocity field
% (typically 0<OPALPHA<5, default 0.5).
opAlpha = 0.5;
% Non-linearity penalty parameter: close to zero will be highly non-linear,
% large values will be approximately linear, resulting in faster
% computations but possibly less accurate flow fields (default 1)
opBeta = 1;
% Flag to calculate amplitude velocity fields rather than phase (default
% false)
useAmplitude = false;

%% Singular value decomposition parameters
% Flag to perform SVD on velocity fields to extract spatiotemporal modes,
% otherwise this will be skipped (default true)
performSVD = true;
% Flag to perform a complex decomposition, allowing spatial modes that can
% rotate over time. Otherwise, real decomposition will be used (default
% false)
useComplexSVD = false;
% Number of modes to plot (default 6)
nSVDmodes = 6;

%% Pattern parameters
% Minimum order parameter for a plane wave to be detected (default 0.85)
planeWaveThreshold = 0.85;
% Maximum velocity field magnitude for synchrony to be detected (default
% 0.85)
synchronyThreshold = 0.85;
% Minimum duration of a pattern for it to be stored (in seconds, default 0)
minDurationSecs = 0.02;
% Maximum duration between critical points (or synchrony/plane waves) for
% them to be counted as the same pattern (in seconds, default 0)
maxTimeGapSecs = 0.005;
% Maxiumum displacement between critical points between time steps for them
% to be counted as the same pattern (measured in grid spaces, default 1)
maxDisplacement = 1;
% Minimum spatial radius for a critical point to occupy for it to be
% counted, quantified by the winding number (measured in grid spaces,
% default 2)
minCritRadius = 2;
% Minimum distance from the edge of the system (in grid spaces, default 2)
minEdgeDistance = 2;
% Boolean paramter to combine node and focus type critical points (default
% false)
combineNodeFocus = false;
% Boolean parameter to combine stable and unstable critical points (default
% false)
combineStableUnstable = false;

else
    % Add in a temporary value for this minimum frequency so that the
    % maximum frequency can be compared against it
    hilbFreqLow = params.hilbFreqLow;

end

%% Verify parameters and add to output structure

% Generate cell array with all parameter names and limitations
% Format: Name, Flag for logical, Minimum value, Maximum value
allParams = {...
    'downsampleScale', false, 1-eps, []; ... % Pre-processing parameters
    'subtractBaseline', true, [], []; ...
    'zscoreChannels', true, [], []; ...
    'filterData', true, [], []; ... % Filtering parameters
    'useHilbert', true, [], []; ...
    'morletCfreq', false, 0, []; ...
    'morletParam', false, 0, []; ...
    'hilbFreqLow', false, 0, []; ...
    'hilbFreqHigh', false, hilbFreqLow-eps, []; ...
    'opAlpha', false, 0, []; ... % Optical flow parameters
    'opBeta', false, 0, []; ...
    'useAmplitude', true, [], []; ...
    'performSVD', true, [], []; ... % SVD parameters
    'useComplexSVD', true, [], []; ...
    'nSVDmodes', false, 0, 20; ...
    'planeWaveThreshold', false, 0-eps, 1+eps; ... % Pattern detection parameters
    'synchronyThreshold', false, 0-eps, 1+eps; ...
    'minDurationSecs', false, 0-eps, []; ...
    'maxTimeGapSecs', false, 0-eps, []; ...
    'maxDisplacement', false, 0-eps, []; ...
    'minCritRadius', false, 0-eps, []; ...
    'minEdgeDistance', false, 0-eps, []; ...
    'combineNodeFocus', true, [], []; ...
    'combineStableUnstable', true, [], [] ...
    };

% If only changing one parameter, find it in the list
if ~setAllParams
    allParams = allParams(cellfun(@(x) strcmp(x,name), allParams(:,1)), :);
end

% Loop over all parameters in the list
for iparam = 1:size(allParams, 1)
    [pname, logicalFlag, minVal, maxVal] = allParams{iparam, :};
    % Choose default or user input value
    if strcmp(pname, name)
        ivalue = value;
    else
        ivalue = eval(pname);
    end
    % Test if parameter value is valid and add to output structure
    if logicalFlag == 1
        params.(pname) = testLogical(ivalue, pname);
    else
        params.(pname) = testNumeric(ivalue, pname, minVal, maxVal);
    end
end

% Set secondary parameters that are defined by the sampling frequency
params.maxTimeGap = floor(params.maxTimeGapSecs * fs / params.downsampleScale);
params.minDuration = max(1, round(params.minDurationSecs * fs / params.downsampleScale));

end

function value = testNumeric(value, varName, minval, maxval)
    % Subfunction to test numeric values for variable VARNAME with an 
    % optional minimum of MINVAL and maximum of MAXVAL, throwing an error
    % if invalid
    if ~isnumeric(value) || ~isscalar(value) || ...
            (nargin>=3 && ~isempty(minval) && value<=minval) || ...
            (nargin>=4 && ~isempty(maxval) && value>=maxval)
        
        % Construct error message to show the limitations on the parameter
        if nargin>=4 && ~isempty(minval) && ~isempty(maxval)
            errorMsgEnd = sprintf(' and must be greater than %0.1f and less than %0.1f.', ...
                minval, maxval);
        elseif nargin>=4 && ~isempty(maxval)
            errorMsgEnd = sprintf(' and must be less than %0.1f.', maxval);
        elseif nargin>=3 && ~isempty(minval)
            errorMsgEnd = sprintf(' and must be greater than %0.1f.', minval);
        else
            errorMsgEnd = '.';
        end
                
        error('Invalid value for parameter %s. Value must be a numeric scalar%s', ...
            varName, errorMsgEnd)
    end
end

function value = testLogical(value, varName)
    % Subfunction to test logical values for variable VARNAME, throwing an
    % error if invalid
    if isscalar(value) && value
        value = true;
    elseif isscalar(value) && ~value
        value = false;
    else
        error('Invalid value for parameter %s. Value must resolve to a logical expression.', varName)
    end
end


