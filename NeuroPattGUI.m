function varargout = NeuroPattGUI(varargin)
% NEUROPATTGUI MATLAB code for NeuroPattGUI.fig
%      NEUROPATTGUI, by itself, creates a new NEUROPATTGUI or raises the existing
%      singleton*.
%
%      H = NEUROPATTGUI returns the handle to a new NEUROPATTGUI or the handle to
%      the existing singleton*.
%
%      NEUROPATTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEUROPATTGUI.M with the given input arguments.
%
%      NEUROPATTGUI('Property','Value',...) creates a new NEUROPATTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NeuroPattGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NeuroPattGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NeuroPattGUI

% Rory Townsend, Aug 2018
% rory.townsend@sydney.edu.au

% Begin initialization code
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NeuroPattGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NeuroPattGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before NeuroPattGUI is made visible.
function NeuroPattGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NeuroPattGUI (see VARARGIN)

% Choose default command line output for NeuroPattGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Check if inputs exist and have the correct format
global data Fs params outputs
params = [];
inputErrorMsg = 'NeuroPatt GUI must be initialised with data and sampling frequency inputs i.e. NeuroPattGUI(DATA, FS), where DATA is a 3D or 4D numeric matrix and FS is a positive scalar.';
if nargin>3 && length(varargin)==2 && ...
        isnumeric(varargin{1}) && ...
        (ndims(varargin{1})==3 || ndims(varargin{1})==4) && ...
        isscalar(varargin{2}) && isnumeric(varargin{2}) && varargin{2}>0
    data = varargin{1};
    Fs = varargin{2};
else
    error(inputErrorMsg)
end

% Determine data class to determine size in bytes
switch(class(data))
    case 'double'
        bytesPerElem = 8;
    case 'single'
        bytesPerElem = 4;
    otherwise
        bytesPerElem = 0;
end

% Generate summary text based on input data
summaryStr = sprintf('%s\n%s', ...
    sprintf('Input data: %i rows, %i columns, %i time steps, %i repetitions,', ...
    size(data,1), size(data,2), size(data,3), size(data,4)), ...
    sprintf('%i Hz sampling frequency, %0.1f s, %i MB', ...
    Fs, size(data,3)/Fs, numel(data)*bytesPerElem/1e6));
set(handles.summaryText, 'String', summaryStr)

% Set outputs to be blank initially
outputs = [];

% UIWAIT makes NeuroPattGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NeuroPattGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global outputs
% Set output as output structure from analysis
varargout{1} = outputs;
% Get default command line output from handles structure
varargout{2} = handles.output;
delete(hObject);

% --- Executes on button press in zscoreCheckbox.
function zscoreCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to zscoreCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% If data is z-scored, the baseline must also be subtracted
if get(hObject, 'Value')==1 && get(handles.subtractBaselineCheckbox, 'Value')==0
    set(handles.subtractBaselineCheckbox, 'Value', 1);
end

% --- Executes on button press in spectraButton.
function spectraButton_Callback(hObject, eventdata, handles)
% hObject    handle to spectraButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data Fs params
updateParameters(handles);
plotSpectra(data(:,:,1:params.downsampleScale:end,:), Fs/params.downsampleScale);

% --- Executes on button press in chooseThresholdButton.
function chooseThresholdButton_Callback(hObject, eventdata, handles)
% hObject    handle to chooseThresholdButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data Fs params
updateParameters(handles);

% --- Executes on button press in continueButton.
function continueButton_Callback(hObject, eventdata, handles)
% hObject    handle to continueButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data Fs params outputs
updateParameters(handles);
% Set up temporary figure to output processing steps
% htemp = figure;
% htext = annotation('textbox',[0.1,0.1,0.8,0.8]);
% Perform all processing and save outputs
outputs = mainProcessingWithOutput(data, Fs, params);
outputs.params.isSurrogate = false;
% Open result summary figure
ResultsGUI(data, outputs);

% --- Function to update parameters based on text box inputs
function updateParameters(handles)
global Fs params
% Parse through all GUI input elements and match them to parameter names
inputHandleNames = {...
    'subtractBaselineCheckbox', 'zscoreCheckbox', 'downsampleScaleEdit', ...
    'filterCheckbox', 'hilbertButton', ...
    'centerFreqEdit', 'morletParamEdit', ...
    'minFreqEdit', 'maxFreqEdit', ...
    'ampButton', 'alphaEdit', 'betaEdit', ...
    'vfDecompCheckbox', 'complexDecompButton', 'nSVDmodesEdit', ...
    'minDurEdit', 'minRadiusEdit', 'pwThresholdEdit', ...
    'maxTimeGapEdit', 'minEdgeDistEdit', 'syncThresholdEdit', ...
    'maxDispEdit', 'combNodeFocusCheckbox', 'combSpiralCheckbox'};
paramNames = {...
    'subtractBaseline', 'zscoreChannels', 'downsampleScale', ...
    'filterData', 'useHilbert', ...
    'morletCfreq', 'morletParam', ...
    'hilbFreqLow', 'hilbFreqHigh', ...
    'useAmplitude', 'opAlpha', 'opBeta', ...
    'performSVD', 'useComplexSVD', 'nSVDmodes', ...
    'minDurationSecs', 'minCritRadius', 'planeWaveThreshold', ...
    'maxTimeGapSecs', 'minEdgeDistance', 'synchronyThreshold', ...
    'maxDisplacement', 'combineNodeFocus', 'combineStableUnstable'};
% Cycle through all parameters and update values
for iparam = 1:length(paramNames)
    thisHandle = handles.(inputHandleNames{iparam});
    if strcmp(get(thisHandle, 'Style'), 'edit')
        newVal = str2double(get(thisHandle, 'String'));
    else
        newVal = get(thisHandle, 'Value');
    end
    params = setNeuroPattParams(params, paramNames{iparam}, newVal, Fs);
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end
