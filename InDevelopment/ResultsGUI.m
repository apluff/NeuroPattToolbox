function varargout = ResultsGUI(varargin)
% RESULTSGUI MATLAB code for ResultsGUI.fig
%      RESULTSGUI, by itself, creates a new RESULTSGUI or raises the existing
%      singleton*.
%
%      H = RESULTSGUI returns the handle to a new RESULTSGUI or the handle to
%      the existing singleton*.
%
%      RESULTSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESULTSGUI.M with the given input arguments.
%
%      RESULTSGUI('Property','Value',...) creates a new RESULTSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ResultsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ResultsGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ResultsGUI

% Rory Townsend, Aug 2018
% rory.townsend@sydney.edu.au

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ResultsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ResultsGUI_OutputFcn, ...
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


% --- Executes just before ResultsGUI is made visible.
function ResultsGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ResultsGUI (see VARARGIN)

% Choose default command line output for ResultsGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Set global variables for other functions
global data results

% Populate header summary text
data = varargin{1};
results = varargin{2};
dataLength = results.nTimeSteps;
dataLengthSecs = dataLength/results.Fs;
overviewStr = sprintf('Processed %0.1f s recording with %i repetitions sampled at %i Hz in %0.1f s, finished at %s.', ...
    dataLengthSecs, length(results.patterns), results.Fs, ...
    seconds(results.processTime), datestr(datetime));
set(handles.overviewText, 'String', overviewStr)

% Populate pattern summary statistics table
summaryTable = zeros(length(results.patternTypes), 4, ...
    length(results.patterns));
% Calculate pattern statistics in each trial repetition
for icell = 1:length(results.patterns)
    thisCell = results.patterns{icell};
    % Loop over each pattern type
    patternArray = makeActivePatternsArray(thisCell, length(results.patternTypes), dataLength);
    for iptype = 1:length(results.patternTypes)
        thisPattern = thisCell(thisCell(:,1)==iptype, :);
        % Number of this pattern detected
        summaryTable(iptype, 1, icell) = size(thisPattern,1);
        % Number detected per second
        summaryTable(iptype, 2, icell) = size(thisPattern,1)/dataLengthSecs;
        % Fraction of total time where at least one of this pattern is
        % active
        summaryTable(iptype, 3, icell) = sum(patternArray(iptype,:)>0) ...
            / dataLength * 100;
        % Mean duration of all of this pattern's occurences
        summaryTable(iptype, 4, icell) = mean(thisPattern(:,4))/results.Fs;
    end
end
set(handles.pattSummaryTable, 'RowName', results.patternTypes, ...
    'Data', num2cell(nanmean(summaryTable,3)))

% Populate pattern transition table
set(handles.pattTransitionsTable, 'RowName', results.patternTypes, ...
    'ColumnName', results.patternTypes, ...
    'Data', nanmean(results.pattTransitionsObs - ...
    results.pattTransitionsExp, 3)) % TO-DO Divided by occurrences?

% Hide surrogate data options if data already is surrogate
if results.params.isSurrogate
    set(handles.headingText, 'String', 'Surrogate Data Results')
    set(handles.surrogateHeadingText, 'Visible', 'off');
    set(handles.surrogateRepsEdit, 'Visible', 'off');
    set(handles.surrogateTrialsText, 'Visible', 'off');
    set(handles.surrogateAnalysisButt, 'Visible', 'off');
end

% UIWAIT makes ResultsGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Executes when user attempts to close figure1.
function varargout = figure1_CloseRequestFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
end
% The GUI is no longer waiting, just close it
delete(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = ResultsGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
%delete(hObject)


% --- Executes on button press in durationDistributionButt.
function durationDistributionButt_Callback(hObject, eventdata, handles)
% hObject    handle to durationDistributionButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global results

% Plot distributions of pattern durations and displacements
figure
colNames = results.patternResultColumns;
Fs = results.Fs;
allPatts = cat(1, results.patterns{:});
unqPatts = unique(allPatts(:,1));
maxDur = max(allPatts(:, strcmp(colNames, 'duration')))/Fs;
maxDisp = max(allPatts(:, strcmp(colNames, 'meanDisplacement')))*Fs;
for ipatt = 1:length(unqPatts)
    thisPatt = allPatts(allPatts(:,1)==unqPatts(ipatt), :);
    % Histogram of durations
    subplot(2, length(unqPatts), ipatt)
    thisDurs = thisPatt(:, strcmp(colNames, 'duration'));
    histogram(thisDurs/Fs, linspace(0, maxDur, 10))
    xlabel('Duration (s)')
    ylabel('Counts')
    title(sprintf('%s, mean %0.3g', ...
        results.patternTypes{unqPatts(ipatt)}, nanmean(thisDurs/Fs)))
    
    % Histogram of displacements
    subplot(2, length(unqPatts), length(unqPatts)+ipatt)
    thisDisp = thisPatt(:, strcmp(colNames, 'meanDisplacement'));
    histogram(thisDisp*Fs, linspace(0, maxDisp, 10))
    xlabel('Displacement (grid/s)')
    ylabel('Counts')
    title(sprintf('Mean %0.3g', nanmean(thisDisp*Fs)))
end


% --- Executes on button press in locationDistributionButt.
function locationDistributionButt_Callback(hObject, eventdata, handles)
% hObject    handle to locationDistributionButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global results
figure
% Plot plane wave direction and pattern center locations over space
plotPatternLocs(results.patternLocs, results.patternTypes, ...
    results.nTimeSteps/results.Fs, length(results.patterns));

% --- Executes on button press in surrogateAnalysisButt.
function surrogateAnalysisButt_Callback(hObject, eventdata, handles)
% hObject    handle to surrogateAnalysisButt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data results
% Generate surrogate data based on input
nreps = round(str2num(get(handles.surrogateRepsEdit, 'String')));
szdata = size(data);
surData = randn([szdata(1:3), nreps]);
surData = bsxfun(@plus, surData, nanmean(data(:,:,:), 3));
surData = bsxfun(@times, surData, nanstd(data(:,:,:), [], 3));

% Repeat analysis with surrogate data
disp('_________________________________________________________')
disp('Repeating all analysis with surrogate data (white noise).')
surResults = mainProcessingWithOutput(surData, results.Fs, results.params,[],1,1);
surResults.params.isSurrogate = true;
results.surReps = nreps;
results.surPatterns = surResults.patterns;
results.surPatternLocs = surResults.patternLocs;
results.surProcessTime = surResults.processTime;
ResultsGUI(data, surResults);


% --- Executes on button press in closeButton.
function closeButton_Callback(hObject, eventdata, handles)
% hObject    handle to closeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all