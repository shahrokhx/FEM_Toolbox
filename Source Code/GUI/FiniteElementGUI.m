function varargout = FiniteElementGUI(varargin)
% FiniteElementGUI MATLAB code for FiniteElementGUI.fig
%      FiniteElementGUI, by itself, creates a new FiniteElementGUI or raises the existing
%      singleton*.
%
%      H = FiniteElementGUI returns the handle to a new FiniteElementGUI or the handle to
%      the existing singleton*.
%
%      FiniteElementGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FiniteElementGUI.M with the given input arguments.
%
%      FiniteElementGUI('Property','Value',...) creates a new FiniteElementGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FiniteElementGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FiniteElementGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FiniteElementGUI

% Last Modified by GUIDE v2.5 17-Sep-2021 05:48:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FiniteElementGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @FiniteElementGUI_OutputFcn, ...
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


% --- Executes just before FiniteElementGUI is made visible.
function FiniteElementGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FiniteElementGUI (see VARARGIN)

% Choose default command line output for FiniteElementGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FiniteElementGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% Main Initialization
initializeGUI(hObject, eventdata, handles)

function varargout = FiniteElementGUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function edit1_Callback(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton7_Callback(hObject, eventdata, handles)
function pushbutton8_Callback(hObject, eventdata, handles)

function edit2_Callback(hObject, eventdata, handles)
function edit2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% +-----------------------------------------------------------------------+
% |                    Input File Section                                 |
% +-----------------------------------------------------------------------+
function edit3_Callback(hObject, eventdata, handles)
function edit3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Browse
function pushbutton9_Callback(hObject, eventdata, handles)
[file,path] = uigetfile({'*.txt';'*.inp';'*.*'},'Choose an Input File...');
if file ~=0
    fileName = [path,file];
    set(handles.edit3,'String',fileName);

    try
        Model = inpFileReader(fileName);
        graphicalView(hObject, eventdata, handles, Model);
        printSummary(Model)
        
        if get(handles.checkbox7,'Value')
            save('Model','Model')
        end
    catch
        errordlg('Invalid Input File!','Error');
    end
    

end

% View/Edit input file
function pushbutton16_Callback(hObject, eventdata, handles)
fileName = get(handles.edit3,'String');
try
    edit (fileName);
catch
    errordlg('Check the Input File','Error');
end

% Re-load
function pushbutton17_Callback(hObject, eventdata, handles)
fileName = get(handles.edit3,'String');
% Model = inpFileReader(fileName);
% printSummary(Model)
try
    Model = inpFileReader(fileName);
    graphicalView(hObject, eventdata, handles, Model);
    printSummary(Model)
    
    
    if get(handles.checkbox7,'Value')
        save('Model','Model')
    end
catch
    errordlg('Check the Input File','Error');
end



function popupmenu1_Callback(hObject, eventdata, handles)
function popupmenu1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%View/Edit Solver
function pushbutton10_Callback(hObject, eventdata, handles)
solverFunc = [];
path(pathdef);
pathThis = mfilename('fullpath');
pathThis(end-length(mfilename):end)=[];
addpath(fullfile(pathThis));

switch get(handles.popupmenu1,'Value')
    case 1
        container  = 'PlaneTruss'; 
        solverFunc = 'trussAnalysisExplicit';
    case 2
        container  = 'PlaneTruss'; 
        solverFunc = 'trussAnalysisNumerical';
    case 3
        container  = 'PlaneFrame'; 
        solverFunc = 'planeFrameExplicit';
    case 4
        container  = 'PlaneFrame'; 
        solverFunc = 'planeFrameNumerical';
    case 5
        container  = 'PlaneStress PlaneStrain'; 
        solverFunc = 'femGeneral';
    case 6
        container  = 'PlaneStress PlaneStrain'; 
        solverFunc = 'femGeneral';
end
if ~isempty(solverFunc)
    pathThis = mfilename('fullpath');
    pathThis(end-length(mfilename):end)=[];
    addpath(fullfile(pathThis,'..',container,'lib'))
    
    edit(solverFunc)
end

function checkbox7_Callback(hObject, eventdata, handles)
function checkbox8_Callback(hObject, eventdata, handles)

% Run Analysis
function pushbutton11_Callback(hObject, eventdata, handles)
path(pathdef);
pathThis = mfilename('fullpath');
pathThis(end-length(mfilename):end)=[];
addpath(fullfile(pathThis));

switch get(handles.popupmenu1,'Value')
    case 1
        solver = @(M)trussAnalysisExplicit(M);
        container  = 'PlaneTruss'; 
        solverFunc = 'trussAnalysisExplicit';
    case 2 
        solver = @(M)trussAnalysisNumerical(M);
        container  = 'PlaneTruss'; 
        solverFunc = 'trussAnalysisNumerical';
    case 3
        solver = @(M)planeFrameExplicit(M);
        container  = 'PlaneFrame'; 
        solverFunc = 'planeFrameExplicit';
    case 4
        solver = @(M)planeFrameNumerical(M);
        container  = 'PlaneFrame'; 
        solverFunc = 'planeFrameNumerical';
    case 5
        solver = @(M)femGeneral(M);
        container  = 'PlaneStress PlaneStrain'; 
        solverFunc = 'femGeneral';
    case 6
        solver = @(M)femGeneral(M);
        container  = 'PlaneStress PlaneStrain'; 
        solverFunc = 'femGeneral';
end

pathThis = mfilename('fullpath');
pathThis(end-length(mfilename):end)=[];
addpath(fullfile(pathThis,'..',container,'lib'))
addpath(fullfile(pathThis,'..',container,'lib_io'))

fileName = get(handles.edit3,'String');
try
    Model = inpFileReader(fileName);
    if get(handles.checkbox7,'Value')
        save('Model','Model')
    end
catch
    errordlg('Check the Input File','Error');
end

%--------------------------
Model.analysis.saveToFile       =    get(handles.checkbox15,'Value');
Model.analysis.outputFileName   =    get(handles.edit6,'String');
Model.analysis.showResults      =    1;
Model.analysis.showDetails      =    get(handles.checkbox8,'Value');
%--------------------------


try
    Solution = solver(Model);
    if get(handles.checkbox7,'Value')
        save('Solution','Solution')
    end
    
    % printResult(Solution)   %%% NOT HERE
    
%     if get(handles.checkbox8,'Value')
%         % printing stiffness matrices
%     end
%     if get(handles.checkbox15,'Value')
%         % save output to file
%         outFileName = get(handles.edit6,'String');
%         if ~isempty(outFileName)
%             try 
%                 outFile = fopen(outFileName, 'w'); 
%                 printSummary(Model,outFile)
%                 printResult(Solution,outFile)
%                 fclose(outFile);
%             catch
%                 errordlg('Check the Output File Name','Error');
%             end
%         end
%     end
    u = Solution.nodalDisplacements;
    coords = Model.geometry.coordinates;
    u = u(:,1:size(coords,2));
    
    maxU = max(max(abs(u)));
    maxL = max(max(abs(coords))) - min(min(coords));
    
    FACTOR_MAGNIFYING = 0.10;  %10 percent of maxL
    exaggFactor = FACTOR_MAGNIFYING * maxL / maxU;
    
    ModelDef = Model;
    ModelDef.geometry.coordinates = coords + exaggFactor.*u;
    
    graphicalView(hObject, eventdata, handles, Model)
    graphicalViewDef(hObject, eventdata, handles, ModelDef)
catch
    errordlg({'Check the Chosen Solver:'; ...
              [solverFunc,'.m']},...
              'Error');
end


% save to file
function checkbox15_Callback(hObject, eventdata, handles)
if get(handles.checkbox15,'Value')
    set (handles.edit6,'Visible','on')
else
    set (handles.edit6,'Visible','off')
end


function edit6_Callback(hObject, eventdata, handles)
function edit6_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% +-----------------------------------------------------------------------+
% |    graphic options           checkbox                                 |
% +-----------------------------------------------------------------------+
% Node Label Checkbox
function checkbox9_Callback(hObject, eventdata, handles)

if get(handles.checkbox9,'Value')
    set (handles.uipanel8,'Visible','on')
else
    set (handles.uipanel8,'Visible','off')
end
liveGraphics(hObject, eventdata, handles);

% Element Label Checkbox
function checkbox10_Callback(hObject, eventdata, handles)
if get(handles.checkbox10,'Value')
    set (handles.uipanel9,'Visible','on')
else
    set (handles.uipanel9,'Visible','off')
end
liveGraphics(hObject, eventdata, handles);

% Node Color
function checkbox11_Callback(hObject, eventdata, handles)
if get(handles.checkbox11,'Value')
    set (handles.pushbutton12,'Visible','on')
else
    set (handles.pushbutton12,'Visible','off')
end
liveGraphics(hObject, eventdata, handles);

%Element Color
function checkbox12_Callback(hObject, eventdata, handles)
if get(handles.checkbox12,'Value')
    set (handles.pushbutton13,'Visible','on')
else
    set (handles.pushbutton13,'Visible','off')
end
liveGraphics(hObject, eventdata, handles);

function checkbox13_Callback(hObject, eventdata, handles)
liveGraphics(hObject, eventdata, handles);

function checkbox14_Callback(hObject, eventdata, handles)
liveGraphics(hObject, eventdata, handles);

% +-----------------------------------------------------------------------+
% |    graphic options      color puchbutton                              |
% +-----------------------------------------------------------------------+
function pushbutton12_Callback(hObject, eventdata, handles)
c = uisetcolor('Select a Color');
set(handles.pushbutton12,'BackgroundColor',c); 
liveGraphics(hObject, eventdata, handles);
function pushbutton13_Callback(hObject, eventdata, handles)
c = uisetcolor('Select a Color');
set(handles.pushbutton13,'BackgroundColor',c); 
liveGraphics(hObject, eventdata, handles);
function pushbutton14_Callback(hObject, eventdata, handles)
c = uisetcolor('Select a Color');
set(handles.pushbutton14,'BackgroundColor',c); 
liveGraphics(hObject, eventdata, handles);
function pushbutton15_Callback(hObject, eventdata, handles)
c = uisetcolor('Select a Color');
set(handles.pushbutton15,'BackgroundColor',c); 
liveGraphics(hObject, eventdata, handles);

% +-----------------------------------------------------------------------+
% |    graphic options      View puchbutton                               |
% +-----------------------------------------------------------------------+
function pushbutton1_Callback(hObject, eventdata, handles)
view(handles.axes1,2) %XY
function pushbutton2_Callback(hObject, eventdata, handles)
view(handles.axes1,[0,0])
function pushbutton3_Callback(hObject, eventdata, handles)
view(handles.axes1,[90,0])
function pushbutton4_Callback(hObject, eventdata, handles)
view(handles.axes1,3) %3D view



function edit4_Callback(hObject, eventdata, handles)
function edit4_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit5_Callback(hObject, eventdata, handles)
function edit5_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% +-----------------------------------------------------------------------+
% |               A U X I L I A R Y    F U N C T I O N S                  |
% |    All functions that are not directly related to the GUI             |
% |    should be added here                                               |
% +-----------------------------------------------------------------------+
function initializeGUI(hObject, eventdata, handles)
cla(handles.axes1)
% initialize the graphical view interface
set(handles.checkbox9 ,'Value',0)
set(handles.checkbox10,'Value',0)
set(handles.checkbox11,'Value',1)
set(handles.checkbox12,'Value',1)
set(handles.checkbox13,'Value',1)  % show bo
set(handles.checkbox14,'Value',1)  % show loads
set(handles.pushbutton12,'BackgroundColor',[1   0   0]); % Node/Edge color
set(handles.pushbutton13,'BackgroundColor',[0   1   1]); % Element color
set(handles.pushbutton15,'BackgroundColor',[1   1 0.2]); % node label
set(handles.pushbutton14,'BackgroundColor',[1   0   1]); % element label
set(handles.uipanel8 , 'Visible', 'off')
set(handles.uipanel9 , 'Visible', 'off')
set(handles.edit5, 'String', num2str(10));
set(handles.edit4, 'String', num2str(10));

set(handles.checkbox7, 'Value',1);

set(handles.edit3,'String','')

set(handles.checkbox15,'Value',0)
set(handles.edit6,'String','outputs.txt')
set(handles.edit6,'Visible','off')


function graphicalView(hObject, eventdata, handles, Model)
options.showNodeLabel = get(handles.checkbox9,'Value');
options.showElementLabel = get(handles.checkbox10,'Value');
if get(handles.checkbox11,'Value')
    options.color = get(handles.pushbutton12,'BackgroundColor');
else
    if isfield(options,'color')
        options = rmfield(options,'color');
    end
end
    
if get(handles.checkbox12,'Value')
    options.faceColor = get(handles.pushbutton13,'BackgroundColor');
else
    if isfield(options,'faceColor')
        options = rmfield(options,'faceColor');
    end
end

 
options.nodeLabelSize = str2num(get(handles.edit5,'String'));
options.elementLabelSize = str2num(get(handles.edit4,'String'));

options.nodeLabelColor = get(handles.pushbutton15,'BackgroundColor');
options.elemLabelColor = get(handles.pushbutton14,'BackgroundColor');

options.showBoundary = get(handles.checkbox13,'Value');
options.showLoad = get(handles.checkbox14,'Value');

cla(handles.axes1)
reset(handles.axes1)
options.figure = handles.axes1;

Model.graphicalOptions = options;
PlotMeshX(Model)


function graphicalViewDef(hObject, eventdata, handles, Model)
options.showNodeLabel = 0;
options.showElementLabel = 0;
options.color = 'g';
options.faceColor = 'none';
options.showBoundary = 0;
options.showLoad = 0;
options.figure = handles.axes1;
Model.graphicalOptions = options;
PlotMeshX(Model)


function liveGraphics(hObject, eventdata, handles)
fileName = get(handles.edit3,'String');
try
    Model = inpFileReader(fileName);
    graphicalView(hObject, eventdata, handles, Model);
catch
    errordlg('Check the Input File','Error');
end


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web("www.sshahi.com")