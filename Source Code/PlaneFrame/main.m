%%                     M A I N      F U N C T I O N 
%__________________________________________________________________________
% 
%                       Finite Element Methods
%                     Developed by SHAHROKH SHAHI 
%                           (www.sshahi.com)
%
%                   Georgia Institute of Technology
%__________________________________________________________________________
% This file is the main function to handle the PLANE BEAM/FRAME analysis 
% by employing the functions in the folder "lib"

% DEVELOPMENT HISTORY:
%   [2018, Nov, 01] Shahrokh Shahi -- transform to function formats
%   [2018, XXX, XX] Shahrokh Shahi -- 

%% Initial Setup  ( You do NOT need to change this section)
clc
clear
close all

format short g
format compact

% adding "lib" folder to the path
path = mfilename('fullpath');
path(end-length(mfilename):end)=[];
addpath(fullfile(path,'lib'))
addpath(fullfile(path,'lib_io'))

%% Reading the Input File

% input file name:
inpFileName = 'input1.txt';
inpFileName = 'input2.txt';

%--------------------------------------------------------------------------
%      E X T R A C T    M O D E L    F R O M    I N P U T    F I L E 
%--------------------------------------------------------------------------
Model = inpFileReader(inpFileName)


%--------------------------------------------------------------------------
%                 P R I N T    I N P U T   S U M M A R Y 
%--------------------------------------------------------------------------
printSummary(Model)


%--------------------------------------------------------------------------
%                    S I M P L E       P L O T
%--------------------------------------------------------------------------
plotmesh(Model,'b');

%% Run Analysis

Solution1 = planeFrameExplicit(Model);
Solution2 = planeFrameNumerical(Model);

printResult(Solution2)