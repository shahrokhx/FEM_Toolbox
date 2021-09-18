%%                    P A C K A G E     S T A R T E R
%__________________________________________________________________________
% 
%                       Finite Element Package
%                             developed by
%                           SHAHROKH  SHAHI
%                   Georgia Institute of Technology
%                              Spring 2019
%__________________________________________________________________________
% 
% Homepage: www.sshahi.com
% Email: shahi@gatech.edu
%__________________________________________________________________________

%% Initialization
clc
clear 
close all

format compact

path (pathdef)

%% Cleaning
delete *.mat

%% Add to path
pathThis = mfilename('fullpath');
pathThis(end-length(mfilename):end)=[];
% addpath(fullfile(path,'PlaneTruss'             ,'lib'))
% addpath(fullfile(path,'PlaneFrame'             ,'lib'))
% addpath(fullfile(path,'PlaneStress PlaneStrain','lib'))
addpath(fullfile(pathThis,'GUI'));

%% Run GUI
FiniteElementGUI()