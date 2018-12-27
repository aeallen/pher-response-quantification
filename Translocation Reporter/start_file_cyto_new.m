clear;
clc;
close all;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% specify folders , and name. 
% pick one file in each folder, type in their names
% exclude the position names and channel numbers
% run('');
%% 45p 15b
i=0;
% %

i = i+1;
Folders{i} = '';
BaseFileNameExample{i} = '';
% i = i+1;
% Folders{i} = '';
% BaseFileNameExample{i} = '';
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% specify channel numbers
% the struct has three fields: channelName, channelNum
% user need to change them according to their FPs
usedChName{1} = 'Phase';
usedChNumb{1} = 1;

usedChName{2} = '';
usedChNumb{2} = ;

usedChName{3} = '';
usedChNumb{3} = ;

usedChName{4} = '';
usedChNumb{4} = ;



%% which channel you want to measure, clearly most of you don't want the phase image
channels2measure = {''};

%% specify which channel is nucleos
nuc_ch_name = '';

phase_ch_name = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% some other stuff
%% test with just one file or not
oneFileTest = false;

%% specify the tracking direction
trackDirection = 'backward'; % or 'backward'

%% registor image or not
registerOrNot = true; % choose true if there is big stage shift, usually not

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save intermidiate results or not, if use v2, they are not saved
saveLabelAndTraj = true;

%% load labels and skip segmentation
skip_seg = true;

%% load traj and skip tracking
% does nothing when choose parallel computation
skip_track = true;

%% load intensities and skip measuring
% does nothing when choose parallel computation
skip_measure = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run('~/yeast_analysis_nuc_cyto_v1.m');
toc




