%% Load and project Lead Field
restoredefaultpath;
clc;
clear all;
close all;
tic
addpath('common_functions');
addpath(genpath('D:\Toolsboxs\brainstorm3'));

preproced_data_file = 'E:\Data\Conectome\175237\meg\175237\MEG\Restin\rmegpreproc\175237_MEG_5-Restin_rmegpreproc.mat';
preproced_data = load(preproced_data_file);

protocol_path = 'C:\Users\Ariosky\.brainstorm\local_db\BigBrain_Template';

leadfield_MEG_file = fullfile(protocol_path,'data','ICBM152_MEG','175237_MEG_5-Restin_rmegpreproc','headmodel_surf_os_meg');
leadfield_MEG = load(leadfield_MEG_file);

MEG_channel_file = fullfile(protocol_path,'data','ICBM152_MEG','175237_MEG_5-Restin_rmegpreproc','channel_4d_acc1.mat');
MEG_channel = load(MEG_channel_file);

disp('-->> Creating channel');
% [L_MEG3D,MEG_channel] = remove_leadfield_channel(MEG_channel,leadfield_MEG,preproced_data);



GridOrient = leadfield_MEG.GridOrient;
GridAtlas = leadfield_MEG.GridAtlas;

clearvars leadfield_MEG

L_MEG = bst_gain_orient(L_MEG3D, GridOrient, GridAtlas);

surface_file = fullfile(protocol_path,'anat','ICBM152_MEG','tess_cortex_pial_8000V.mat');
surface = load(surface_file);

%% Load MEG trials in single file 
disp('-->> Loading MEG trials in single file');
Ntpoints    = size(preproced_data.data.trial{1,1},2);
Nsegments   = length(preproced_data.data.trial);
Data        = zeros(size(L_MEG,1),Ntpoints,Nsegments);
Data10Hz    = zeros(size(L_MEG,1),Ntpoints,Nsegments);
Fs          = 508.6275;
deltaf      = Fs/Ntpoints;
F           = 0:deltaf:(Ntpoints*deltaf);

Svv         = 0;
Svv10Hz     = 0;
for count = 1:Nsegments
    disp(strcat('-->> loading and filtering segments: ', num2str(count)));
    tmp       = preproced_data.data.trial{1,count};
    tmp10Hz   = bandpass(tmp,[8 12],Fs);
    Svv       = Svv + tmp*tmp'/Ntpoints;
    Svv10Hz   = Svv10Hz + tmp10Hz*tmp10Hz'/Ntpoints;
    Data(:,:,count) = tmp;
    Data10Hz(:,:,count) = tmp10Hz;
end
clearvars preproced_data MEG_channel

Svv     = Svv/Nsegments;
Svv10Hz = Svv10Hz/Nsegments;
disp('-->> Computing lcmv');
%% Inverse solution filter
p                   = size(L_MEG,1);
Ip                  = eye(p);
scaleL_MEG          = sqrt(sum(abs(diag(L_MEG*L_MEG')))/p);
L_MEG               = L_MEG/scaleL_MEG;

Svv0                = Svv10Hz;
scaleV              = (sum(abs(diag(Svv0)))/p);
Svv0                = Svv0/scaleV;

gamma               = sum(abs(diag(Svv0)))/(length(Svv0)*100);
[T_MEG,T1_MEG,Wout] = mkfilt_lcmv(L_MEG,Svv0,gamma);
T_MEG               = T_MEG';
%% Source time series 
Data0   = Data10Hz;
Source  = zeros(size(L_MEG,2),size(Data,2),size(Data,3));
for count = 1:size(Source,3)    
    Source(:,:,count) = T_MEG*Data0(:,:,count);
end

clearvars Data Data0 Data10Hz

%% Saving trials 

for count3 = 1:Nsegments
    disp(strcat("-->> Saving trial: ",num2str(count3)));
    trial = squeeze(Source(:,:,count3));
    file_name = strcat('trials\trial_',num2str(count3),'.mat');
    save( file_name, 'trial', '-v7.3');
end
clearvars Source;