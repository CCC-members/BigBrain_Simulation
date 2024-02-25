%% Load and project Lead Field
restoredefaultpath;
clc;
clear all;
close all;
tic
addpath('common_functions');
addpath('D:\Toolsboxs\brainstorm3');

preproced_data_file = 'E:\Data\Conectome\175237\meg\175237\MEG\Restin\rmegpreproc\175237_MEG_5-Restin_rmegpreproc.mat';
preproced_data = load(preproced_data_file);

protocol_path = 'C:\Users\Ariosky\.brainstorm\local_db\Simulation';

leadfield_MEG_file = fullfile(protocol_path,'data','175237_MEG_4D_Neuroimaging','@raw5-Restin_c_rfDC','headmodel_surf_os_meg');
leadfield_MEG = load(leadfield_MEG_file);

MEG_channel_file = fullfile(protocol_path,'data','175237_MEG_4D_Neuroimaging','@raw5-Restin_c_rfDC','channel_4d_acc1.mat');
MEG_channel = load(MEG_channel_file);

disp('-->> Creating channel');
[L_MEG3D,MEG_channel] = remove_leadfield_channel(MEG_channel,leadfield_MEG,preproced_data);

GridOrient = leadfield_MEG.GridOrient;
GridAtlas = leadfield_MEG.GridAtlas;

clearvars leadfield_MEG

% L_MEG = bst_gain_orient(L_MEG3D, GridOrient, GridAtlas);
L_MEG = L_MEG3D;

surface_file = fullfile(protocol_path,'anat','175237_MEG_4D_Neuroimaging','tess_cortex_concat_8000V.mat');
surface = load(surface_file);

%% Data spectral analysis
disp('-->> Loading MEG trials in single file');
Ntpoints    = size(preproced_data.data.trial{1,1},2);
Nsegments   = length(preproced_data.data.trial);
Data_FC     = zeros(size(L_MEG,1),Ntpoints,Nsegments);
Fs          = preproced_data.data.fsample;
deltaf      = Fs/Ntpoints;
F           = 0:deltaf:(Ntpoints*deltaf);
F1          = 8;
F2          = 12;
[val,idf1]   = min(abs(F-F1));
[val,idf2]   = min(abs(F-F2));
Svv_MEG         = 0;
for count3 = 1:Nsegments
    disp(strcat('-->> loading and filtering segments: ', num2str(count3)));
    tmp                 = preproced_data.data.trial{1,count3};
    tmp_FC              = fft(tmp,[],2);
    Data_FC(:,:,count3) = tmp_FC;
    Svv_MEG                 = Svv_MEG + tmp_FC(:,idf1:idf2)*tmp_FC(:,idf1:idf2)';
end
Nsamples = Nsegments*(idf2-idf1);
Svv_MEG = Svv_MEG/Nsamples;
%% Inverse solution filter
nonovgroups           = [];
counter               = 1;
for ii = 1:3:length(L_MEG)
    nonovgroups{counter} = [ii;ii+1;ii+2];
    counter = counter + 1;
end
[miu,sigma_post,T_MEG3D] = cross_nonovgrouped_enet_ssbl({Svv_MEG},{L_MEG},Nsamples,nonovgroups);
miu_stat              = sqrt(sum(reshape(miu,3,length(L_MEG)/3),1))';
sigma_post_stat       = sqrt(sum(reshape(sigma_post,3,length(L_MEG)/3),1))';
stat                  = abs(miu_stat)./abs(sigma_post_stat);
indms                 = find(stat > 1);
J                     = zeros(length(L_MEG)/3,1);
J(indms)              = stat(indms); J = J/max(abs(J));
%%
figure;
title('HCP Subject #175237');
patch('Faces',surface.Faces,'Vertices',surface.Vertices,'FaceVertexCData',J,'FaceColor','interp',...
    'EdgeColor',[0.313725501298904 0.313725501298904 0.313725501298904],'FaceAlpha',.99);
axis equal;

OptionZ.FrameRate=25;
OptionZ.Duration=50;
OptionZ.Periodic=true;
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], [pwd,filesep,'HCP_Subject_175237'],OptionZ);
toc
disp('done');


