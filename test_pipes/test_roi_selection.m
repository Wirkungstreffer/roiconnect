% Test connectivity estimation on a selection of regions (ROIs).
%% Run pipeline
clear
eeglab

eeglabp = fileparts(which('eeglab.m'));
EEG = pop_loadset('filename','eeglab_data_epochs_ica.set','filepath',fullfile(eeglabp, 'sample_data/'));
EEG = pop_resample( EEG, 100);
EEG = pop_epoch( EEG, { }, [-0.5 1.5], 'newname', 'EEG Data epochs epochs', 'epochinfo', 'yes');
EEG = pop_select( EEG, 'trial',1:30);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
eeglab redraw;

EEG = pop_dipfit_settings( EEG, 'hdmfile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','standard_vol.mat'), ...
    'coordformat','MNI','mrifile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','standard_mri.mat'),...
    'chanfile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','elec', 'standard_1005.elc'),...
    'coord_transform',[0.83215 -15.6287 2.4114 0.081214 0.00093739 -1.5732 1.1742 1.0601 1.1485] ,'chansel',[1:32] );

EEG = pop_leadfield(EEG, 'sourcemodel',fullfile(eeglabp,'functions','supportfiles','head_modelColin27_5003_Standard-10-5-Cap339.mat'), ...
    'sourcemodel2mni',[0 -24 -45 0 0 -1.5708 1000 1000 1000] ,'downsample',1);
% EEG = pop_leadfield(EEG, 'sourcemodel',fullfile(eeglabp,'plugins','dipfit','LORETA-Talairach-BAs.mat'), ...
%     'sourcemodel2mni',[0 -24 -45 0 0 -1.5708 1000 1000 1000] ,'downsample',1);

EEG = pop_roi_activity(EEG, 'leadfield',EEG.dipfit.sourcemodel,'model','LCMV','modelparams',{0.05},'atlas','LORETA-Talairach-BAs','nPCA',3, 'chansel', EEG.dipfit.chansel);

%% Test PAC for a selection of ROIs
low = 10;
high = 50;

fcomb.low = low;
fcomb.high = high;
roi_selection = {1, 5, 11};

tic
EEG1 = pop_roi_connect(EEG, 'methods', {'PAC'}, 'fcomb', fcomb); % default option (PAC computed on all ROIs)
toc

tic
EEG2 = pop_roi_connect(EEG, 'methods', {'PAC'}, 'fcomb', fcomb, 'roi_selection', roi_selection); % PAC computed on a selection of ROIs
toc

%% Test if PAC values are (exactly) the same
% add 'test_pipes/cmp_roi_selection' to path first
cmp_roi_selection(roi_selection, EEG1.roi.PAC.b_orig_norm, EEG2.roi.PAC.b_orig_norm)

%% Test FC options for a selection of ROIs
% Tested for: MIM, MIC, CS, wPLI, GC, TRGC, COH, PDC, TRPDC, DTF, TRDTF
roi_selection = {1, 5, 11};
tic
EEG3 = pop_roi_connect(EEG, 'methods', {'MIM'}); % default option (MIM computed on all ROIs)
toc

tic
EEG4 = pop_roi_connect(EEG, 'methods', {'MIM'}, 'roi_selection', roi_selection); % default option (MIM computed on all ROIs)
toc

%% Test if FC values are (exactly) the same (for some reason, MATLAB introduces some rounding errors sometimes..)
% add 'test_pipes/cmp_roi_selection' to path first
cmp_roi_selection(roi_selection, squeeze(EEG3.roi.MIM(10, :, :)), squeeze(EEG4.roi.MIM(10, :, :))) 

%% Plot connectivity 
pop_roi_connectplot(EEG3, 'measure', 'mim', 'plotcortex', 'off', 'plotmatrix', 'on', 'freqrange', [8 13]);
pop_roi_connectplot(EEG3, 'measure', 'mim', 'plotcortex', 'off', 'plotmatrix', 'on', 'freqrange', [8 13], 'region', 'central'); % region option disabled

pop_roi_connectplot(EEG4, 'measure', 'mim', 'plotcortex', 'off', 'plotmatrix', 'on', 'freqrange', [8 13]);
pop_roi_connectplot(EEG4, 'measure', 'mim', 'plotcortex', 'off', 'plotmatrix', 'on', 'freqrange', [8 13], 'region', 'central'); % region option disabled