% Test PAC implementation for different inputs: single frequencies and frequency bands.
%% Run pipeline
clear
eeglab

% eeglabp = fileparts(which('eeglab.m'));
% EEG_orgin = pop_loadset('filename','eeglab_data_epochs_ica.set','filepath',fullfile(eeglabp, 'sample_data/'));
% EEG_orgin = pop_resample( EEG_orgin, 100);
% EEG_orgin = pop_epoch( EEG_orgin, { }, [-0.5 1.5], 'newname', 'EEG_orgin Data epochs epochs', 'epochinfo', 'yes');
% EEG_orgin = pop_select( EEG_orgin, 'trial',1:30);
% % [ALLEEG_orgin, EEG_orgin, CURRENTSET] = eeg_store(ALLEEG_orgin, EEG_orgin);
% eeglab redraw;
% 
% EEG_orgin = pop_dipfit_settings( EEG_orgin, 'hdmfile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','standard_vol.mat'), ...
%     'coordformat','MNI','mrifile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','standard_mri.mat'),...
%     'chanfile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','elec', 'standard_1005.elc'),...
%     'coord_transform',[0.83215 -15.6287 2.4114 0.081214 0.00093739 -1.5732 1.1742 1.0601 1.1485] );
% 
% % EEG_orgin = pop_leadfield(EEG_orgin, 'sourcemodel',fullfile(eeglabp,'plugins','dipfit','LORETA-Talairach-BAs.mat'), ...
% %     'sourcemodel2mni',[0 -24 -45 0 0 -1.5708 1000 1000 1000] ,'downsample',1);
% 
% EEG_orgin = pop_leadfield(EEG_orgin, 'sourcemodel',fullfile(eeglabp,'functions','supportfiles','head_modelColin27_5003_Standard-10-5-Cap339.mat'), ...
%     'sourcemodel2mni',[0 -24 -45 0 0 -1.5708 1000 1000 1000] ,'downsample',1);
% 
% EEG_orgin = pop_roi_activity(EEG_orgin, 'leadfield',EEG_orgin.dipfit.sourcemodel,'model','LCMV','modelparams',{0.05},'atlas','LORETA-Talairach-BAs','nPCA', 3, 'chansel', EEG_orgin.dipfit.chansel);

%% Test simulation inputs
% addpath /home/ben/eeglab/plugins/roiconnect/simulations

% roi_index1 = 27; %lingual L
% roi_index2 = 59; % superiorparietal L

% low = [5 7]; %in Hz
% high = [30 34]; %in Hz
sim_case = 3; %   sim_case - [1, 2, 3] univariate, bivariate or univariate + bivariate
n_univ = 2;
n_biv = 2;
isnr = 0.5;
% iroi_pac = [4 4];
% iroi_pac = [[13 19]; [22 28]];
iroi_pac = [[4 4]; [13 19]; [22 28];[45 45]]; 
% add functions not only PAC, but also other connectivity methods in between ROIs. 
% Add the time and frequency interaction of different ROIs.

% iroi_pac = [];

[signal_roi, signal_sensor, fs, source, filt, L] = sim_wholebrain_pac(sim_case, n_univ, n_biv, isnr, iroi_pac);

EEG_syn = EEG_syn_generator(signal_roi, signal_sensor, fs);


%% Test bispectrum for single frequency inputs
low = 8;
high = 32;

if any(low > fs/2)
    error('All elements of low must be less than or equal to fs/2.');
end

if any(high > fs/2)
    error('All elements of high must be less than or equal to fs/2.');
end

fcomb.low = low;
fcomb.high = high;

EEG_syn_1 = pop_roi_connect(EEG_syn, 'methods', {'PAC'}, 'fcomb', fcomb); % test all 3 connectivity functions (data2spwctrgc, data2strgcmim, roi_pac)
% EEG_orgin_1 = pop_roi_connect(EEG_orgin, 'methods', {'PAC', 'MIM', 'COH'}, 'fcomb', fcomb); % test all 3 connectivity functions (data2spwctrgc, data2strgcmim, roi_pac)
tic
EEG_syn_2 = pop_roi_connect(EEG_syn, 'methods', {'PAC'}, 'fcomb', fcomb, 'conn_stats', 'on', 'nshuf', 3, 'poolsize', 12); % compute only b_anti, b_anti_norm
toc
% EEG_syn_3 = pop_roi_connect(EEG_syn, 'methods', {'PAC'}, 'fcomb', fcomb, 'bs_outopts', 5); 

% test if first shuffle equals the true value
tolerance = 1e-7; 
if ~all(ismembertol(EEG_syn_1.roi.PAC.b_orig, squeeze(EEG_syn_2.roi.PAC.b_orig(:, :, 1)), tolerance), 'all')
    error 'The first shuffle in the surrogate connectivity array is not the true matrix.'
end

if ~all(ismembertol(EEG_syn_1.roi.PAC.b_anti, squeeze(EEG_syn_2.roi.PAC.b_anti(:, :, 1)), tolerance), 'all')
    error 'The first shuffle in the surrogate connectivity array is not the true matrix.'
end

%% Test bispectrum for frequency band inputs
low = [4 8];
high = [48 50];

if any(low > fs/2)
    error('All elements of low must be less than or equal to fs/2.');
end

if any(high > fs/2)
    error('All elements of high must be less than or equal to fs/2.');
end

fcomb.low = low;
fcomb.high = high;

tic
EEG_syn_4 = pop_roi_connect(EEG_syn, 'methods', {'PAC'}, 'fcomb', fcomb); % test all 3 connectivity functions (data2spwctrgc, data2strgcmim, roi_pac)toc
toc
% tic
% EEG_syn_5 = pop_roi_connect(EEG_syn, 'methods', {'PAC'}, 'fcomb', fcomb, 'conn_stats', 'off', 'nshuf', 3, 'bs_outopts', 5); 
% toc

% test if first shuffle equals the true value
tolerance = 1e-7; 
if ~all(ismembertol(EEG_syn_4.roi.PAC.b_anti, squeeze(EEG_syn_5.roi.PAC.b_anti(:, :, 1)), tolerance), 'all')
    error 'The first shuffle in the surrogate connectivity array is not the true matrix.'
end

%% Test PAC plotting
% Test for single frequency inputs
pop_roi_connectplot(EEG_syn_1, 'measure', 'PAC', 'plotmatrix', 'on');
% pop_roi_connectplot(EEG_syn_1, 'measure', 'PAC', 'plotmatrix', 'on', 'plotcortex', 'off');
pop_roi_connectplot(EEG_syn_1, 'measure', 'PAC_anti', 'plotmatrix', 'on');

% Provoke errors by plotting bispectral tensors that do not exist
pop_roi_connectplot(EEG_syn_2, 'measure', 'PAC_anti', 'plotmatrix', 'on'); % bs_outopts 4 means only original bispectra are computed, so cannot plot anti
pop_roi_connectplot(EEG_syn_2, 'measure', 'PAC', 'plotmatrix', 'on'); % bs_outopts 5 means only antisymm. bispectra are computed, so cannot plot normal bispectrum

% Test for frequency bands
pop_roi_connectplot(EEG_syn_4, 'measure', 'PAC', 'plotmatrix', 'on');
pop_roi_connectplot(EEG_syn_4, 'measure', 'PAC_anti', 'plotmatrix', 'on');

% plot MIM and aCOH as a sanity check
pop_roi_connectplot(EEG_syn_1, 'measure', 'MIM', 'plotmatrix', 'on');
pop_roi_connectplot(EEG_syn_1, 'measure', 'aCOH', 'plotmatrix', 'on');

% Statistic test plot
pop_roi_statsplot(EEG_syn_2, 'measure', 'PAC', 'bispec', 'b_anti');
pop_roi_statsplot(EEG_syn_2, 'measure', 'PAC', 'bispec', 'b_orig');
