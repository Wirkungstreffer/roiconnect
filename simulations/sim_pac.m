% Pipeline to run the decomposition method on simulated PAC data with a
% specified number of univariate and bivariate interactions. 
%
% Input:
%   n_shuf - [integer] number of shuffles
%   sim_case - [1, 2, 3] univariate, bivariate or univariate + bivariate
%   n_univ   - [integer] number of univariate interactions, default is 1
%   n_biv    - [integer] number of bivariate interactions, default is 2
%   isnr     - [float] total signal-to-noise ratio
%   iroi_pac - [integer] array of coupled regions, where the first column corresponds to the phase region and the second column to the amplitude
%              region. For example, if [30 40], then region 30 will be the phase region and 40 will be the amplitude region. If [[4 5]; [8 9]], 
%              then [4 8] will be the phase regions and [5 9] the amplitude regions. If left empty, the univariate/bivariate interactions are randomly placed.
%
% Optional inputs:
%   alpha  - [float] significance level, default is 0.05.

function [signal_roi, signal_sensor, fs, source, filt, L] = sim_pac(sim_case, n_univ, n_biv, isnr, iroi_pac, varargin)
    
    % local_path = '/home/ben/eeglab/plugins/roiconnect/';
    % DIROUT = [local_path 'simulations/sim_pac/figures/'];
    % 
    % if ~exist(DIROUT, 'dir')
    %     mkdir(DIROUT)
    % end

    % setup
    % eeglab

    %% Simulate one bivariate PAC interaction between an occopital and a parietal region
    roi_idx1 = 27; %lingual L
    roi_idx2 = 59; % superiorparietal L

    g = finputcheck(varargin, { 
        'alpha'          'float'         { }              0.05;...
        'sim_case'      'integer'        { }              2;
        'n_univ'        'integer'        { }              0; ...
        'n_biv'         'integer'        { }              1;
        'isnr'          'float'          { }              0.5 ;
        'iroi_pac'      'integer'        { }              [roi_idx1 roi_idx2];
        }, 'sim_pac'); 
    if ischar(g), error(g); end


    if exist('signal_roi', 'var') && exist('sim_case', 'var') && exist('n_univ', 'var') && exist('n_biv', 'var') && exist('isnr', 'var') && exist('iroi_pac', 'var')
        % If inputs are given
        [signal_roi, signal_sensor, fs, source, filt, L] = sim_wholebrain_pac(sim_case, n_univ, n_biv, isnr, iroi_pac);
    else
        % If inputs are not given
         [signal_roi, signal_sensor, fs, source, filt, L] = sim_wholebrain_pac(g.sim_case, g.n_univ, g.n_biv, g.isnr, g.iroi_pac);
    end
    %[signal_sensor, fs, source, filt, L] = sim_wholebrain_pac(g.sim_case, g.n_univ, g.n_biv, g.isnr, g.iroi_pac);
    %[signal_sensor, fs, source, filt, L] = sim_wholebrain_pac(sim_case, n_univ, n_biv, isnr, iroi_pac);
    
    % % sampling frequency
    % fres = fs; 
    % frqs = sfreqs(fres, fs); % freqs in Hz
    % 
    % % set parameter values for (cross)-bispectrum estimation
    % freqinds = [mean(filt.low) mean(filt.high)]; % in Hz
    % len_epochs = 2; % 2-second epochs
    % segleng = fs * len_epochs; 
    % segshift = floor(segleng/2);
    % epleng = fs * len_epochs; % create epochs of [e.epleng] seconds 
    % 
    % % run decomposition without antisymmetrization
    % n = [1 2 3];
    % err_colors = ['r', 'b', 'k'];
    % [P_source_fdr, P_source, F, F_moca, A_hat, A_demixed, D_hat, D_demixed, errors] = bsfit_stats(signal_sensor, freqinds(1), ...
    %     freqinds(2), n, n_shuf, frqs, segleng, segshift, epleng, g.alpha, L);
    % plot_error(errors, 1, n, err_colors, '', '', DIROUT)

end