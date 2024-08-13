% Script to simulate whole-brain signals with univariate and bivariate PAC interactions. 
% Based on Franziska Pellegrini's PAC repository https://github.com/fpellegrini/PAC/blob/master/fp_pac_sim.m.
%
% Inputs:
%   sim_case - [1, 2, 3] univariate, bivariate or univariate + bivariate
%   n_univ   - [integer] number of univariate interactions
%   n_biv    - [integer] number of bivariate interactions
%   isnr     - [float] total signal-to-noise ratio
%   iroi_pac - [integer] array of coupled regions, where the first column corresponds to the phase region and the second column to the amplitude
%              region. For example, if [30 40], then region 30 will be the phase region and 40 will be the amplitude region. If [[4 5]; [8 9]], 
%              then [4 8] will be the phase regions and [5 9] the amplitude regions. If left empty, the univariate/bivariate interactions are randomly placed.
%   iReg     - number of active voxels per region
%   iss      - brain noise to sensor noise ration (19 dB)
%   t        - [0, 1] true voxel pipeline or not
%   fs       - frequency rate 
%   epleng   - epoch length in sec
%   N        - total sample number
%   low      - the lower frequency of variate
%   high     - the higher frequency of variate
% Output:
%   signal_roi - (n_roi x epleng x n_epochs) simulated ROI data
%   fs         - [integer] sampling frequency
%   source     - ([epleng * n_epochs] x [n_univ + n_biv]) source data containing source PAC
%   filt       - filter settings

function [signal_roi, signal_sensor, fs, sources, filt, L] = sim_wholebrain_pac(sim_case, n_univ, n_biv, isnr, iroi_pac, varargin)
%% Signal generation
    g = finputcheck(varargin, { 
        %'alpha'         'float'          { }              0.05;...
        %'sim_case'      'integer'        { }              2;
        %'n_univ'        'integer'        { }              0; ...
        %'n_biv'         'integer'        { }              1;
        %'isnr'          'float'          { }              0.5 ;
        %'iroi_pac'      'integer'        { }              [27 59]; % [lingual L superiorparietal L]
        'iReg'          'integer'        { }              1;
        'iss'           'float'          { }              0.9;
        't'             'integer'        { }              0;
        'fs'            'integer'        { }              100;
        'epleng'        'float'          { }              2;
        'N'             'integer'        { }              12000;
        'low'           'float'          { }              [5 7];
        'high'          'float'          { }              [30 34];
    }, 'sim_wholebrain_pac'); 
    if ischar(g), error(g); end    

    if isempty(iroi_pac)
        params.iroi_phase = [];
        params.iroi_amplt = [];
    else
        warning('You have manually specified interacting regions. Make sure that the number of interactions, the case, and the number of passed ROIs coincide.')
        params.iroi_phase = iroi_pac(:, 1);
        params.iroi_amplt = iroi_pac(:, 2);
    end

    % Error checks for simulation cases and number of interactions
    if (sim_case == 1 || sim_case == 3) && n_univ == 0
        warning('The number of univariate interactions does not coincide with the simulation case.');
    end
    if (sim_case == 2 || sim_case == 3) && n_biv == 0
        warning('The number of bivariate interactions does not coincide with the simulation case.');
    end

    % Check the consistency of iroi_pac with n_univ and n_biv
    if ~isempty(iroi_pac)
        num_univ = 0;
        num_biv = 0;

        % Iterate through each row of roi_pac
        for i = 1:size(iroi_pac, 1)
            if iroi_pac(i, 1) == iroi_pac(i, 2)
                % If both elements in the row are equal, increment univariate counter
                num_univ = num_univ + 1;
            else
                % If elements are different, increment bivariate counter
                num_biv = num_biv + 1;
            end
        end
        
        if num_univ ~= n_univ
                warning('The number of univariate interactions in the input iroi_pac does not coincide with the input n_univ.');
        end
        if num_biv ~= n_biv
            warning('The number of bivariate interactions in the input iroi_pac does not coincide with the input n_biv.');
        end
    end

    % set parameters (see fp_pac_signal.m for documentation)
    params.case = sim_case; % univariate + bivariate case
    params.iReg = g.iReg; % number of active voxels per region
    params.iss = g.iss; % brain noise to sensor noise ration (19 dB)
    params.isnr = isnr; % total SNR 
    params.t = g.t; % [0, 1] true voxel pipeline or not
    params.fs = g.fs;
    params.epleng = g.epleng; % Epoch length in sec
    params.N = g.N;
    params.low = g.low;
    params.high = g.high;

    % varify the simulation case
    if params.case == 3
        assert(~any([n_univ n_biv] == 0), 'Indicate the number of uni and bivariate interactions in the mixed case.')
        params.iInt = [n_univ n_biv]; % number of univariate interactions, number of bivariate interactions
    else
        params.iInt = sum([n_univ n_biv]); 
    end

    % high and low <= fs/2 

    % To do :validating the low and high = epleng * fs
    % if (abs(fcomb.low(:,2) - fcomb.low(:,1)) * abs(fcomb.high(:,2) - fcomb.high(:,1))) * fs ~= N
    %     warning('Make sure your data is coincide.')

    % get atlas, voxel and roi indices; active voxel of each region is aleady selected here
    iReg = 1; 
    fprintf('Getting atlas positions... \n')
    D = fp_get_Desikan(iReg); 

    % signal generation
    fprintf('Signal generation... \n')
    [sig, brain_noise, sensor_noise, L, iroi_phase, iroi_amplt, D, fs, n_trials, filt, sources] = fp_pac_signal(params, D);

    % combine noise sources
    noise = params.iss*brain_noise + (1-params.iss)*sensor_noise;
    noise = noise ./ norm(noise(:),'fro');

    % combine signal and noise
    signal_sensor1 = params.isnr*sig + (1-params.isnr)*noise;
    signal_sensor1 = signal_sensor1 ./ norm(signal_sensor1(:), 'fro');

    % high-pass filter signal
    signal_sensor = (filtfilt(filt.bhigh, filt.ahigh, signal_sensor1'))';
    signal_sensor = signal_sensor / norm(signal_sensor, 'fro');
    
    % reshape (create epochs)
    signal_sensor = reshape(signal_sensor,[],size(signal_sensor,2)/n_trials,n_trials);
    [n_sensors, l_epoch, n_trials] = size(signal_sensor);

    %% Source reconstruction

    % select only voxels that belong to any roi
    L_backward = L(:, D.ind_cortex, :);
    A = fp_get_lcmv(signal_sensor, L_backward);

    % dimensionality reduction
    signal_roi = fp_dimred(signal_sensor,D,A,params.t);
end