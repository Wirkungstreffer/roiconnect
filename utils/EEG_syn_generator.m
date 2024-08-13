% Script to specifically load and process the head model procedural from Brainstorm. 
% Based on pop_importdata.m.
%
% Inputs:
%   signal_roi     -  (n_roi x epleng x n_epochs) simulated ROI data
%   signal_sensor  -  (n_channel x epleng x n_epochs) simulated sensor data
%   fs             -  [integer] sampling frequency
%   source         -  ([epleng * n_epochs] x [n_univ + n_biv]) source data containing source PAC
% 
% Output:
%   EEGOUT         - modified EEG dataset structure
%   

function [EEG_syn, com] = EEG_syn_generator(signal_roi, signal_sensor, fs, varargin)
    % load the forward solution of bs_result.mat
    load bs_results.mat

%% Signal generation
    g = finputcheck(varargin, { 
        'dataformat'        'string'             { }              'array';
        'nbchan'            'integer'            { }              0;
        'setname'           'string'             { }              'PAC_Simulation';
        'pnts'              'integer'            { }              0;
        'xmin'              'float'              { }              0;
        'nPCA'              'integer'            { }              1;
    }, 'EEG_syn_generator'); 


    EEG_syn = pop_importdata('dataformat',g.dataformat,'nbchan',g.nbchan,'data',signal_sensor,'setname',g.setname,'srate',fs,'pnts',g.pnts,'xmin',g.xmin);

    % add EEG fields
    EEG_syn.roi.source_roi_data = signal_roi;
    EEG_syn.roi.nPCA = g.nPCA;
    EEG_syn.roi.nROI =  size(signal_roi, g.nPCA);

    EEG_syn.roi.cortex = cortex;
    EEG_syn.roi.atlas = cortex.Atlas(3);

    % Access the required structure
    atlas = cortex.Atlas(3);
    
    % Extract the fields into a cell array
    fields = fieldnames(atlas.Scouts);
    
    % Identify the columns to keep and their order
    keepFields = {'Label', 'Vertices'};
    
    % Reorder the fields
    [~, idx] = ismember(keepFields, fields);
    reorderedFields = fields(idx);
    
    % Create a new structure with only the required fields
    newScouts = rmfield(atlas.Scouts, setdiff(fields, keepFields));
    newScouts = orderfields(newScouts, reorderedFields);
    
    % Replace the original structure with the new one
    EEG_syn.roi.atlas.Scouts = newScouts;
    EEG_syn.roi.cortex.Atlas = EEG_syn.roi.atlas;
    
    fres = fs; 
    freqs = sfreqs(fres, fs); % freqs in Hz
    EEG_syn.roi.freqs = freqs;
    
    if ~isfield(EEG_syn.roi, 'pnts')
        EEG_syn.roi.pnts = EEG_syn.pnts;
    end
    
    if ~isfield(EEG_syn.roi, 'srate')
        EEG_syn.roi.srate = EEG_syn.srate;
    end

