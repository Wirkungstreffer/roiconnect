% surface only
function [sourcemodelout, transform] = transform_move_inward(sourcemodel, headmodel, transform)

if ischar(headmodel)
    headmodel = load('-mat', headmodel);
    if isfield(headmodel, 'vol')
        headmodel = headmodel.vol;
        headmodel.unit = 'mm';
    end
end
if ischar(sourcemodel)
    try
        sourcemodel = load('-mat', sourcemodel);
    catch
        % Likely a volume atlas
        sourcemodelout = sourcemodel;
        return
    end
    if isfield(sourcemodel, 'cortex')
        sourcemodel = sourcemodel.cortex;
    end
end
    
newsourcemodel = [];
if isfield(sourcemodel, 'Vertices')
    newsourcemodel.pos = sourcemodel.Vertices;
    newsourcemodel.tri = sourcemodel.Faces;
elseif isfield(sourcemodel, 'vertices')
    newsourcemodel.pos = sourcemodel.vertices;
    newsourcemodel.tri = sourcemodel.faces;
else
    newsourcemodel.pos = sourcemodel.pos;
    newsourcemodel.tri = sourcemodel.tri;
end

cfg = [];
pos = [newsourcemodel.pos ones(size(newsourcemodel.pos,1),1) ];
if ~isempty(transform)
    pos = traditionaldipfit(transform)*pos';
end
pos(4,:) = [];
cfg.sourcemodel.pos = pos';
cfg.sourcemodel.tri = newsourcemodel.tri;
cfg.sourcemodel.unit = headmodel.unit;
cfg.moveinward = 1;
cfg.headmodel = headmodel;
disp('moving source model inward if necessary');
sourcemodelout = ft_prepare_sourcemodel(cfg);
transform = [];

% plot3dmeshalign(headmodel);
% 
% hold on;
% plot3dmeshalign(tmp2, [], [1 0 0])
