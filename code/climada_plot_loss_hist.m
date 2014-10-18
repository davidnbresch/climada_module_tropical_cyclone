function climada_plot_loss_hist(event_loss, nametag)
% TC footprint figure
% NAME:
%   climada_plot_windfield
% PURPOSE:
% create footprint figure
% CALLING SEQUENCE:
%   [contr t_handle] = climada_plot_windfield(hazard, tc_track, track_no)
% EXAMPLE:
%   climada_plot_windfield
% INPUTS:
%   hazard: hazard.intensity with wind intensities per centroid
%   tc_track: a structure with the track information:
%   track_no: number of track to show footprint
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   figure with footprint
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Lea Mueller, david.bresch@gmail.com, 20121205
%-

global climada_global
if ~climada_init_vars, return; end
if ~exist('event_loss'   , 'var'), event_loss   = []; end
if ~exist('nametag'      , 'var'), nametag      = ''; end


if any(full(event_loss))
    event_loss = full(event_loss(event_loss>0));
    max_loss = max(event_loss);
    min_loss = min(event_loss);
    cmap     = get(gcf,'colormap');
    no_edges = size(cmap,1);
    edges    = logspace(log10(min_loss),log10(max_loss),no_edges);
    no       = hist(event_loss, edges);
    %bar(xo,no)
    for ii = 1:no_edges-1    
        patch([edges(ii) edges(ii) edges(ii+1) edges(ii+1)], ...
              [0 no(ii) no(ii) 0],...
              cmap(ii,:),'edgecolor','white');
    end
    set(gca, 'Xscale', 'log','layer','top');
    ylim([0 max(no)*1.2])
else
    no = 0;
end
ylabel('Number of pixels')
xlabel('Loss (USD)')
titlestr = sprintf('%d loss pixels, total loss: %2.1G',sum(no),full(sum(event_loss)));
title(titlestr)

end
    
