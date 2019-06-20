%function climada_tc_damage(param1,param2)
% climada tc damage calculation
% MODULE:
%   _LOCAL
% NAME:
%   climada_tc_damage
% PURPOSE:
%   given a single track file and an entity, calculate damage and plot it
%
% CALLING SEQUENCE:
%   climada_tc_damage
% EXAMPLE:
%   climada_tc_damage
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20190504, initial
%-

%global climada_global
if ~climada_init_vars,return;end % init/import global variables

% PARAMETERS
markersize=1;

tc_track=climada_tc_read_unisys_track;
entity=climada_entity_load;
hazard=climada_tc_hazard_set(tc_track,'NOSAVE',entity);
% poor man's encoding:
entity.assets.hazard.filename='NOSAVE';
entity.assets.hazard.comment='on the fly';
entity.assets.centroid_index=1:length(entity.assets.lon);

EDS=climada_EDS_calc(entity,hazard);

%climada_entity_plot(entity_1km,1);hold on;
climada_hazard_plot(hazard,-1,1);hold on;
plot(tc_track.lon,tc_track.lat,'-k');hold on;
entity_dmg=entity; % make a copy to plot damage
entity_dmg.assets.Value=EDS.ED_at_centroid';
entity_dmg.assets.Value(entity_dmg.assets.Value<=1)=NaN;
params.unit_scale=1e6;params.cbar_ylabel='damage';
climada_entity_plot(entity_dmg,1,params);hold on;
title_str=sprintf('%s, damage %3.2g USD',tc_track.name,EDS.ED);
title(title_str);box on

hazard_ts=climada_ts_hazard_set(hazard,'NOSAVE');
EDS(2)=climada_EDS_calc(entity,hazard_ts);

figure('Name','storm surge');
%climada_entity_plot(entity_1km,markersize);hold on;
climada_hazard_plot(hazard_ts,-1,markersize);hold on;
plot(tc_track.lon,tc_track.lat,'-k');hold on;
entity_dmg=entity; % make a copy to plot damage
entity_dmg.assets.Value=EDS(2).ED_at_centroid';
entity_dmg.assets.Value(entity_dmg.assets.Value<=1)=NaN;
params.unit_scale=1e6;params.cbar_ylabel='damage';
climada_entity_plot(entity_dmg,markersize,params);hold on;
title_str=sprintf('%s, surge damage %3.2g USD',tc_track.name,EDS(2).ED);
title(title_str);box on

fprintf('combined damage: %3.2g USD\n',EDS(1).ED+EDS(2).ED);
