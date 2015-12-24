function hazard=climada_ts_hazard_set(hazard_TC,hazard_set_file,save_bathymetry_flag,check_plots)
% climada storm surge TS hazard event set
% NAME:
%   climada_ts_hazard_set
% PURPOSE:
%   create a storm surge (TS) hazard event, based on an existing
%   tropical cyclone (TC) hazard event set
%
%   just a caller for tc_surge_hazard_create, as there are routines for
%   - climada_tc_hazard_set
%   - climada_tr_hazard_set
%
%   the special thing is that TS only needs an existing TC hazard set, it
%   does not start from tc tracks again (as TC and TR do)
%
%   see climada_ts_hazard_set
%
% CALLING SEQUENCE:
%   hazard=climada_ts_hazard_set(hazard,hazard_set_file,save_bathymetry_flag,check_plots)
% EXAMPLE:
%   hazard=climada_ts_hazard_set
% INPUTS:
%   hazard_TC: an already existing tropical cyclone (TC) hazard event set (a
%       TC hazard structure, see climada_tc_hazard_set)
%       > prompted for if not given (for .mat file containing a TC hazard event set)
%   hazard_set_file: the name of the newly created storm surge (TS) hazard
%       event set
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   save_bathymetry_flag: if =1, save bathymetry in a .mat file (speeds up
%       subsequent calls), =0 do not save bathymetry (default)
%   check_plots: =1, do show check plots, =0: no plots (default)
% OUTPUTS:
%   hazard: a hazard event set, see climada doc
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20140922 (over the Atlantic, LX016)
% David N. Bresch, david.bresch@gmail.com, 20141026, save_bathymetry_flag
%-
hazard=[]; % init output
if ~exist('hazard_TC','var'),hazard=[];end
if ~exist('hazard_set_file','var'),hazard_set_file=[];end
if ~exist('save_bathymetry_flag','var'),save_bathymetry_flag=0;end
if ~exist('check_plots','var'),check_plots=0;end

hazard=tc_surge_hazard_create(hazard_TC,hazard_set_file,save_bathymetry_flag,check_plots);

return


