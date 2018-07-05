function [screw, basin] = climada_hazard_climate_set_screw(reference, hazard_type, basin)
%%%
%   
% MODULE:
%   tropical_cyclone
% NAME:
%   climada_hazard_climate_set_screw
% PURPOSE:
%    Set struct screw for TC based on scientific reference key and ocean basin for climada_hazard_climate_screw
% CALLING SEQUENCE:
%   screw = climada_hazard_climate_set_screw('Knutson2015', 'TC','NA');
%   hazard = climada_hazard_load(hazard_filename);
%   hazard = climada_hazard_climate_screw(hazard, hazard_set_file, reference_year, screw);
% EXAMPLE:
%   screw = climada_hazard_climate_set_screw('IPCC_SREX',[],'GL');
% INPUTS:
%   reference:   CHAR/STRING linked to specific scientific reference. So far implemented:
%                        'IPCC_SREX' :: IPCC SREX, Chapter 3,
%                            https://www.ipcc.ch/pdf/special-reports/srex/SREX-Chap3_FINAL.pdf (default)
%                        'Knutson2015' :: Knutson et al., 2015. Global projections of intense tropical cyclone activity for the late twenty-first century from dynamical downscaling of CMIP5/RCP4.5 scenarios
%                            (CMIP5 RCP4.5 2081?2100 (i.e., late twenty-first century) and 2001?20 (i.e., present day))
% OPTIONAL INPUT PARAMETERS:
%   hazard_type: CHAR/STRING defining hazard type:
%                         'TC' :: Tropical Cyclone (default)
%                         'TR' :: Tropical Cyclone Rainfall
%                Note: base new TS hazard (storm surge) on srewed TC hazard. No extra TS screw.
%   basin:       CHAR/STRING defining specific ocean basin (at the moment only for Knutson2015):
%                         'GL' :: Global (default)
%                         'EP' :: North East Pacific Ocean
%                         'NA' :: North Atlantic Ocean
%                         'NI' :: North Indian Ocean
%                         'SA' :: South Atlantic Ocean
%                         'SI' :: South Indian Ocean
%                         'SP' :: South Pacific Ocean
%                         'WP' :: North West Pacific Ocean
% OUTPUTS:
%   screw:      defines the climate change scenario. A 1xN structure with fields:
%                   .hazard_fld     defines the hazard field to be changed
%                   .change         extent of the change at time horizon
%                   .year           time horizon of change
%                   .hazard_crit    hazard field to which criteria apply
%                   .criteria       criteria for events/locations to change
%                   .bsxfun_op      operation of change (e.g. @times,@plus) (function handle)
%               specifying N transformations to the original hazard set.
%
%  basin:       same as input parameter
%
% MODIFICATION HISTORY:
% Samuel Eberenz, eberenz@posteo.eu, 20171114, initial
% Samuel Eberenz, eberenz@posteo.eu, 20180705, init at tropical_cyclone module
%-

if ~exist('reference', 'var')   , reference = 'IPCC_SREX'; end
if ~exist('hazard_type', 'var') , hazard_type = 'TC'; end
if ~exist('basin', 'var')       , basin = 'GL'; end
i = 1; 

if ~iscell(basin)
    basin{1} = basin;
end

for i_basin = 1:length(basin)
    switch reference % SWITCH between input parameter reference, then check for basin and hazard_type using IF
        case 'IPCC_SREX'
        % Seneviratne et al., 2012
        % IPCC SREX, Chapter 3, https://www.ipcc.ch/pdf/special-reports/srex/SREX-Chap3_FINAL.pdf
            %enhance frequency of category 4 and 5 storms
            if isequal(hazard_type,'TC') && isequal(basin{i_basin},'GL')
                screw(i).hazard_fld         = 'frequency';
                screw(i).change             = 1.8;
                screw(i).year               = 2100;
                screw(i).hazard_crit        = 'category';
                screw(i).criteria           = [4 5];
                screw(i).bsxfun_op          = @times;
                i=i+1;

                %decrease overall frequency of all storms (0 up to category 5)
                screw(i).hazard_fld         = 'frequency';
                screw(i).change             = 0.72;
                screw(i).year               = 2100;
                screw(i).hazard_crit        = 'category';
                screw(i).criteria           = [0 1 2 3];
                screw(i).bsxfun_op          = @times;
                i=i+1;

                %increase global mean maximum wind speed (intensity)
                screw(i).hazard_fld         = 'intensity';
                screw(i).change             = 1.11;
                screw(i).year               = 2100;
                screw(i).hazard_crit        = 'category';
                screw(i).criteria           = [0 1 2 3 4 5];
                screw(i).bsxfun_op          = @times;
                i=i+1;
            else
                cprintf([1 0 0], sprintf('ERROR: No %s screw defined for given input parameters: %s, %s \n',reference, hazard_type, basin{i_basin}));     
            end
        case 'Knutson2015'
        % Knutson et al. 2015
        % Global projections of intense tropical cyclone activity for the late twenty-first century from dynamical downscaling of CMIP5/RCP4.5 scenarios
        % Global mean results only. But big differences between ocean basins
        % --> implement differentiation of ocean basins within hazard-struct and screw!
        if isequal(hazard_type,'TC') && isequal(basin{i_basin},'GL')
            %increase global frequency  (u > 58m/s ~ category 4 and 5) {!basins vary strongly}
            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 1.28; 
            screw(i).year               = 2100;
            screw(i).hazard_crit        = 'category';
            screw(i).criteria           = 4;
            screw(i).bsxfun_op          = @times;
            i=i+1;

            %increase global frequency  (u > 65m/s ~ category 5) {not southwest pacific}
            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 1.59;
            screw(i).year               = 2100;
            screw(i).hazard_crit        = 'category';
            screw(i).criteria           = 5;
            screw(i).bsxfun_op          = @times;
            i=i+1;

            %decrease overall frequency of all storms (0 up to category 5)
            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 0.83; % -17%
            screw(i).year               = 2100;
            screw(i).hazard_crit        = 'category';
            screw(i).criteria           = [0 1 2 3];
            screw(i).bsxfun_op          = @times;
            i=i+1;

            %increase global mean maximum wind speed (intensity) for Tropical Storms (maxwnd_ts)
            screw(i).hazard_fld         = 'intensity';
            screw(i).change             = 1.036;
            screw(i).year               = 2100;
            screw(i).hazard_crit        = 'category';
            screw(i).criteria           = 0;
            screw(i).bsxfun_op          = @times;   
            i=i+1;

            %increase global mean maximum wind speed (intensity) for Tropical Cyclones (maxwnd_hur)
            screw(i).hazard_fld         = 'intensity';
            screw(i).change             = 1.041;
            screw(i).year               = 2100;
            screw(i).hazard_crit        = 'category';
            screw(i).criteria           = [1 2 3 4 5];
            screw(i).bsxfun_op          = @times;    
            i=i+1;


        elseif isequal(hazard_type,'TC') && isequal(basin{i_basin},'NA') % Knutson 2015
            %increase North Atlantic frequency  (u > 58m/s ~ category 4 and 5) {!basins vary strongly}
    %         screw(i).hazard_fld         = 'frequency';
    %         screw(i).change             = 1.42; % not significant!
    %         screw(i).year               = 2100;
    %         screw(i).hazard_crit        = {'category','basin'};
    %         screw(i).criteria           = {4 'NA'};
    %         screw(i).bsxfun_op          = @times;
    %         i=i+1;
    %
    %         %increase  frequency  (u > 65m/s ~ category 5) {not southwest pacific}
    %         screw(i).hazard_fld         = 'frequency';
    %         screw(i).change             = 2.25; % additional 125% % not significant!
    %         screw(i).year               = 2100;
    %         screw(i).hazard_crit        = {'category','basin'};
    %         screw(i).criteria           = {5 'NA'};
    %         screw(i).bsxfun_op          = @times;
    %         i=i+1;
    %         
    %         %decrease frequency of weak storms (0)
    %         screw(i).hazard_fld         = 'frequency';
    %         screw(i).change             = 0.91; % -9.4% % not significant!
    %         screw(i).year               = 2100;
    %         screw(i).hazard_crit        = {'category','basin'};
    %         screw(i).criteria           = {0 'NA'};
    %         screw(i).bsxfun_op          = @times;
    %         i=i+1;
    %
    %         %decrease frequency of weak TC (1-3)
    %         screw(i).hazard_fld         = 'frequency';
    %         screw(i).change             = 0.82; % -17.5% % not significant!
    %         screw(i).year               = 2100;
    %         screw(i).hazard_crit        = {'category','basin'};
    %         screw(i).criteria           = {[1 2 3] 'NA'};
    %         screw(i).bsxfun_op          = @times;
    %         i=i+1;
    %
    %         %increase mean maximum wind speed (intensity) for Tropical Storms (maxwnd_ts)
    %         screw(i).hazard_fld         = 'intensity';
    %         screw(i).change             = 1.004; % not significant!
    %         screw(i).year               = 2100;
    %         screw(i).hazard_crit        = {'category','basin'};
    %         screw(i).criteria           = {0 'NA'};
    %         screw(i).bsxfun_op          = @times;   
    %         i=i+1;

            %increase mean maximum wind speed (intensity) for Tropical Cyclones (maxwnd_hur)
            screw(i).hazard_fld         = 'intensity';
            screw(i).change             = 1.045;
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {[1 2 3 4 5] basin{i_basin}};
            screw(i).bsxfun_op          = @times;    
            i=i+1;

         elseif isequal(hazard_type,'TC') && isequal(basin{i_basin},'EP') % Knutson 2015
            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 1.163;
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {0 basin{i_basin}};
            screw(i).bsxfun_op          = @times;
            i=i+1;

            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 1.193; 
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {[1 2] basin{i_basin}};
            screw(i).bsxfun_op          = @times;
            i=i+1;

            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 1.837; 
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {3 basin{i_basin}};
            screw(i).bsxfun_op          = @times;
            i=i+1;

            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 3.375; 
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {[4 5] basin{i_basin}};
            screw(i).bsxfun_op          = @times;
            i=i+1;        

            % maxwnd_ts
            screw(i).hazard_fld         = 'intensity';
            screw(i).change             = 1.082; 
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {0 basin{i_basin}};
            screw(i).bsxfun_op          = @times;   
            i=i+1;

            % maxwnd_hur
            screw(i).hazard_fld         = 'intensity';
            screw(i).change             = 1.078;
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {[1 2 3 4 5] basin{i_basin}};
            screw(i).bsxfun_op          = @times;    
            i=i+1;

         elseif isequal(hazard_type,'TC') && isequal(basin{i_basin},'WP') % Knutson 2015
            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 1-.345;
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {0 basin{i_basin}};
            screw(i).bsxfun_op          = @times;
            i=i+1;

            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 1-.316; 
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {[1 2] basin{i_basin}};
            screw(i).bsxfun_op          = @times;
            i=i+1;

            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 1-.169; 
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {[3 4 5] basin{i_basin}};
            screw(i).bsxfun_op          = @times;
            i=i+1;        

            % maxwnd_ts
            screw(i).hazard_fld         = 'intensity';
            screw(i).change             = 1.074; 
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {0 basin{i_basin}};
            screw(i).bsxfun_op          = @times;   
            i=i+1;

            % maxwnd_hur
            screw(i).hazard_fld         = 'intensity';
            screw(i).change             = 1.055;
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {[1 2 3 4 5] basin{i_basin}};
            screw(i).bsxfun_op          = @times;    
            i=i+1;

        elseif isequal(hazard_type,'TC') && isequal(basin{i_basin},'NI') % Knutson 2015

            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 1.256; 
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {[1 2 3 4 5] basin{i_basin}};
            screw(i).bsxfun_op          = @times;
            i=i+1;

        elseif isequal(hazard_type,'TC') && isequal(basin{i_basin},'SI') % Knutson 2015
            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 1-.261;
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {0 basin{i_basin}};
            screw(i).bsxfun_op          = @times;
            i=i+1;

            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 1-.284; 
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {[1 2 3 4 5] basin{i_basin}};
            screw(i).bsxfun_op          = @times;
            i=i+1;

            % maxwnd_hur
            screw(i).hazard_fld         = 'intensity';
            screw(i).change             = 1.033;
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {[1 2 3 4 5] basin{i_basin}};
            screw(i).bsxfun_op          = @times;    
            i=i+1;

         elseif isequal(hazard_type,'TC') && isequal(basin{i_basin},'SP') % Knutson 2015
            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 1-.366;
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {0 basin{i_basin}};
            screw(i).bsxfun_op          = @times;
            i=i+1;

            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 1-.406; 
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {[1 2] basin{i_basin}};
            screw(i).bsxfun_op          = @times;
            i=i+1;

            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 1-.506; 
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {3 basin{i_basin}};
            screw(i).bsxfun_op          = @times;
            i=i+1;  

            screw(i).hazard_fld         = 'frequency';
            screw(i).change             = 1-.583; 
            screw(i).year               = 2100;
            screw(i).hazard_crit        = {'category','basin'};
            screw(i).criteria           = {[4 5] basin{i_basin}};
            screw(i).bsxfun_op          = @times;
            i=i+1;   

        elseif isequal(hazard_type,'TR') && isequal(basin{i_basin},'GL')
        % Knutson et al. 2015
        % Global projections of intense tropical cyclone activity for the late twenty-first century from dynamical downscaling of CMIP5/RCP4.5 scenarios
        % Global mean results only. But big differences between ocean basins
        % --> implement differentiation of ocean basins within hazard-struct and screw!

            %increase global rain rate
            screw(i).hazard_fld         = 'intensity';
            screw(i).change             = 1.143;
            screw(i).year               = 2100;
            screw(i).hazard_crit        = 'category';
            screw(i).criteria           = 0;
            screw(i).bsxfun_op          = @times;
            i=i+1;

            screw(i).hazard_fld         = 'intensity';
            screw(i).change             = 1.134;
            screw(i).year               = 2100;
            screw(i).hazard_crit        = 'category';
            screw(i).criteria           = [1 2];
            screw(i).bsxfun_op          = @times;
            i=i+1;

            screw(i).hazard_fld         = 'intensity';
            screw(i).change             = 1.088; 
            screw(i).year               = 2100;
            screw(i).hazard_crit        = 'category';
            screw(i).criteria           = 3;
            screw(i).bsxfun_op          = @times;
            i=i+1;

            screw(i).hazard_fld         = 'intensity';
            screw(i).change             = 1.077;
            screw(i).year               = 2100;
            screw(i).hazard_crit        = 'category';
            screw(i).criteria           = [4 5];
            screw(i).bsxfun_op          = @times;  
            i=i+1;
        else
            cprintf([1 0 0], sprintf('ERROR: No %s screw defined for given input parameters: %s, %s \n',reference, hazard_type, basin{i_basin}));     
        end
    otherwise
        cprintf([1 0 0], sprintf('ERROR: No %s screw defined for given input parameters: %s, %s \n',reference, hazard_type, basin{i_basin}));     
    end
end
for j=1:i-1, screw(j).reference = reference; end
end
