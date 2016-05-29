function [rainrate] = climada_RCLIPER(fmaxwind_kn,inreach,Radius_km)                                  
% calculate rainrate based on RCLIPER in mm/h
% NAME:
%   climada_RCLIPER
% PURPOSE:
%   given the windspeed (kn) at a specific node calculate the rainrate at 
%   all centroids according to RCLIPER (symmetric rainfield)
%
%   usually called from: climada_tr_rainfield (see there)
% CALLING SEQUENCE:
%   climada_RCLIPER(fmaxwind_kn, inreach, Radius_km)
% EXAMPLE:
%   rainrate=climada_RCLIPER(tc_track.MaxSustainedWind(i),inreach,fRadius_km);      
% INPUTS:
%   fmaxwind_kn: maximum sustained wind at specific node (array) 
%   inreach:     logical vector of centroids length, containing 1 if centroid
%                is within (3, see climada_tr_rainfield) deg of node, otherwise 0 
%   Radius_km:   vector of centroids length, containing distance to node
%                for every centroid
% OPTIONAL INPUT PARAMETERS:
%   none
% OUTPUTS:
%   res.G: the rainrate [mm/h] at all centroids
%       the single-character variables refer to the Pioneer offering circular
%       that's why we kept these short names (so one can copy the OC for
%       documentation)
%   res.lat: the latitude of the centroids
%   res.lon: the longitude of the centroids
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Lea Mueller, 201106038
% david.bresch@gmail.com, 20140804, GIT update
% david.bresch@gmail.com, 20160529, header edited
%-

rainrate = zeros(1,length(inreach));

% Calculate CLIPER rain field
% ---------------------------

% Define Coefficients (CLIPER NHC bias adjusted (Tuleya, 2007))
a1 =  -1.1 ; %inch per day
a2 =  -1.6 ; %inch per day
a3 =  64   ; %km
a4 = 150   ; %km

b1 =   3.96; %inch per day
b2 =   4.8 ; %inch per day
b3 = -13   ; %km
b4 = -16   ; %km

% Calculate the normalized windspeed U in knots
U_norm_kn = (1 + (fmaxwind_kn-35)/33); 

% Calculate parameters dependent on normalized windspeed U
T0 = a1 + b1*U_norm_kn;
Tm = a2 + b2*U_norm_kn;
rm = a3 + b3*U_norm_kn;
r0 = a4 + b4*U_norm_kn;

%Logical vector (inreach & specific distance for R-CLIPER Input)
i  = Radius_km <= rm & inreach;
ii = Radius_km >  rm & inreach;

% Calculate R-CLIPER Symmetric rain rate in mm /h
rainrate(i)    = (T0+(Tm-T0)*(Radius_km(i)/rm))   /24 * 25.4 ;
rainrate(ii)   = (Tm*exp(-(Radius_km(ii)-rm)/r0)) /24 * 25.4 ;

rainrate(isnan(rainrate)) = 0;
rainrate(rainrate < 0)    = 0; 
    
end % climada_RCLIPER