function hazard = climada_hazard_merge(hazard1, hazard2)

% hazard1 = hazard_wkn;
% hazard2 = hazard_ext;

if hazard1.reference_year ~= hazard2.reference_year   
    fprintf(['Warning: Reference years are not equal: ' int2str(hazard1.reference_year) ' and ' int2str(hazard2.reference_year) '\n']);
end

if hazard1.peril_ID ~= hazard2.peril_ID   
    fprintf(['Warning: Peril IDs are not equal: ' hazard1.peril_ID ' and ' hazard2.peril_ID '\n']);
end

if hazard1.orig_years ~= hazard2.orig_years   
    fprintf(['Warning: orig years are not equal: ' int2str(hazard1.orig_years) ' and ' int2str(hazard2.orig_years) '\n']);
end

if hazard1.event_count ~= hazard2.event_count   
    fprintf(['Warning: event counts are not equal: ' int2str(hazard1.event_count) ' and ' int2str(hazard2.event_count) '\n']);
end

if size(hazard1.lon,2) ~= size(hazard2.lon,2)
    fprintf('Merge two hazards with different centroids\n');
    
    hazard     = hazard2;
    hazard.lon = [hazard1.lon hazard2.lon];
    hazard.lat = [hazard1.lat hazard2.lat];
    no_cen     = size(hazard.lon,2);
    hazard.centroid_ID = 1:no_cen;
    hazard.intensity = [hazard1.intensity hazard2.intensity];
end



