
function ELS = climada_ELS_scaleup_factor(ELS, factor_)

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('ELS'    ,'var'), ELS     = []; end
if ~exist('factor_','var'), factor_ = []; end

if isempty(ELS)
    ELS = climada_ELS_load;
end

if isempty(factor_)
    fprintf('No scaleup factor given. Unable to proceed.\n')
    return
end

% % scale up ELS with given factor
if ~isempty(ELS)
    ELS.loss        = factor_ * ELS.loss;
    ELS.loss_per_cu = factor_ * ELS.loss_per_cu;
    ELS.EL_per_cu   = factor_ * ELS.EL_per_cu;
    ELS.Value       = factor_ * ELS.Value;
    ELS.comment     = sprintf('Scaled up ELS with factor %2.2f, %s',...
                               factor_, ELS.comment);
    ELS.annotation_name = sprintf('Scaled up ELS with factor %2.2f, %s',...
                                   factor_, ELS.annotation_name);
    ELS.EL          = factor_ * ELS.EL;                           
end




