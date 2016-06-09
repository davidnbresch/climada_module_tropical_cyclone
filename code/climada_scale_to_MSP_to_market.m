function [entity,ELS] = climada_scale_to_MSP_to_market(entity, ELS, MSP_AEL, market_TIV)
% UNDOCUMENTED
%
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('entity'    ,'var'), entity     = []; end
if ~exist('ELS'       ,'var'), ELS        = []; end
if ~exist('MSP_AEL'   ,'var'), MSP_AEL    = []; end
if ~exist('market_TIV','var'), market_TIV = []; end


if isempty(MSP_AEL) && ~isempty(ELS)
    if isfield(entity.assets, 'MSP_Loss')
        MSP_AEL = sum(entity.assets.MSP_Loss);
    end
    factor_AEL_climada_to_MSP = MSP_AEL/ELS.EL;
else
    fprintf('No AEL of MSP is given. Not scaling factor applied to AEL.\n')
    factor_AEL_climada_to_MSP = 1;
end

if isempty(market_TIV)
    factor_GDP_to_mrkt = 4; 
    fprintf('No market TIV is given. Average GDP to market factor 4 applied.\n')
else
    factor_GDP_to_mrkt = market_TIV/sum(entity.assets.Value);
end

% scale up entity from GDP to market
entity.assets.excel_file_name = sprintf('GDP to market scaled with factor %2.3f, %s' , factor_GDP_to_mrkt, entity.assets.excel_file_name);
entity.assets.Value           = factor_GDP_to_mrkt * entity.assets.Value ;
entity.assets.Deductible      = factor_GDP_to_mrkt * entity.assets.Deductible;
entity.assets.Cover           = factor_GDP_to_mrkt * entity.assets.Cover;
entity.assets.Value_2012      = factor_GDP_to_mrkt * entity.assets.Value_2012;
if isfield(entity.assets, 'MSP_Loss')
    entity.assets.MSP_Loss    = factor_GDP_to_mrkt * entity.assets.MSP_Loss;
end

% % calculate new ELS
if ~isempty(ELS)
    ELS.loss        = factor_GDP_to_mrkt * ELS.loss;
    ELS.loss_per_cu = factor_GDP_to_mrkt * ELS.loss_per_cu;
    ELS.EL_per_cu   = factor_GDP_to_mrkt * factor_AEL_climada_to_MSP * ELS.EL_per_cu;
    ELS.Value       = factor_GDP_to_mrkt * ELS.Value;
    ELS.comment     = sprintf('Scaled up ELS with factor %2.2f to market and %2.2f to MSP, %s',...
                               factor_GDP_to_mrkt, factor_AEL_climada_to_MSP, ELS.comment);
    ELS.annotation_name = sprintf('Scaled up ELS with factor %2.2f to market and %2.2f to MSP, %s',...
                                   factor_GDP_to_mrkt, factor_AEL_climada_to_MSP, ELS.annotation_name);
    ELS.EL          = factor_GDP_to_mrkt * factor_AEL_climada_to_MSP * ELS.EL;                           
else
    hazard = '';
    fprintf('Hazard to be loaded\n')
    annotation_name = sprintf('Scaled up entity with factor %2.2f', factor_GDP_to_mrkt);
    ELS = climada_ELS_calc(entity, hazard, annotation_name);
end




% ELS_scaled = ELS;
% for i = 1:length(ELS)
%     ELS_scaled(i).Value = factor_climada_to_MSP * factor_GDP_to_mrkt * ELS(i).Value;
%     %ELS_scaled(i).loss  = factor_climada_to_MSP * factor_GDP_to_mrkt * ELS(i).loss;
%     ELS_scaled(i).loss  = factor_GDP_to_mrkt * ELS(i).loss;
%     ELS_scaled(i).EL    = factor_climada_to_MSP * factor_GDP_to_mrkt * ELS(i).EL;
% end


% 
% % Loss_Market_Portfolio = 17*10^9;
% 
% % AEL_MSP               = 4732256617.75385;
% % AEL_MSP               = 16789587264.449;
% % TIV_Market_Portfolio  = 27360*10^9;
% TIV_Market_Portfolio    = 30810*10^9;
% 
% % 1) climada to MSP
% % take AEL from TC market portfolio 2012
% AEL_MSP               = 0.0544917132427529/100 * ELS(1).Value;
% factor_climada_to_MSP = AEL_MSP/ELS(1).EL;
% % 2) GDP prtf to market portfolio
% factor_GDP_to_mrkt    = TIV_Market_Portfolio/ELS(1).Value;
% 
% ELS_scaled = ELS;
% for i = 1:length(ELS)
%     ELS_scaled(i).Value = factor_climada_to_MSP * factor_GDP_to_mrkt * ELS(i).Value;
%     %ELS_scaled(i).loss  = factor_climada_to_MSP * factor_GDP_to_mrkt * ELS(i).loss;
%     ELS_scaled(i).loss  = factor_GDP_to_mrkt * ELS(i).loss;
%     ELS_scaled(i).EL    = factor_climada_to_MSP * factor_GDP_to_mrkt * ELS(i).EL;
% end
% % AEL in percentage
% % ELS_scaled(i).EL/ELS_scaled(i).Value*1000
% 
% 
% % % % compare with Multisnap TC market portfolio LFC
% R_MSP = [1 2 5 10 20 50 100 200 500 1000 2000 5000 10000];
% L_MSP_mrkt = [2.070940E+09	6.331243E+09	1.779480E+10	3.276317E+10	5.170143E+10	8.657718E+10	1.194006E+11	1.572744E+11	2.155808E+11	2.602951E+11	3.239819E+11	3.833480E+11	4.434115E+11];
% L_MSP = [830000000,2400000000.00000,5800000000.00000,9600000000.00000,14000000000.0000,23000000000.0000,30000000000.0000,38000000000.0000,49000000000.0000,58000000000.0000,67000000000.0000,78000000000.0000,91000000000.0000];
% 
% 
% %% market view
% % % % compare with Multisnap TC market portfolio LFC
% for i = 1:length(ELS_scaled)
%     ELS_scaled_(i) = climada_ELS_stats(ELS_scaled(i),[],R_MSP,0);
% end
% % b = ELS_scaled_(1).loss_fit';
% % b = ELS_scaled_(1).EL;
% 
% color_     = [255 215 0 ;...   %today
%               255 127 0 ;...   %eco 
%                238 64 0 ;...  %clim
%               205 0 0 ;...   %total risk
%               120 120 120]/256; %dotted line]/255;
% % color_(1:4,:) = brighten(color_(1:4,:),0.3);     
% digits = 9;
% 
% fig = climada_figuresize(0.5,0.8);
% hold on
% for i = length(ELS_scaled_):-1:1
%     plot(ELS_scaled_(i).R_fit, ELS_scaled_(i).loss_fit'*10^-digits,'o-','color',color_(i,:),'linewidth',2,'markersize',5)
% end
% plot(R_MSP,L_MSP_mrkt*10^-digits,'d--k','linewidth',2,'markersize',3)
% plot(R_MSP,L_MSP * factor_GDP_to_mrkt * 10^-digits,'d-k','linewidth',2,'markersize',3)
% xlabel('Return period (years)','fontsize',12)
% ylabel(sprintf('Loss 10^{%d} USD',digits),'fontsize',12)
% xlim([0 550])
% ytick_ = get(gca,'ytick');
% ytick_ = ytick_(1:2:end);
% set(gca,'ytick',ytick_,'ygrid','on','fontsize',12)
% L = legend('climada cc 2030','climada eco 2030','climada 2012','MSP market portfolio','MSP GDP upscaled','location','southeast');
% % legend('Multisnap, AEL: 4.7E+09',sprintf('climada, AEL: %2.2g', ELS(1).EL),'location','southeast')
% legend('boxoff')
% set(L,'fontsize',12)
% title('USA: market view','fontsize',16)
% print(fig,'-dpdf',[climada_global.data_dir filesep 'results' filesep 'TC_global' filesep 'LFC_2012_eco_cc_2030_USA.pdf'])  
% 
% 
% 
% % %% GDP view.
% % for i = 1:length(ELS)
% %     ELS_(i) = climada_ELS_stats(ELS(i),[],R_MSP,0);
% % end
% % 
% % fig = climada_figuresize(0.5,0.8);
% % hold on
% % for i = length(ELS_):-1:1
% %     plot(ELS_(i).R_fit, ELS_(i).loss_fit'*10^-digits,'o-','color',color_(i,:),'linewidth',2,'markersize',5)
% % end
% % % plot(R_MSP,L_MSP_mrkt*10^-digits,'d--k','linewidth',2,'markersize',3)
% % plot(R_MSP,L_MSP*10^-digits,'d-k','linewidth',2,'markersize',3)
% % xlabel('Return period (years)','fontsize',12)
% % ylabel(sprintf('Loss 10^{%d} USD',digits),'fontsize',12)
% % xlim([0 550])
% % ytick_ = get(gca,'ytick');
% % ytick_ = ytick_(1:2:end);
% % set(gca,'ytick',ytick_,'ygrid','on','fontsize',12)
% % L = legend('climada cc 2030','climada eco 2030','climada 2012','MSP GDP','location','southeast');
% % % legend('Multisnap, AEL: 4.7E+09',sprintf('climada, AEL: %2.2g', ELS(1).EL),'location','southeast')
% % legend('boxoff')
% % set(L,'fontsize',12)
% % title('USA: GDP view','fontsize',16)
% % print(fig,'-dpdf',[climada_global.data_dir filesep 'results' filesep 'TC_global' filesep 'LFC_2012_eco_cc_2030_USA_GDP_view.pdf'])  
% % 
% % %%
% % 
% % 
% % % Percentage_Of_Value_Flag = 0;
% % % fig = climada_ELS_LFC(ELS_scaled(1), ELS_scaled(2:3), Percentage_Of_Value_Flag);
% % % %title off
% % % set(get(gca,'title'),'string','')
% % % xlim([0 550])
% % % print(fig,'-dpdf',[climada_global.data_dir filesep 'results' filesep 'TC_global' filesep 'LFC_2012_eco_cc_2030_US.pdf'])  
% % 
% % % waterfall graph
% % fig = climada_waterfall_graph(ELS_scaled(1),ELS_scaled(2), ELS_scaled(3), 'AEL');
% % print(fig,'-dpdf',[climada_global.data_dir filesep 'results' filesep 'TC_global' filesep 'waterfall_AEL_US.pdf'])   
% % close
% % 
% % fig = climada_waterfall_graph(ELS_scaled(1),ELS_scaled(2), ELS_scaled(3), 100);
% % print(fig,'-dpdf',[climada_global.data_dir filesep 'results' filesep 'TC_global' filesep 'waterfall_100year_US.pdf'])   
% % close
% % 
% % fig = climada_waterfall_graph(ELS_scaled(1),ELS_scaled(2), ELS_scaled(3), 250);
% % print(fig,'-dpdf',[climada_global.data_dir filesep 'results' filesep 'TC_global' filesep 'waterfall_250year_US.pdf'])   
% % close
% % 
% % fig = climada_waterfall_graph(ELS_scaled(1),ELS_scaled(2), ELS_scaled(3), 500);
% % print(fig,'-dpdf',[climada_global.data_dir filesep 'results' filesep 'TC_global' filesep 'waterfall_500year_US.pdf'])
% % close
% % 
% % 
% % 
% % 
% % %% FINAL calibration (A = 0.85, A = 0.9)
% % % % wind decay
% % ELS(1) = climada_ELS_calc(entity, hazard_unwkn,'2012 US unweakened');
% % ELS(2) = climada_ELS_calc(entity, hazard_wkn1,'2012 US wkn ori');
% % 
% % hazard_wkn_v2 = climada_hazard_distance_to_coast_USA(hazard_unwkn, centroids, tc_track, 1);
% % 
% % 
% % % % Multisnap US
% % R_MSP = [1 2 5 10 20 50 100 200 500 1000 2000 5000 10000];
% % L_MSP = [830000000,2400000000.00000,5800000000.00000,9600000000.00000,14000000000.0000,23000000000.0000,30000000000.0000,38000000000.0000,49000000000.0000,58000000000.0000,67000000000.0000,78000000000.0000,91000000000.0000];
% % % 
% % % % strong wind decay (version 2)
% % % hazard_wkn = climada_hazard_distance_to_coast(hazard, centroids, tc_track);
% % % ELS(1)     = climada_ELS_calc(entity, hazard_wkn,'2012 US strong wind decay over land');
% % figure
% % ELS_2      = climada_ELS_stats(ELS(1),[],R_MSP,1);
% % % % b          = ELS_2.loss_fit';
% % % % ELS.EL
% % % 
% % climada_figuresize(0.5,0.8)
% % plot(R_MSP,L_MSP,'.-k')
% % hold on
% % plot(ELS_2.R_fit, ELS_2.loss_fit','*-')
% % xlabel('Return period (years)')
% % ylabel('Loss')
% % xlim([0 550])
% % legend('Multisnap, AEL: 4.7E+09',sprintf('climada, AEL: %2.2g', ELS(1).EL),'location','southeast')
% % legend('boxoff')
% % 
% % hazard_wkn_v2 = climada_hazard_distance_to_coast_USA(hazard_unwkn, centroids, tc_track, 1);
% % ELS(2) = climada_ELS_calc(entity, hazard_wkn_v2,'2012 US wkn v2');
% % ELS_2  = climada_ELS_stats(ELS(2),[],[],0);
% % plot(ELS_2.R_fit, ELS_2.loss_fit','*--r')
% % b      = ELS_2.loss_fit';
% % ELS_2.EL
% % 
