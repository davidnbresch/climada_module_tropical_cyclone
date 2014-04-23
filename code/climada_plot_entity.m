
function climada_plot_entity(entity)



% set colormap
% miv = min(arrayfun(@(x) (min(x.Value)),assets));
% mav = max(arrayfun(@(x) (max(x.Value)),assets));
% miv = []; mav = [];
beginColor  = [232 232 232 ]/255;
middleColor = [105 105 105 ]/255;
cmap1 = makeColorMap(beginColor, middleColor, 4);
cmap2 = makeColorMap([255 236 139]/255, [255 97 3 ]/255, 6); %[255 153 18]/255 yellow
cmap3 = makeColorMap([255 64 64 ]/255, [176 23 31 ]/255, 2); %[255 153 18]/255 yellow
% cmap3 = makeColorMap([205 150 205 ]/255, [93 71 139 ]/255, 2); %[255 153 18]/255 yellow

cmap  = [cmap1; cmap2; cmap3];


% plot the assets
markersize    = 2;

fig = climada_figuresize(0.5,0.8);
climada_plot_world_borders(0.7);
hold on
xlabel('Longitude')
ylabel('Latitude')
set(gca,'layer','top')

[cbar asset_handle]= plotclr(entity.assets.Longitude, entity.assets.Latitude, entity.assets.Value, 's',markersize, 1,0,[],cmap,1,0); 
set(get(cbar,'ylabel'),'string','USD (exponential)','fontsize',12)  

titlestr = sprintf('%s\n %d assets: %5.0f Bn USD, ', entity.assets.excel_file_name, size(entity.assets.Value,2), sum(entity.assets.Value)*10^-9);
title(titlestr,'fontsize',11)
d = 0;
x_range = [min(entity.assets.Longitude)-d max(entity.assets.Longitude)+d];
y_range = [min(entity.assets.Latitude)-d max(entity.assets.Latitude)+d];

hb = get(gca,'PlotBoxAspectRatio');
hb = hb(1)/hb(2);
if hb/(diff(x_range)/diff(y_range))<1
    dif     = ( diff(x_range)/hb-diff(y_range) )/2;
    y_range = [y_range(1)-dif y_range(2)+dif];
    set(gca,'xlim',x_range,'ylim',y_range)
else
    dif     = ( diff(y_range)*hb-diff(x_range) )/2;
    x_range = [x_range(1)-dif x_range(2)+dif]; 
    set(gca,'xlim',x_range,'ylim',y_range)
end
% axis equal            
            
% for r_i = 1:length(region_str)
%     [cbar asset_handles{r_i}]= plotclr(assets(r_i).Longitude, assets(r_i).Latitude, assets(r_i).Value, 's',markersize, 0,0,[],cmap,1,0); 
%     miv(r_i) = min(assets(r_i).Value);
%     mav(r_i) = max(assets(r_i).Value);
% end
% set([asset_handles{:}],'markersize',3)
% set([asset_handles{:}],'HandleVisibility','off')


% colormap(cmap)
% cbar_assets = colorbar;
% % freezeColors(handles.axes5) %freezeColors(cbar_assets)
% cbar_assets = colorbar('location','north'); %southoutside
% caxis([0 1])
% set(cbar_assets,'xlim',[0 1],'xtick',[])
% cbar_assets = cbfreeze(cbar_assets);
