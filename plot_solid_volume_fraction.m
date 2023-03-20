
% % plot slice of volume with an overlay of volume fraction
clear
filename = ('C:\Users\20220428\OneDrive - Murdoch University\Documents\MATLAB\Sea_Urchin\data\greyscale\3_72_interambulacral_plate_filtered.tif');
greyscale = tiffreadVolume(filename);
binary = imbinarize(greyscale,0.46);
temp = binary(1:1000,1:1000,1:1000);
size_subvolume = 25;
[subvolume] = subvolume_extraction(temp, size_subvolume);
[volume_fraction] = solid_volume_fraction(subvolume);

% x_slice = 20;
% z_slice = 22;
y_slice = 21;

figure(1);
ax1 = axes;
x_z_vol_frac = squeeze(volume_fraction(1:size(volume_fraction,1),y_slice,1:size(volume_fraction,3))); %example for a slice of y
% x_y_vol_frac = squeeze(volume_fraction(1:size(x,1),1:size(x,1),1)); %example for a slice of z
volume_frac_map = imresize(x_z_vol_frac,size_subvolume,"bicubic"); %scales volume fraction matrix up by 
% 50 to fit image volume and uses bicubic interpolation by default
% "nearest"
% volume_frac_map = imresize(x_y_vol_frac,size(x,2),"bicubic");
imagesc(volume_frac_map,'AlphaData',0.6); axis equal tight; %displays volume fraction for that slice of the volume
set(ax1,'YDir','normal')
colormap(ax1,"jet")
caxis([0 0.6])
box off

ax2 = axes;
imagesc(squeeze(temp(:,526,:)),'AlphaData',0.3); axis equal tight; %displays slice through y
% imagesc(squeeze(binary(:,:,1)),'AlphaData',0.5) %displays slice through z
set(ax2,'YDir','normal')
box off

linkaxes([ax1,ax2])

ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

colormap(ax2,"gray")

colorbar(ax1,'Position',[0.85 0.15 0.025 0.7])
