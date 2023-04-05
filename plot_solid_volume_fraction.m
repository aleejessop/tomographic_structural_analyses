
%% load data, binarize and calculate solid volume fraction
clear
% filename = 'C:\Users\20220428\OneDrive - Murdoch University\Documents\Dragonfly\1_69_inner_interambulacral_plate_newly_derived.tiff';;
filename = 'C:\Users\20220428\OneDrive - Murdoch University\Documents\Dragonfly\732nm_inner_interambulacral_plate_newly_derived.tiff';
% filename = ('C:\Users\20220428\OneDrive - Murdoch University\Documents\Sea Urchin\Figures\Figure3\Diamond_sample1.tiff');
greyscale = tiffreadVolume(filename);
binary = imbinarize(greyscale,0.47);
%padding = padarray(temp,[2 0 8],0,'post');
% temp = binary(1:1400,1:1750,1:1700);
temp = binary(1:1100,1:1100,1:1050);
size_subvolume = 50;
[subvolume] = subvolume_extraction(temp, size_subvolume);
[volume_fraction] = solid_volume_fraction(subvolume);

%% plot solid volume fraction overlaid onto binary image
%z_slice = 22;
y_slice = 21;

figure(1);
ax1 = axes;
x_z_vol_frac = squeeze(volume_fraction(1:size(volume_fraction,1),y_slice,1:size(volume_fraction,3))); %example for a slice of y
% x_y_vol_frac = squeeze(volume_fraction(1:size(x,1),1:size(x,1),1)); %example for a slice of z
volume_frac_map = imresize(x_z_vol_frac,size_subvolume,"nearest"); %scales volume fraction matrix up by 
% 50 to fit image volume and uses bicubic interpolation by default
% "nearest"
volume_frac_map_norm = rescale(volume_frac_map);
% volume_frac_map = imresize(x_y_vol_frac,size(x,2),"bicubic");
imagesc(volume_frac_map_norm,'AlphaData',0.6); axis equal tight; %displays volume fraction for that slice of the volume
set(ax1,'YDir','normal')
colormap(ax1,"jet")
caxis([0 0.6])
box off

ax2 = axes;
imagesc(squeeze(temp(:,525,:)),'AlphaData',0.3); axis equal tight; %displays slice through y
% imagesc(squeeze(binary(:,:,1)),'AlphaData',0.5) %displays slice through z
set(ax2,'YDir','normal')
box off

linkaxes([ax1,ax2])

ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

colormap(ax2,"gray")

colorbar(ax1,'Position',[0.85 0.15 0.025 0.7])

%% histogram through sample
bins = 0:0.1:1;

sample_through_z1 = discretize(squeeze(volume_frac_map_norm(355,1:1125)),bins);
sample_through_z2 = discretize(squeeze(volume_frac_map_norm(555,1:1125)),bins);
sample_through_z3 = discretize(squeeze(volume_frac_map_norm(755,1:1125)),bins);

figure;
z_axis = 1:1125;
bar(z_axis,bins(sample_through_z1))
box off
hold on




