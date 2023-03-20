
% % plot slice of volume with an overlay of volume fraction
% clear
% filename = ('C:\Users\20220428\OneDrive - Murdoch University\Documents\MATLAB\Sea_Urchin\data\greyscale\3_72_interambulacral_plate_filtered.tif');
% greyscale = tiffreadVolume(filename);
% temp = binary;
% % padding = padarray(temp,[2 0 8],0,'post');
% % temp = padding(1:1150,1:1050,1:1150);
% % imagesc(squeeze(padding(800,:,:))); axis equal
% binary = imbinarize(greyscale,0.5);
% temp = binary(101:1100,26:1025,101:1100);
% % imshow(squeeze(padding(151:1050,500,151:1050)))
imshow(squeeze(padding(151:1050,26:1025,500)))
% temp = padding;

% pixelregionx = [1 1050];
% pixelregiony = [1 1000];
% pixelregionz = [1 1000];
% binarize = 1;
% 
% [greyscale, binary] = readTIF(filename,pixelregionx,pixelregiony,pixelregionz,binarize);

temp = binary(1:1000,1:1000,1:1000);
size_subvolume = 25;
[subvolume] = subvolume_extraction(temp, size_subvolume);
[volume_fraction] = solid_volume_fraction(subvolume);

z_slice = 22;

for i = 1:40

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

% F(y_slice) = getframe;

end


hold on
rectangle('Position',[800 400 25 25],'EdgeColor','w','FaceColor','w') %disordered
rectangle('Position',[280 800 25 25],'EdgeColor','w','FaceColor','w') %ordered
rectangle('Position',[466 550 25 25],'EdgeColor','w','FaceColor','w')
rectangle('Position',[600 690 25 25],'EdgeColor','w','FaceColor','w')
rectangle('Position',[560 511 25 25],'EdgeColor','w','FaceColor','w')
rectangle('Position',[433 110 25 25],'EdgeColor','w','FaceColor','w')


dis1 = mean(volume_frac_map(412:462,490:540),'all') %disordered
dis2 = mean(volume_frac_map(800:825,280:305),'all') % ordered
dis3 = mean(volume_frac_map(550:575,466:491),'all')
dis4 = mean(volume_frac_map(690:715,600:625),'all')
dis5 = mean(volume_frac_map(511:536,560:585),'all')
dis6 = mean(volume_frac_map(110:135,433:458),'all')


plot([90 90], [100 200],'LineWidth',2,'Color','w')



rectangle('Position',[490 412 50 50],'LineWidth',2)


video = VideoWriter('CRT_y.avi'); %create the video object
open(video); %open the file for writing
for ii=1:10 %where N is the number of images
  writeVideo(video,frames(ii)); %write the image to file
end
close(video); %close the file

%X = 845:946
%Y = 759:860
%Z = 797:898
