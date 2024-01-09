%% plot solid volume fraction overlaid onto binary image for every slice

y_slice = 1:42;

for i = 1:numel(y_slice)
    x_z_vol_frac(:,:,i) = squeeze(volume_fraction(1:size(volume_fraction,1),y_slice(i),1:size(volume_fraction,3)));
    volume_frac_map(:,:,i) = imresize(x_z_vol_frac(:,:,i),size_subvolume,"bicubic");
    volume_frac_map_norm(:,:,i) = rescale(volume_frac_map(:,:,i));
end


figure(1);
ax1 = axes;

imagesc(volume_frac_map_norm(:,:,21),'AlphaData',0.6); axis equal tight; %displays volume fraction for that slice of the volume
set(ax1,'YDir','normal')
colormap(ax1,"jet")
caxis([0 0.6])
box off

ax2 = axes;
imagesc(squeeze(diamond(:,525,:)),'AlphaData',0.3); axis equal tight; %displays slice through y
% imagesc(squeeze(binary(:,:,1)),'AlphaData',0.5) %displays slice through z
set(ax2,'YDir','normal')
box off

linkaxes([ax1,ax2])

ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

colormap(ax2,"gray")

colorbar(ax1,'Position',[0.85 0.15 0.025 0.7])

histogram(volume_frac_map_norm)
diamond_volume_fraction = x_z_vol_frac(x_z_vol_frac>0.20&x_z_vol_frac<0.35);
percentage_diamond = nnz(diamond_volume_fraction)/nnz(x_z_vol_frac);

disordered_volume_fraction = x_z_vol_frac(x_z_vol_frac>0.44&x_z_vol_frac<0.55);
percentage_disordered = nnz(disordered_volume_fraction)/nnz(x_z_vol_frac);

primitive_volume_fraction = x_z_vol_frac(x_z_vol_frac>0.34&x_z_vol_frac<0.44);
percentage_primitive = nnz(primitive_volume_fraction)/nnz(x_z_vol_frac);

nnz(diamond_volume_fraction)
nnz(x_z_vol_frac)
numel(x_z_vol_frac)


%% percentage of volume with specific CRT values
image01 = image0.*3.72;

%diamond
diamond = image01;
diamond(diamond<4.36|diamond>6.48) = 0; %need to do this first!
diamond(diamond>4.36&diamond<6.48) = 1;
imagesc(squeeze(diamond(500,:,:))); axis equal off
colormap("gray")
hold on

%disordered
disordered = image01;
disordered(disordered<10.41|disordered>11.47) = 0;
disordered(disordered>10.41&disordered<11.47) = 2;
dis = imagesc(squeeze(disordered(500,:,:))); axis equal off
alpha(dis,'color')

%primitive
primitive = image01;
primitive(primitive<6.31|primitive>8.44) = 0;
primitive(primitive>6.31&primitive<8.44) = 3;
prim = imagesc(squeeze(primitive(500,:,:))); axis equal off
alpha(prim,'color')

diamond_trabecluae_width = image01(image01>4.36&image01<6.48);
percentage_diamond = nnz(diamond_trabecluae_width)/numel(image01(image01~=-3.72));

disordered_trabecluae_width = image01(image01>10.41&image01<11.47);
percentage_disordered = nnz(disordered_trabecluae_width)/numel(image01(image01~=-3.72));

primitive_trabecluae_width = image01(image01>6.314&image01<8.44);
percentage_primitive = nnz(primitive_trabecluae_width)/numel(image01(image01~=-3.72));

figure(1);
ax1 = axes;
imagesc(squeeze(primitive(500,:,:)),'AlphaData',0.5); axis equal off %displays volume fraction for that slice of the volume
set(ax1,'YDir','normal')
colormap(ax1,"jet")
caxis([0 0.6])
box off

ax2 = axes;
imagesc(squeeze((:,550,:)),'AlphaData',0.3); axis equal tight; %displays slice through y
% imagesc(squeeze(binary(:,:,1)),'AlphaData',0.5) %displays slice through z
set(ax2,'YDir','normal')
box off

linkaxes([ax1,ax2])

ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

colormap(ax2,"gray")

colorbar(ax1,'Position',[0.85 0.15 0.025 0.7])


