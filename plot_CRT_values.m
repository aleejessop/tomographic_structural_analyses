%plotting CRT maps and histograms
clear
load("CRT_solid_phase.mat")
cropped_image0 = image0(1:1125,1:1050,1:1125);
figure;
imagesc(squeeze(cropped_image0(:,525,:))); axis equal
box off
colormap('hot')
colorbar
clim([0 5])

bins = 0:6;

sample_through_z1 = discretize(squeeze(image0(355,525,1:1125)),bins);
sample_through_z2 = discretize(squeeze(image0(555,525,1:1125)),bins);
sample_through_z3 = discretize(squeeze(image0(755,525,1:1125)),bins);

figure;
z_axis = 1:1125;
bar(z_axis,sample_through_z1)
box off
hold on
