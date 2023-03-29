%plotting CRT maps and histograms
clear
load("CRT_pore_phase.mat")
cropped_image0 = image0(1:1125,1:1050,1:1125);
figure;
imagesc(squeeze(image0(:,497,:))); axis equal
box off
colormap('hot')
colorbar
clim([0 5])
hold on
pixel_res = 3.72;
diam_sample1_pos = round([1737.41 2067.07 2202.68]/pixel_res);
rectangle('Position',[475 510 43 43])



