
% clear
% filename = ('C:\Users\20220428\OneDrive - Murdoch University\Documents\MATLAB\Sea_Urchin\data\interambulacral_plate.tif');
% greyscale = tiffreadVolume(filename);
% pixelregionx = [1 1024];
% pixelregiony = [1 1003];
% pixelregionz = [1 1016];
% 
% [greyscale, binary] = readTIF(filename,pixelregionx,pixelregiony,pixelregionz,1);
s=size(data);

figure; %2D slices of binary_data
for i = 1:s(3)
    clf;
%     imagesc(squeeze(skeleton(:,:,i))); axis equal tight
%     hold on
    imagesc(squeeze(data(:,i,:)),'AlphaData',0.5); axis equal tight
    colormap hot
    colorbar
    drawnow;
    hold on
    pause
end


