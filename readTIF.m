
% filename = 'C:\Users\20220428\OneDrive - Murdoch University\Documents\MATLAB\Sea_Urchin\data\greyscale\Sea Urchin Diamond_111.tiff';
% pixelregionx = [1 100];
% pixelregiony = [1 100];
% pixelregionz = [1 100];
% binarize = 1;

function [greyscale, binary] = readTIF(filename,pixelregionx,pixelregiony,pixelregionz,binarize)
% addpath('C:\Users\20220428\OneDrive - Murdoch University\Documents\MATLAB\Sea_Urchin\data')
greyscale = tiffreadVolume(filename, 'PixelRegion',{pixelregionx, pixelregiony, pixelregionz});
if binarize == 1
    binary = imbinarize(greyscale,0.51); %0.5059 %0.48
else
    binary = [];
end
end
    
   





