%filename = 'xxxx.tif';
%pixelregion = [26 75];
%binarize = 1

function [greyscale, binary] = readTIF(filename,pixelregionx,pixelregiony,pixelregionz,binarize)
% addpath('C:\Users\20220428\OneDrive - Murdoch University\Documents\MATLAB\Sea_Urchin\data')
greyscale = tiffreadVolume(filename, 'PixelRegion',{pixelregionx, pixelregiony, pixelregionz});
if binarize == 1
    binary = imbinarize(greyscale,0.51); %0.5059 %0.48
else
    binary = [];
end
end
    
   





