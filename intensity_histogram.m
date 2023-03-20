clear
%read in tif volume and show
greyscale = tiffreadVolume('C:\Users\20220428\OneDrive - Murdoch University\Documents\MATLAB\Sea_Urchin\data\greyscale\primitive.tiff'); 
volshow(greyscale);
%convert uint16 to double
intensity_values = double(greyscale)./(2^16);
%plot intensity values as histogram
figure, histogram(greyscale)
% xlim([0 1]);
% ylim([0 6*10^6])

%remove zeros from the matrix but turns into single vector!!
intensity_values(intensity_values==0)=[];
figure, histogram(intensity_values)
box off
xlim([0 1])

%computes global theshold T from grayscale image using Otsu's method:
%Otsu's method chooses a threshold that minimizes the intraclass variance of the thresholded black and white pixels.
%In other words it searches for the threshold intensity (level) which
%maximizes the between class variance
[level effectiveness] = graythresh(intensity_values)
%for the sea urchin data level = 0.5059
%plot histogram with all pixels greater than level as white and less than
%level black
white_pixels = intensity_values(intensity_values>=level);
black_pixels = intensity_values(intensity_values<=level);

figure, histogram(white_pixels,'FaceColor',[1 1 1])
hold on
histogram(black_pixels,'FaceColor',[0 0 0])
%now plot a line showing the threshold level
xline(0.46,'LineWidth',1)
box off

% binarize image using this level and show
binary = imbinarize(greyscale,0.51);
figure, volshow(binary)
%binarize image using different levels and show
binary_data2 = imbinarize(intensity_values,0.51);
figure, volshow(binary_data2)

intensityvaluesofthevalley = intensity_values>0.45 & intensity_values<0.5;
intensityvaluesofthevalley = intensity_values(intensityvaluesofthevalley);
figure,  h = histogram(intensityvaluesofthevalley)