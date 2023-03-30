
%% simulated data

N=15;
PixelSize=1/N;
Threshold=0.;
NVoxels=4*round([sqrt(3)*N,sqrt(2)*N,sqrt(6)*N]);
binary=createNodalSurface("p",NVoxels,PixelSize,1,[1,1,1],[1 -1 0],[0,0,0],Threshold);

% muck up data set to cut off corner
for i=1:size(binary,1)
    for j=1:size(binary,2)
        for k=1:size(binary,3)
            if i<j
                binary(i,j,k)=0;
            end
            if (j-NVoxels(2)/2)^2+(i-NVoxels(3)/4)^2 > (0.5*mean(NVoxels))^2
                binary(i,j,k)=0;
            end
        end
    end
end 

disp("Data set created.");

%% sea urchin data
clear
filename = 'C:\Users\20220428\OneDrive - Murdoch University\Documents\MATLAB\Sea_Urchin\data\greyscale\1_69_inner_interambulacral_plate.tif';
greyscale = tiffreadVolume(filename);
binary = imbinarize(greyscale,0.51);
binary = 1-binary;


%% run algorithm

paramsCRT.maximalCRTValueToCompute = 11;
paramsCRT.calculateCRTOnlyForRectangularSubset=1;
paramsCRT.rectangularSubsetForCRTCalculation=[[3,3,3]',(size(binary)-7)'];
[image0,image1,D,statusCRTComputation]=maxCoveringDistanceTransform(binary,paramsCRT);

%data=statusCRTComputation; sliceview
%data=image0; sliceview

