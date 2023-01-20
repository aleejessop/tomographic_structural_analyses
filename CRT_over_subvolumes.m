clear
filename = ('C:\Users\20220428\OneDrive - Murdoch University\Documents\MATLAB\Sea_Urchin\data\greyscale\1_69_inner_interambulacral_plate_transformed_cropped_for_CRT.tiff');
greyscale = tiffreadVolume(filename);
binary = imbinarize(greyscale,0.51);
% binary = binary(201:1200,181:1180,301:1300);
binary = binary(1:400,1:1200,1:1000);
binary = 1-binary;
imagesc(squeeze(binary(200,:,:))); axis equal
size_subvolume = 50;
[subvolume] = subvolume_extraction(binary, size_subvolume);
volshow(squeeze(subvolume(1,1,2,:,:,:)))

image0 = zeros(20,20,20,50,50,50);
image1 = zeros(20,20,20,50,50,50);

for i = 1:size(subvolume,1)
    i
    for j = 1:size(subvolume,2)
        j
        for k = 1:size(subvolume,3)
            k
            [image0(i,j,k,:,:,:),image1(i,j,k,:,:,:)] = maxCoveringDistanceTransform(squeeze(subvolume(i,j,k,:,:,:)));
        end
    end
end

CRT_data.pore_sizes_solid = image0;
CRT_data.filename = filename;
CRT_data.subvolume_size = 50;
CRT_data_pore_sizes_solid_max_centre = image1;
CRT_data.voxels = [1:400,1:1200,1:1000];
CRT_data.pore_sizes_solid_max = image1;

save('CRT_cropped_data.mat','CRT_data','-v7.3')
