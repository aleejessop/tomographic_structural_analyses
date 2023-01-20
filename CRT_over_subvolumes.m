clear
filename = ('C:\Users\20220428\OneDrive - Murdoch University\Documents\MATLAB\Sea_Urchin\data\greyscale\1_69_inner_interambulacral_plate_transformed.tif');
greyscale = tiffreadVolume(filename);
binary = imbinarize(greyscale,0.51);
binary = binary(201:1200,181:1180,301:1300);
binary = 1-binary;
imagesc(squeeze(binary(800,:,:))); axis equal
size_subvolume = 50;
[subvolume] = subvolume_extraction(binary, size_subvolume);
volshow(squeeze(subvolume(1,1,2,:,:,:)))

for a = 1:size(subvolume,1)
    for b = 1:size(subvolume,2)
        for c = 1:size(subvolume,3)
            if mean(subvolume(a,b,c,:,:,:)) == 0
                subvolume(a,b,c,:,:,:) = 1;
            else
                subvolume(a,b,c,:,:,:) = subvolume(a,b,c,:,:,:);
            end
        end
    end
end

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

volshow(squeeze(image0(1,1,3,:,:,:)))

maxCoveringDistanceTransform(squeeze(subvolume(5,5,5,:,:,:)));