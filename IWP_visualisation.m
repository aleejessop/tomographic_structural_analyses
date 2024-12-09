
filename = 'C:\Users\20220428\OneDrive - Murdoch University\IWP\Sas_male_greenblue_isotropic15_ordered.tif';
ordered = tiffreadVolume(filename);
ordered = double(ordered);
volshow(ordered)
imshow(squeeze(ordered(:,:,10)))

