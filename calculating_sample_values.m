%values for table of samples

clear
load("CRT_pore_phase.mat")
load('solid_volume_fraction.mat')
pixel_res = 3.72;
diam_sample1_pos = round([1737.41 2067.07 2202.68]/pixel_res); %sample position in image coordinates
disordered_sample_pos = round([1896.14 2100.36 1765.71]/pixel_res);
size_subvolume = 25;
resized_volume_fraction = imresize3(volume_fraction,size_subvolume,"linear");
padded_volume_fraction = padarray(resized_volume_fraction,[23 7 17],0,'post');
imagesc(squeeze(padded_volume_fraction(:,497,:)),'AlphaData',0.3); axis equal tight;

sample_size = 42;
%diamond 1
diamsample1_CRT_pore_phase = image0(diam_sample1_pos(1):diam_sample1_pos(1)+sample_size,...
    diam_sample1_pos(2):diam_sample1_pos(2)+sample_size,diam_sample1_pos(3):diam_sample1_pos(3)+sample_size);
diamsample1_D_pore_phase = D(diam_sample1_pos(1):diam_sample1_pos(1)+sample_size,...
    diam_sample1_pos(2):diam_sample1_pos(2)+sample_size,diam_sample1_pos(3):diam_sample1_pos(3)+sample_size);
diamsample1_solid_volume_pore_phase = padded_volume_fraction(diam_sample1_pos(1):diam_sample1_pos(1)+sample_size,...
    diam_sample1_pos(2):diam_sample1_pos(2)+sample_size,diam_sample1_pos(3):diam_sample1_pos(3)+sample_size);

solid_volume_fraction = mean(sample1_solid_volume_pore_phase,'all');
maxEDM = max(sample1_D_pore_phase(:))*3.72; % in microns
sample1_CRT_pore_phase(sample1_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(sample1_CRT_pore_phase(:),'omitnan')*3.72; %mean in microns
sdCRT = std(sample1_CRT_pore_phase(:),'omitnan')*3.72; %standard deviation in microns

%disordered 1
dissample1_CRT_pore_phase = image0(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample1_D_pore_phase = D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample1_solid_volume_pore_phase = padded_volume_fraction(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

solid_volume_fraction = mean(dissample1_solid_volume_pore_phase,'all');
maxEDM = max(dissample1_D_pore_phase(:))*3.72; % in microns
sample1_CRT_pore_phase(dissample1_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample1_CRT_pore_phase(:),'omitnan')*3.72; %mean in microns
sdCRT = std(dissample1_CRT_pore_phase(:),'omitnan')*3.72; %standard deviation in microns