%values for table of samples

clear
load("CRT_void_space_169.mat")
load('solid_volume_fraction.mat')

size_subvolume = 25;
resized_volume_fraction = imresize3(volume_fraction,size_subvolume,"nearest");
padded_volume_fraction = padarray(resized_volume_fraction,[13 15 19],0,'post');
imagesc(squeeze(padded_volume_fraction(:,669,:)),'AlphaData',0.3); axis equal tight;

diam_sample1_pos = [574 669 232]; %check you are using the correct values here sometimes x and y are swapped
hold on
rectangle('Position',[232 574 118 118])

sample_size = 117;
pixel_res = 1.69;

%diamond 1
diamsample1_CRT_pore_phase = image0(diam_sample1_pos(1):diam_sample1_pos(1)+sample_size,...
    diam_sample1_pos(2):diam_sample1_pos(2)+sample_size,diam_sample1_pos(3):diam_sample1_pos(3)+sample_size);
diamsample1_D_pore_phase = D(diam_sample1_pos(1):diam_sample1_pos(1)+sample_size,...
    diam_sample1_pos(2):diam_sample1_pos(2)+sample_size,diam_sample1_pos(3):diam_sample1_pos(3)+sample_size);
diamsample1_solid_volume_pore_phase = padded_volume_fraction(diam_sample1_pos(1):diam_sample1_pos(1)+sample_size,...
    diam_sample1_pos(2):diam_sample1_pos(2)+sample_size,diam_sample1_pos(3):diam_sample1_pos(3)+sample_size);

solid_volume_fraction = mean(diamsample1_solid_volume_pore_phase,'all');
maxEDM = max(diamsample1_D_pore_phase(:))*pixel_res; % in microns
diamsample1_CRT_pore_phase(diamsample1_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample1_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample1_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%disordered 1
dissample1_CRT_pore_phase = image0(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample1_D_pore_phase = D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample1_solid_volume_pore_phase = padded_volume_fraction(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

solid_volume_fraction = mean(dissample1_solid_volume_pore_phase,'all');
maxEDM = max(dissample1_D_pore_phase(:))*pixel_res; % in microns
sample1_CRT_pore_phase(dissample1_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample1_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample1_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns