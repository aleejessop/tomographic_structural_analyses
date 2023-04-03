%values for table of samples

clear
load("CRT_solid_phase_169.mat")
padded_CRT = padarray(image0,[105 279 238],0,'post');
padded_CRT = padarray(padded_CRT,[249 299 249],0,'pre');
padded_D = padarray(D,[105 279 238],0,'post');
padded_D = padarray(padded_D,[249 299 249],0,'pre');
imagesc(squeeze(padded_CRT(:,669,:)),'AlphaData',0.3); axis equal tight;
set(gca,'YDir','normal')

load('solid_volume_fraction.mat')

size_subvolume = 50;
resized_volume_fraction = imresize3(volume_fraction,size_subvolume,"nearest");
padded_volume_fraction = padarray(resized_volume_fraction,[5 29 38],0,'post');
imagesc(squeeze(padded_volume_fraction(:,669,:)),'AlphaData',0.3); axis equal tight;

 %check you are using the correct values here sometimes x and y are swapped

hold on
rectangle('Position',[667 1166 95 95])

sample_size = 94;
pixel_res = 1.69;

%% diamond 1
diam_sample1_pos = [834 1779-614 667];
%pore phase
diamsample1_CRT_pore_phase = padded_CRT(diam_sample1_pos(1):diam_sample1_pos(1)+sample_size,...
    diam_sample1_pos(2):diam_sample1_pos(2)+sample_size,diam_sample1_pos(3):diam_sample1_pos(3)+sample_size);
diamsample1_D_pore_phase = padded_D(diam_sample1_pos(1):diam_sample1_pos(1)+sample_size,...
    diam_sample1_pos(2):diam_sample1_pos(2)+sample_size,diam_sample1_pos(3):diam_sample1_pos(3)+sample_size);
diamsample1_solid_volume = padded_volume_fraction(diam_sample1_pos(1):diam_sample1_pos(1)+sample_size,...
    diam_sample1_pos(2):diam_sample1_pos(2)+sample_size,diam_sample1_pos(3):diam_sample1_pos(3)+sample_size);

solid_volume_fraction = mean(diamsample1_solid_volume,'all');
maxEDM = max(diamsample1_D_pore_phase(:))*pixel_res; % in microns
diamsample1_CRT_pore_phase(diamsample1_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample1_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample1_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%solid phase
diamsample1_CRT_solid_phase = padded_CRT(diam_sample1_pos(1):diam_sample1_pos(1)+sample_size,...
    diam_sample1_pos(2):diam_sample1_pos(2)+sample_size,diam_sample1_pos(3):diam_sample1_pos(3)+sample_size);
diamsample1_D_solid_phase = padded_D(diam_sample1_pos(1):diam_sample1_pos(1)+sample_size,...
    diam_sample1_pos(2):diam_sample1_pos(2)+sample_size,diam_sample1_pos(3):diam_sample1_pos(3)+sample_size);
data = diamsample1_D_solid_phase; sliceview

imagesc(squeeze(diamsample1_D_solid_phase(10,:,:)))
maxEDM = max(diamsample1_D_solid_phase(:))*pixel_res; % in microns
diamsample1_CRT_solid_phase(diamsample1_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample1_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample1_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% diamond 2
diam_sample2_pos = [760 1779-832 628];

%pore phase
diamsample2_CRT_pore_phase = padded_CRT(diam_sample2_pos(1):diam_sample2_pos(1)+sample_size,...
    diam_sample2_pos(2):diam_sample2_pos(2)+sample_size,diam_sample2_pos(3):diam_sample2_pos(3)+sample_size);
diamsample2_D_pore_phase = padded_D(diam_sample2_pos(1):diam_sample2_pos(1)+sample_size,...
    diam_sample2_pos(2):diam_sample2_pos(2)+sample_size,diam_sample2_pos(3):diam_sample2_pos(3)+sample_size);
diamsample2_solid_volume = padded_volume_fraction(diam_sample2_pos(1):diam_sample2_pos(1)+sample_size,...
    diam_sample2_pos(2):diam_sample2_pos(2)+sample_size,diam_sample2_pos(3):diam_sample2_pos(3)+sample_size);

solid_volume_fraction = mean(diamsample2_solid_volume,'all');
maxEDM = max(diamsample2_D_pore_phase(:))*pixel_res; % in microns
diamsample2_CRT_pore_phase(diamsample2_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample2_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample2_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%solid phase
diamsample2_CRT_solid_phase = padded_CRT(diam_sample2_pos(1):diam_sample2_pos(1)+sample_size,...
    diam_sample2_pos(2):diam_sample2_pos(2)+sample_size,diam_sample2_pos(3):diam_sample2_pos(3)+sample_size);
diamsample2_D_solid_phase = padded_D(diam_sample2_pos(1):diam_sample2_pos(1)+sample_size,...
    diam_sample2_pos(2):diam_sample2_pos(2)+sample_size,diam_sample2_pos(3):diam_sample2_pos(3)+sample_size);
data = diamsample2_D_solid_phase; sliceview

imagesc(squeeze(diamsample2_D_solid_phase(10,:,:)))
maxEDM = max(diamsample2_D_solid_phase(:))*pixel_res; % in microns
diamsample2_CRT_solid_phase(diamsample2_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample2_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample2_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns


%% disordered 1
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