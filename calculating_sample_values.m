%values for Table 1

clear
load("CRT_solid_phase_169.mat")
padded_CRT = padarray(image0,[105 279 238],0,'post');
padded_CRT = padarray(padded_CRT,[249 299 249],0,'pre');
padded_D = padarray(D,[105 279 238],0,'post');
padded_D = padarray(padded_D,[249 299 249],0,'pre');
imagesc(squeeze(padded_D(:,:,724)),'AlphaData',0.3); axis equal tight;
set(gca,'YDir','normal')

load('solid_volume_fraction.mat')

size_subvolume = 50;
resized_volume_fraction = imresize3(volume_fraction,size_subvolume,"nearest");
padded_volume_fraction = padarray(resized_volume_fraction,[5 29 38],0,'post');
imagesc(squeeze(padded_volume_fraction(:,:,760)),'AlphaData',1); axis equal tight;

data = image0; sliceview

ss =  padded_CRT(730:1250,796:1332,535:1015);
max(ss(:))
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

data = diamsample1_D_pore_phase; sliceview

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

imagesc(squeeze(diamsample1_D_pore_phase(10,:,:)))
maxEDM = max(diamsample1_D_solid_phase(:))*pixel_res; % in microns
diamsample1_CRT_solid_phase(diamsample1_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample1_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample1_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% diamond 1 high resolution
sample_size = 218;
pixel_res = 0.732;
diam_sample1_pos = [196 1142-339 247];
imagesc(squeeze(binary(:,669,:)))
padded_CRT = image0;
padded_D = D;

size_subvolume = 50;
resized_volume_fraction = imresize3(volume_fraction,size_subvolume,"nearest");
padded_volume_fraction = padarray(resized_volume_fraction,[46 42 3],0,'post');
imagesc(squeeze(padded_volume_fraction(:,669,:)),'AlphaData',0.3); axis equal tight;

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

imagesc(squeeze(diamsample1_D_solid_phase(100,:,:)))
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

data = diamsample2_D_pore_phase; sliceview

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

imagesc(squeeze(diamsample2_D_solid_phase(5,:,:)))
maxEDM = max(diamsample2_D_solid_phase(:))*pixel_res; % in microns
diamsample2_CRT_solid_phase(diamsample2_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample2_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample2_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% diamond 3
diam_sample3_pos = [733 1779-566 483];

%pore phase
diamsample3_CRT_pore_phase = padded_CRT(diam_sample3_pos(1):diam_sample3_pos(1)+sample_size,...
    diam_sample3_pos(2):diam_sample3_pos(2)+sample_size,diam_sample3_pos(3):diam_sample3_pos(3)+sample_size);
diamsample3_D_pore_phase = padded_D(diam_sample3_pos(1):diam_sample3_pos(1)+sample_size,...
    diam_sample3_pos(2):diam_sample3_pos(2)+sample_size,diam_sample3_pos(3):diam_sample3_pos(3)+sample_size);
diamsample3_solid_volume = padded_volume_fraction(diam_sample3_pos(1):diam_sample3_pos(1)+sample_size,...
    diam_sample3_pos(2):diam_sample3_pos(2)+sample_size,diam_sample3_pos(3):diam_sample3_pos(3)+sample_size);

data = diamsample2_D_pore_phase; sliceview

solid_volume_fraction = mean(diamsample3_solid_volume,'all');
maxEDM = max(diamsample3_D_pore_phase(:))*pixel_res; % in microns
diamsample3_CRT_pore_phase(diamsample3_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample3_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample3_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%solid phase
diamsample3_CRT_solid_phase = padded_CRT(diam_sample3_pos(1):diam_sample3_pos(1)+sample_size,...
    diam_sample3_pos(2):diam_sample3_pos(2)+sample_size,diam_sample3_pos(3):diam_sample3_pos(3)+sample_size);
diamsample3_D_solid_phase = padded_D(diam_sample3_pos(1):diam_sample3_pos(1)+sample_size,...
    diam_sample3_pos(2):diam_sample3_pos(2)+sample_size,diam_sample3_pos(3):diam_sample3_pos(3)+sample_size);
data = diamsample3_D_solid_phase; sliceview

imagesc(squeeze(diamsample3_D_solid_phase(5,:,:)))
maxEDM = max(diamsample3_D_solid_phase(:))*pixel_res; % in microns
diamsample3_CRT_solid_phase(diamsample3_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample3_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample3_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% diamond 4
diam_sample4_pos = [834 1779-746 1144];

%pore phase
diamsample4_CRT_pore_phase = padded_CRT(diam_sample4_pos(1):diam_sample4_pos(1)+sample_size,...
    diam_sample4_pos(2):diam_sample4_pos(2)+sample_size,diam_sample4_pos(3):diam_sample4_pos(3)+sample_size);
diamsample4_D_pore_phase = padded_D(diam_sample4_pos(1):diam_sample4_pos(1)+sample_size,...
    diam_sample4_pos(2):diam_sample4_pos(2)+sample_size,diam_sample4_pos(3):diam_sample4_pos(3)+sample_size);
diamsample4_solid_volume = padded_volume_fraction(diam_sample4_pos(1):diam_sample4_pos(1)+sample_size,...
    diam_sample4_pos(2):diam_sample4_pos(2)+sample_size,diam_sample4_pos(3):diam_sample4_pos(3)+sample_size);

data = diamsample4_D_pore_phase; sliceview

solid_volume_fraction = mean(diamsample4_solid_volume,'all');
maxEDM = max(diamsample4_D_pore_phase(:))*pixel_res; % in microns
diamsample4_CRT_pore_phase(diamsample4_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample4_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample4_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%solid phase
diamsample4_CRT_solid_phase = padded_CRT(diam_sample4_pos(1):diam_sample4_pos(1)+sample_size,...
    diam_sample4_pos(2):diam_sample4_pos(2)+sample_size,diam_sample4_pos(3):diam_sample4_pos(3)+sample_size);
diamsample4_D_solid_phase = padded_D(diam_sample4_pos(1):diam_sample4_pos(1)+sample_size,...
    diam_sample4_pos(2):diam_sample4_pos(2)+sample_size,diam_sample4_pos(3):diam_sample4_pos(3)+sample_size);
data = diamsample4_D_solid_phase; sliceview

imagesc(squeeze(diamsample4_CRT_solid_phase(5,:,:)))
maxEDM = max(diamsample4_D_solid_phase(:))*pixel_res; % in microns
diamsample4_CRT_solid_phase(diamsample4_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample4_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample4_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% diamond 5
diam_sample5_pos = [834 1779-943 463];

%pore phase
diamsample5_CRT_pore_phase = padded_CRT(diam_sample5_pos(1):diam_sample5_pos(1)+sample_size,...
    diam_sample5_pos(2):diam_sample5_pos(2)+sample_size,diam_sample5_pos(3):diam_sample5_pos(3)+sample_size);
diamsample5_D_pore_phase = padded_D(diam_sample5_pos(1):diam_sample5_pos(1)+sample_size,...
    diam_sample5_pos(2):diam_sample5_pos(2)+sample_size,diam_sample5_pos(3):diam_sample5_pos(3)+sample_size);
diamsample5_solid_volume = padded_volume_fraction(diam_sample5_pos(1):diam_sample5_pos(1)+sample_size,...
    diam_sample5_pos(2):diam_sample5_pos(2)+sample_size,diam_sample5_pos(3):diam_sample5_pos(3)+sample_size);

data = diamsample5_D_pore_phase; sliceview

solid_volume_fraction = mean(diamsample5_solid_volume,'all');
maxEDM = max(diamsample5_D_pore_phase(:))*pixel_res; % in microns
diamsample5_CRT_pore_phase(diamsample5_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample5_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample5_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%solid phase
diamsample5_CRT_solid_phase = padded_CRT(diam_sample5_pos(1):diam_sample5_pos(1)+sample_size,...
    diam_sample5_pos(2):diam_sample5_pos(2)+sample_size,diam_sample5_pos(3):diam_sample5_pos(3)+sample_size);
diamsample5_D_solid_phase = padded_D(diam_sample5_pos(1):diam_sample5_pos(1)+sample_size,...
    diam_sample5_pos(2):diam_sample5_pos(2)+sample_size,diam_sample5_pos(3):diam_sample5_pos(3)+sample_size);
data = diamsample5_D_solid_phase; sliceview

imagesc(squeeze(diamsample5_CRT_solid_phase(5,:,:)))
maxEDM = max(diamsample5_D_solid_phase(:))*pixel_res; % in microns
diamsample5_CRT_solid_phase(diamsample5_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample5_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample5_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% diamond 6*
diam_sample6_pos = [747 1779-462 574];

%pore phase
diamsample6_CRT_pore_phase = padded_CRT(diam_sample6_pos(1):diam_sample6_pos(1)+sample_size,...
    diam_sample6_pos(2):diam_sample6_pos(2)+sample_size,diam_sample6_pos(3):diam_sample6_pos(3)+sample_size);
diamsample6_D_pore_phase = padded_D(diam_sample6_pos(1):diam_sample6_pos(1)+sample_size,...
    diam_sample6_pos(2):diam_sample6_pos(2)+sample_size,diam_sample6_pos(3):diam_sample6_pos(3)+sample_size);
diamsample6_solid_volume = padded_volume_fraction(diam_sample6_pos(1):diam_sample6_pos(1)+sample_size,...
    diam_sample6_pos(2):diam_sample6_pos(2)+sample_size,diam_sample6_pos(3):diam_sample6_pos(3)+sample_size);

data = diamsample6_D_pore_phase; sliceview

solid_volume_fraction = mean(diamsample6_solid_volume,'all');
maxEDM = max(diamsample6_D_pore_phase(:))*pixel_res; % in microns
diamsample6_CRT_pore_phase(diamsample6_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample6_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample6_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%solid phase
diamsample6_CRT_solid_phase = padded_CRT(diam_sample6_pos(1):diam_sample6_pos(1)+sample_size,...
    diam_sample6_pos(2):diam_sample6_pos(2)+sample_size,diam_sample6_pos(3):diam_sample6_pos(3)+sample_size);
diamsample6_D_solid_phase = padded_D(diam_sample6_pos(1):diam_sample6_pos(1)+sample_size,...
    diam_sample6_pos(2):diam_sample6_pos(2)+sample_size,diam_sample6_pos(3):diam_sample6_pos(3)+sample_size);
data = diamsample6_D_solid_phase; sliceview

imagesc(squeeze(diamsample6_CRT_solid_phase(5,:,:)))
maxEDM = max(diamsample6_D_solid_phase(:))*pixel_res; % in microns
diamsample6_CRT_solid_phase(diamsample6_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample6_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample6_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% diamond 7
diam_sample7_pos = [805 1779-515 724];

%pore phase
diamsample7_CRT_pore_phase = padded_CRT(diam_sample7_pos(1):diam_sample7_pos(1)+sample_size,...
    diam_sample7_pos(2):diam_sample7_pos(2)+sample_size,diam_sample7_pos(3):diam_sample7_pos(3)+sample_size);
diamsample7_D_pore_phase = padded_D(diam_sample7_pos(1):diam_sample7_pos(1)+sample_size,...
    diam_sample7_pos(2):diam_sample7_pos(2)+sample_size,diam_sample7_pos(3):diam_sample7_pos(3)+sample_size);
diamsample7_solid_volume = padded_volume_fraction(diam_sample7_pos(1):diam_sample7_pos(1)+sample_size,...
    diam_sample7_pos(2):diam_sample7_pos(2)+sample_size,diam_sample7_pos(3):diam_sample7_pos(3)+sample_size);

data = diamsample7_D_pore_phase; sliceview

solid_volume_fraction = mean(diamsample7_solid_volume,'all');
maxEDM = max(diamsample7_D_pore_phase(:))*pixel_res; % in microns
diamsample7_CRT_pore_phase(diamsample7_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample7_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample7_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%solid phase
diamsample7_CRT_solid_phase = padded_CRT(diam_sample7_pos(1):diam_sample7_pos(1)+sample_size,...
    diam_sample7_pos(2):diam_sample7_pos(2)+sample_size,diam_sample7_pos(3):diam_sample7_pos(3)+sample_size);
diamsample7_D_solid_phase = padded_D(diam_sample7_pos(1):diam_sample7_pos(1)+sample_size,...
    diam_sample7_pos(2):diam_sample7_pos(2)+sample_size,diam_sample7_pos(3):diam_sample7_pos(3)+sample_size);
data = diamsample7_D_pore_phase; sliceview

imagesc(squeeze(diamsample7_solid_volume(5,:,:)))
maxEDM = max(diamsample7_D_solid_phase(:))*pixel_res; % in microns
diamsample7_CRT_solid_phase(diamsample7_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample7_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample7_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% diamond 8*
diam_sample8_pos = [819 1779-604 1203];

%pore phase
diamsample8_CRT_pore_phase = padded_CRT(diam_sample8_pos(1):diam_sample8_pos(1)+sample_size,...
    diam_sample8_pos(2):diam_sample8_pos(2)+sample_size,diam_sample8_pos(3):diam_sample8_pos(3)+sample_size);
diamsample8_D_pore_phase = padded_D(diam_sample8_pos(1):diam_sample8_pos(1)+sample_size,...
    diam_sample8_pos(2):diam_sample8_pos(2)+sample_size,diam_sample8_pos(3):diam_sample8_pos(3)+sample_size);
diamsample8_solid_volume = padded_volume_fraction(diam_sample8_pos(1):diam_sample8_pos(1)+sample_size,...
    diam_sample8_pos(2):diam_sample8_pos(2)+sample_size,diam_sample8_pos(3):diam_sample8_pos(3)+sample_size);

data = diamsample8_D_pore_phase; sliceview

solid_volume_fraction = mean(diamsample8_solid_volume,'all');
maxEDM = max(diamsample8_D_pore_phase(:))*pixel_res; % in microns
diamsample8_CRT_pore_phase(diamsample8_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample8_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample8_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%solid phase
diamsample8_CRT_solid_phase = padded_CRT(diam_sample8_pos(1):diam_sample8_pos(1)+sample_size,...
    diam_sample8_pos(2):diam_sample8_pos(2)+sample_size,diam_sample8_pos(3):diam_sample8_pos(3)+sample_size);
diamsample8_D_solid_phase = padded_D(diam_sample8_pos(1):diam_sample8_pos(1)+sample_size,...
    diam_sample8_pos(2):diam_sample8_pos(2)+sample_size,diam_sample8_pos(3):diam_sample8_pos(3)+sample_size);
data = diamsample8_D_pore_phase; sliceview

imagesc(squeeze(diamsample8_solid_volume(5,:,:)))
maxEDM = max(diamsample8_D_solid_phase(:))*pixel_res; % in microns
diamsample8_CRT_solid_phase(diamsample8_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(diamsample8_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(diamsample8_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% disordered 1

disordered_sample_pos = [575 1779-690 666];

% pore phase
dissample1_CRT_pore_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample1_D_pore_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample1_solid_volume_pore_phase = padded_volume_fraction(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

data = dissample1_CRT_pore_phase; sliceview

solid_volume_fraction = mean(dissample1_solid_volume_pore_phase,'all');
maxEDM = max(dissample1_D_pore_phase(:))*pixel_res; % in microns
dissample1_CRT_pore_phase(dissample1_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample1_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample1_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
dissample1_CRT_solid_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample1_D_solid_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample1_solid_volume_solid_phase = padded_volume_fraction(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

data = dissample1_CRT_solid_phase; sliceview

solid_volume_fraction = mean(dissample1_solid_volume_solid_phase,'all');
maxEDM = max(dissample1_D_solid_phase(:))*pixel_res; % in microns
dissample1_CRT_solid_phase(dissample1_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample1_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample1_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% disordered 2

sample_size = 94;
pixel_res = 1.69;
disordered_sample_pos = [470 1779-690 1036];

% pore phase
dissample2_CRT_pore_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample2_D_pore_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample2_solid_volume_pore_phase = padded_volume_fraction(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

data = dissample2_CRT_solid_phase; sliceview

solid_volume_fraction = mean(dissample2_solid_volume_pore_phase,'all');
maxEDM = max(dissample2_D_pore_phase(:))*pixel_res; % in microns
dissample2_CRT_pore_phase(dissample2_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample2_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample2_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
dissample2_CRT_solid_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample2_D_solid_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

maxEDM = max(dissample2_D_solid_phase(:))*pixel_res; % in microns
dissample2_CRT_solid_phase(dissample2_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample2_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample2_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% disordered 3

sample_size = 94;
pixel_res = 1.69;
disordered_sample_pos = [473 1779-1176 1036];

% pore phase
dissample3_CRT_pore_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample3_D_pore_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample3_solid_volume_pore_phase = padded_volume_fraction(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

data = dissample3_CRT_pore_phase; sliceview

solid_volume_fraction = mean(dissample3_solid_volume_pore_phase,'all');
maxEDM = max(dissample3_D_pore_phase(:))*pixel_res; % in microns
dissample3_CRT_pore_phase(dissample3_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample3_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample3_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
dissample3_CRT_solid_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample3_D_solid_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

maxEDM = max(dissample3_D_solid_phase(:))*pixel_res; % in microns
dissample3_CRT_solid_phase(dissample3_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample3_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample3_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% disordered 4

sample_size = 94;
pixel_res = 1.69;
disordered_sample_pos = [490 1779-1112 547];

% pore phase
dissample4_CRT_pore_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample4_D_pore_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample4_solid_volume_pore_phase = padded_volume_fraction(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

data = dissample4_CRT_solid_phase; sliceview

solid_volume_fraction = mean(dissample4_solid_volume_pore_phase,'all');
maxEDM = max(dissample4_D_pore_phase(:))*pixel_res; % in microns
dissample4_CRT_pore_phase(dissample4_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample4_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample4_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
dissample4_CRT_solid_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample4_D_solid_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

maxEDM = max(dissample4_D_solid_phase(:))*pixel_res; % in microns
dissample4_CRT_solid_phase(dissample4_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample4_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample4_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% disordered 5

sample_size = 94;
pixel_res = 1.69;
disordered_sample_pos = [509 1779-900 587];

% pore phase
dissample5_CRT_pore_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample5_D_pore_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample5_solid_volume_pore_phase = padded_volume_fraction(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

data = dissample5_CRT_solid_phase; sliceview

solid_volume_fraction = mean(dissample5_solid_volume_pore_phase,'all');
maxEDM = max(dissample5_D_pore_phase(:))*pixel_res; % in microns
dissample5_CRT_pore_phase(dissample5_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample5_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample5_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
dissample5_CRT_solid_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample5_D_solid_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

maxEDM = max(dissample5_D_solid_phase(:))*pixel_res; % in microns
dissample5_CRT_solid_phase(dissample5_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample5_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample5_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% disordered 6

sample_size = 94;
pixel_res = 1.69;
disordered_sample_pos = [490 1779-663 483];

% pore phase
dissample5_CRT_pore_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample5_D_pore_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample5_solid_volume_pore_phase = padded_volume_fraction(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

data = dissample5_CRT_solid_phase; sliceview

solid_volume_fraction = mean(dissample5_solid_volume_pore_phase,'all');
maxEDM = max(dissample5_D_pore_phase(:))*pixel_res; % in microns
dissample5_CRT_pore_phase(dissample5_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample5_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample5_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
dissample5_CRT_solid_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample5_D_solid_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

maxEDM = max(dissample5_D_solid_phase(:))*pixel_res; % in microns
dissample5_CRT_solid_phase(dissample5_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample5_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample5_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% disordered 7

sample_size = 94;
pixel_res = 1.69;
disordered_sample_pos = [481 1779-1142 803];

% pore phase
dissample5_CRT_pore_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample5_D_pore_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample5_solid_volume_pore_phase = padded_volume_fraction(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

data = dissample5_CRT_solid_phase; sliceview

solid_volume_fraction = mean(dissample5_solid_volume_pore_phase,'all');
maxEDM = max(dissample5_D_pore_phase(:))*pixel_res; % in microns
dissample5_CRT_pore_phase(dissample5_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample5_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample5_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
dissample5_CRT_solid_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample5_D_solid_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

maxEDM = max(dissample5_D_solid_phase(:))*pixel_res; % in microns
dissample5_CRT_solid_phase(dissample5_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample5_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample5_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% disordered 8

sample_size = 94;
pixel_res = 1.69;
disordered_sample_pos = [462 1779-1192 368];

% pore phase
dissample5_CRT_pore_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample5_D_pore_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample5_solid_volume_pore_phase = padded_volume_fraction(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

data = dissample5_CRT_solid_phase; sliceview

solid_volume_fraction = mean(dissample5_solid_volume_pore_phase,'all');
maxEDM = max(dissample5_D_pore_phase(:))*pixel_res; % in microns
dissample5_CRT_pore_phase(dissample5_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample5_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample5_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
dissample5_CRT_solid_phase = padded_CRT(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);
dissample5_D_solid_phase = padded_D(disordered_sample_pos(1):disordered_sample_pos(1)+sample_size,...
    disordered_sample_pos(2):disordered_sample_pos(2)+sample_size,disordered_sample_pos(3):disordered_sample_pos(3)+sample_size);

maxEDM = max(dissample5_D_solid_phase(:))*pixel_res; % in microns
dissample5_CRT_solid_phase(dissample5_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(dissample5_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(dissample5_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% primitive 1

sample_size = 94;
pixel_res = 1.69;
primitive_sample_pos = [1022 1779-906 494];

% pore phase
primitive1_CRT_pore_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive1_D_pore_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive1_solid_volume_pore_phase = padded_volume_fraction(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

data = primitive1_CRT_solid_phase; sliceview

solid_volume_fraction = mean(primitive1_solid_volume_pore_phase,'all');
maxEDM = max(primitive1_D_pore_phase(:))*pixel_res; % in microns
primitive1_CRT_pore_phase(primitive1_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive1_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive1_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
primitive1_CRT_solid_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive1_D_solid_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

maxEDM = max(primitive1_D_solid_phase(:))*pixel_res; % in microns
primitive1_CRT_solid_phase(primitive1_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive1_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive1_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% primitive 2

sample_size = 94;
pixel_res = 1.69;
primitive_sample_pos = [1031 1779-1052 643];

% pore phase
primitive2_CRT_pore_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_pore_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_solid_volume_pore_phase = padded_volume_fraction(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

data = primitive2_CRT_pore_phase; sliceview

solid_volume_fraction = mean(primitive2_solid_volume_pore_phase,'all');
maxEDM = max(primitive2_D_pore_phase(:))*pixel_res; % in microns
primitive2_CRT_pore_phase(primitive2_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
primitive2_CRT_solid_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_solid_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

maxEDM = max(primitive2_D_solid_phase(:))*pixel_res; % in microns
primitive2_CRT_solid_phase(primitive2_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% primitive 3

sample_size = 94;
pixel_res = 1.69;
primitive_sample_pos = [1044 1779-961 813];

% pore phase
primitive2_CRT_pore_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_pore_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_solid_volume_pore_phase = padded_volume_fraction(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

data = primitive2_solid_volume_pore_phase; sliceview

solid_volume_fraction = mean(primitive2_solid_volume_pore_phase,'all');
maxEDM = max(primitive2_D_pore_phase(:))*pixel_res; % in microns
primitive2_CRT_pore_phase(primitive2_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
primitive2_CRT_solid_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_solid_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

maxEDM = max(primitive2_D_solid_phase(:))*pixel_res; % in microns
primitive2_CRT_solid_phase(primitive2_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% primitive 4

sample_size = 94;
pixel_res = 1.69;
primitive_sample_pos = [952 1779-1014 978];

% pore phase
primitive2_CRT_pore_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_pore_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_solid_volume_pore_phase = padded_volume_fraction(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

data = primitive2_CRT_solid_phase; sliceview

solid_volume_fraction = mean(primitive2_solid_volume_pore_phase,'all');
maxEDM = max(primitive2_D_pore_phase(:))*pixel_res; % in microns
primitive2_CRT_pore_phase(primitive2_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
primitive2_CRT_solid_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_solid_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

maxEDM = max(primitive2_D_solid_phase(:))*pixel_res; % in microns
primitive2_CRT_solid_phase(primitive2_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% primitive 5

sample_size = 94;
pixel_res = 1.69;
primitive_sample_pos = [903 1779-766 378];

% pore phase
primitive2_CRT_pore_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_pore_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_solid_volume_pore_phase = padded_volume_fraction(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

data = primitive2_CRT_solid_phase; sliceview

solid_volume_fraction = mean(primitive2_solid_volume_pore_phase,'all');
maxEDM = max(primitive2_D_pore_phase(:))*pixel_res; % in microns
primitive2_CRT_pore_phase(primitive2_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
primitive2_CRT_solid_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_solid_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

maxEDM = max(primitive2_D_solid_phase(:))*pixel_res; % in microns
primitive2_CRT_solid_phase(primitive2_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% primitive 6

sample_size = 94;
pixel_res = 1.69;
primitive_sample_pos = [859 1779-1051 374];

% pore phase
primitive2_CRT_pore_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_pore_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_solid_volume_pore_phase = padded_volume_fraction(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

data = primitive2_CRT_solid_phase; sliceview

solid_volume_fraction = mean(primitive2_solid_volume_pore_phase,'all');
maxEDM = max(primitive2_D_pore_phase(:))*pixel_res; % in microns
primitive2_CRT_pore_phase(primitive2_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
primitive2_CRT_solid_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_solid_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

maxEDM = max(primitive2_D_solid_phase(:))*pixel_res; % in microns
primitive2_CRT_solid_phase(primitive2_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% primitive 7

sample_size = 94;
pixel_res = 1.69;
primitive_sample_pos = [840 1779-1358 985];

% pore phase
primitive2_CRT_pore_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_pore_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_solid_volume_pore_phase = padded_volume_fraction(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

data = primitive2_CRT_solid_phase; sliceview

solid_volume_fraction = mean(primitive2_solid_volume_pore_phase,'all');
maxEDM = max(primitive2_D_pore_phase(:))*pixel_res; % in microns
primitive2_CRT_pore_phase(primitive2_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
primitive2_CRT_solid_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_solid_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

maxEDM = max(primitive2_D_solid_phase(:))*pixel_res; % in microns
primitive2_CRT_solid_phase(primitive2_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% primitive 8

sample_size = 94;
pixel_res = 1.69;
primitive_sample_pos = [906 1779-1251 985];

% pore phase
primitive2_CRT_pore_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_pore_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_solid_volume_pore_phase = padded_volume_fraction(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

data = primitive2_CRT_solid_phase; sliceview

solid_volume_fraction = mean(primitive2_solid_volume_pore_phase,'all');
maxEDM = max(primitive2_D_pore_phase(:))*pixel_res; % in microns
primitive2_CRT_pore_phase(primitive2_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
primitive2_CRT_solid_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_solid_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

maxEDM = max(primitive2_D_solid_phase(:))*pixel_res; % in microns
primitive2_CRT_solid_phase(primitive2_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

%% DiamondtoPrimitive1

sample_size = 94;
pixel_res = 1.69;
primitive_sample_pos = [781 1779-999 628];

% pore phase
primitive2_CRT_pore_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_pore_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_solid_volume_pore_phase = padded_volume_fraction(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

data = primitive2_CRT_solid_phase; sliceview

solid_volume_fraction = mean(primitive2_solid_volume_pore_phase,'all');
maxEDM = max(primitive2_D_pore_phase(:))*pixel_res; % in microns
primitive2_CRT_pore_phase(primitive2_CRT_pore_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_pore_phase(:),'omitnan')*pixel_res; %standard deviation in microns


% solid phase
primitive2_CRT_solid_phase = padded_CRT(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);
primitive2_D_solid_phase = padded_D(primitive_sample_pos(1):primitive_sample_pos(1)+sample_size,...
    primitive_sample_pos(2):primitive_sample_pos(2)+sample_size,primitive_sample_pos(3):primitive_sample_pos(3)+sample_size);

maxEDM = max(primitive2_D_solid_phase(:))*pixel_res; % in microns
primitive2_CRT_solid_phase(primitive2_CRT_solid_phase==-1)=NaN; %change all negative 1 values to nans
meanCRT = mean(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %mean in microns
sdCRT = std(primitive2_CRT_solid_phase(:),'omitnan')*pixel_res; %standard deviation in microns

