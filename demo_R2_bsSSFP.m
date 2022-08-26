%% Set source directories
computer_type = computer;
if strcmp(computer_type, 'MACI64') % MAC
    src_directory = '/Users/bli/MREL/LungT2/Code/lowfield_lung/mcse/src';
    thirdparty_directory = '/Users/bli/MREL/LungT2/Code/lowfield_lung/mcse/thirdparty';
    ismrmrd_directory = '/Users/bli/ismrmrd';
elseif strcmp(computer_type, 'PCWIN64') % Windows
    src_drectory = 'E:\lowfield_lung\src';
    thirdparty_directory = 'E:\gradient_nonlinearity\thirdparty';
    ismrmrd_directory = 'D:\ismrmrd\ismrmrd';
elseif strcmp(computer_type, 'GLNXA64') % Server
    %src_directory = '/server/home/nlee/GNL_lowrank_recon_3d_gpu';
    %thirdparty_directory = 'E:\gradient_nonlinearity\thirdparty';
    ismrmrd_directory = '/server/home/nlee/ismrmrd';
end

%% Add source directories to search path
addpath(genpath(src_directory));
addpath(genpath(thirdparty_directory));
addpath(genpath(ismrmrd_directory));


% This scripts estimate T2 from T2-prep bSSFP (DICOM)
DICOM_bssfp_path = '/Users/bli/MREL/LungT2/Data/vol416/DICOM/bssfp_EE';
DICOM_files = dir(fullfile(DICOM_bssfp_path,'*.IMA'));

for idx = 1:length(DICOM_files)
    filepath = fullfile(DICOM_bssfp_path,DICOM_files(idx).name);
    info_bssfp{idx,:} = dicominfo(filepath);
    img_bssfp(:,:,idx) = single(dicomread(filepath));
    TE_bssfp(idx) = info_bssfp{idx}.EchoTime;
end

%% Smooth DICOM
for k = 1: size(img_bssfp,3)
    img_bssfp_sm(:,:,k) = imgaussfilt(img_bssfp(:,:,k),1);   
end

%% Fitting
clear signal_TE
T2prep = [0 25 55].';

idx = 0;
for k = 1: length(T2prep)
    idx = idx + 1;
    temp = abs(img_bssfp_sm(:,:,k));
    signal_TE(:,idx) = temp(bw_roi); % each columne for one TE; each row for one voxel in ROI
end
signal_TE = double(signal_TE);

tic
parfor p = 1:size(signal_TE,1)
    p
    y = abs(signal_TE(p,:)).';
    x0 = [abs(signal_TE(p,1)), 1/60];
    E = @(x) rho_R2bSSFP_fitting(x,T2prep,y);
    options = optimoptions('lsqnonlin');
    [x_bSSFP(p,:),resnorm,residual,exitflag,output]= lsqnonlin(E,x0,[0,0],[],options);
end
toc

%%
R2map_bssfp = zeros(size(img_bssfp_sm,[1 2]));
R2map_bssfp(bw_roi)       = x_bSSFP(:,2)*1000;    % [Hz]
S0_bssfp = zeros(size(img_bssfp_sm,[1 2]));
S0_bssfp(bw_roi)    = x_bSSFP(:,1);

im_file = sprintf('/R2_bssfp_EE');
fprintf('Saving R2 images in %s...\n',[DICOM_bssfp_path,im_file])
save([DICOM_bssfp_path,im_file], 'R2map_bssfp','bw_roi','-v7.3');
