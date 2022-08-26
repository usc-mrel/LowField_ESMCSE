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

%% This scripts estimate T2* from MEGRE (DICOM)
DICOM_MEGRE_path = '/Users/bli/MREL/LungT2/Data/vol416/DICOM/MEGRE_EE';
DICOM_files = dir(fullfile(DICOM_MEGRE_path,'*.IMA'));

for idx = 1:length(DICOM_files)
    filepath = fullfile(DICOM_MEGRE_path,DICOM_files(idx).name);
    info{idx,:} = dicominfo(filepath);
    if idx <= length(DICOM_files)/2
        img_mag(:,:,idx) = single(dicomread(filepath));
    else
        y = info{idx,:}.RescaleIntercept;
        x = info{idx,:}.RescaleSlope;
        img_ph(:,:,idx-length(DICOM_files)/2) = x * single(dicomread(filepath)) + y;
        % Convert to radians
        img_rad = ((single(img_ph) - y) / 8191 * 360 -180) /180 *pi ;
    end
    TE_MEGRE(idx) =  info{idx}.EchoTime;
end
img_MEGRE_EE = img_mag .* exp(1i*img_rad);

%% Smooth DICOM
for k = 1: size(img_MEGRE_EE,3)
    img_MEGRE_EE_real_sm(:,:,k) = imgaussfilt(real(img_MEGRE_EE(:,:,k)),1);
    img_MEGRE_EE_imag_sm(:,:,k) = imgaussfilt(imag(img_MEGRE_EE(:,:,k)),1);
    img_MEGRE_EE_abs_sm(:,:,k)  = imgaussfilt(abs(img_MEGRE_EE(:,:,k)),1);
end
img_MEGRE_EE_sm = img_MEGRE_EE_real_sm + 1i*img_MEGRE_EE_imag_sm;



%% Fitting
clear signal_TE
megre_idx = 1:4;
TE_MEGRE_used = TE_MEGRE(megre_idx).';

idx = 0;
for k = megre_idx
    idx = idx + 1;
    temp = (img_MEGRE_EE_abs_sm(:,:,k));
    signal_TE(:,idx) = temp(bw_roi); % each columne for one TE; each row for one voxel in ROI
end
signal_TE = double(signal_TE);

tic
parfor p = 1:size(signal_TE,1)
    p
    y = abs(signal_TE(p,:)).';
    %y = [real(signal_TE(p,:)).'; imag(signal_TE(p,:)).'];
    x0 = [abs(signal_TE(p,1)), 1/150];
    %x0 = [real(signal_TE(p,1)),imag(signal_TE(p,1)), 1/150,0];
    E = @(x) rho_R2StarMEGRE_fitting(x,y,TE_MEGRE_used);
    %E = @(x) rho_T2StarMEGRE_complexfitting(x,y,TE_MEGRE_used);
    options = optimoptions('lsqnonlin');
    [x_megre_abs_sm_4(p,:),resnorm,residual,exitflag,output]= lsqnonlin(E,x0,[-Inf,0],[Inf,1],options);
    %[x_megre_complex_4te(p,:),resnorm,residual,exitflag,output]= lsqnonlin(E,x0,[-Inf,-Inf,0,-Inf],[Inf,Inf,1,Inf],options);
end
toc

%%
R2Starmap_megre_abs_sm_4te = zeros(size(img_MEGRE_EE_abs_sm,[1 2]));
R2Starmap_megre_abs_sm_4te(bw_roi)     = x_megre_abs_sm_4(:,2)*1000;    % [Hz]

im_file = sprintf('/R2Star_MEGRE_EE');
fprintf('Saving R2* images in %s...\n',[DICOM_MEGRE_path,im_file])
save([DICOM_MEGRE_path,im_file], 'R2Starmap_megre_abs_sm_4te','bw_roi','-v7.3');

