% demo_siemens_mcse_lowrank_subspace_recon.m
% Written by Bochao Li
% Email: bochaoli@usc.edu
% Started: 07/15/2021, Last modified: 08/26/2022

%% Clean slate
close all; clear; clc;
start_time = tic;

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

%% Define data directory
computer_type = computer;
if strcmp(computer_type, 'MACI64')
    data_directory = '/Users/bli/MREL/LungT2/Data/vol375/Raw';
elseif strcmp(computer_type, 'PCWIN64')
    data_directory = 'D:\mri_data\disc\lung\S1';
elseif strcmp(computer_type, 'GLNXA64')
    data_directory = '/mnt/sdata_new/nlee/mri_data/disc/brain';
end

%% Define user options
user_opts.remove_ROoversampling = 1; % remove readout oversampling: 1=yes, 0=no
user_opts.remove_PEoversampling = 1; % remove readout oversampling: 1=yes, 0=no
R = 1;

%% Read a multi-contrast SE dataset (5 in 1)
data_filename          = 'meas_MID00061_FID38893_estse5in1_cor_shift024_esp20_mcse8_R1_exhale';
ismrmrd_noise_fullpath = fullfile(data_directory, sprintf('/noise/noise_%s.h5', data_filename));
ismrmrd_data_fullpath  = fullfile(data_directory, sprintf('h5/%s.h5', data_filename));
siemens_dat_fullpath   = fullfile(data_directory, sprintf('%s.dat', data_filename));
[kspace_all, header, sampling_pattern] = siemens_read_navigated_mcse5in1_kspace(siemens_dat_fullpath, ismrmrd_data_fullpath, ismrmrd_noise_fullpath, user_opts);

%% Get TE
TE = header.sequenceParameters.TE.' * 1e-3; % [msec] * [sec/1e3msec] => [sec]

%% Calculate the dimensions of k-space data
% Nk x Nky x Nkz x Nc x M => Nk x Nky x Nkz x M x Nc
kspace_all  = permute(kspace_all,  [1 2 3 5 6 4]);
[Nk,Nky,Nkz,Necho,Nshift,Nc] = size(kspace_all);
N1 = Nk;
N2 = Nky;
N3 = Nkz;
N = N1 * N2 * N3;

%% Calculate coil sensitivity maps using Walsh's method
%--------------------------------------------------------------------------
% Calculate the size of an autocalibration region
%--------------------------------------------------------------------------
calib_size = 24;
cal_shape = ones(1, 3, 'double');
if size(kspace_all,3) == 1
    nr_dims = 2;
else
    nr_dims = 3;
end
cal_shape(1:nr_dims) = calib_size;

%--------------------------------------------------------------------------
% Calculate the calibration region of k-space
%--------------------------------------------------------------------------
tstart = tic; fprintf('Calculating the calibration region of k-space (size = %d)... ', calib_size);
cal_data = crop(reshape(kspace_all(:,:,:,1,3,:), [Nk Nky Nkz Nc]), [cal_shape Nc]);
hanning_window = bsxfun(@times, hanning(cal_shape(1)) * hanning(cal_shape(2)).', reshape(hanning(cal_shape(3)), [1 1 cal_shape(3)]));
cal_data = bsxfun(@times, cal_data, hanning_window);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Apply FFT along the slice direction
%--------------------------------------------------------------------------
tstart = tic; fprintf('Applying FFT along the slice direction... ');
cal_data = zpad(cal_data, [cal_shape Nc]);
cal_data = 1 / sqrt(N3) * fftshift(fft(ifftshift(cal_data, 3), [], 3), 3);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Calculate coil sensitivity maps (k-space <=> image-space)
%--------------------------------------------------------------------------
csm = complex(zeros(Nk, Nky, Nkz, Nc, 'double'));
for idx3 = 1:Nkz
    tstart = tic; fprintf('Calculating csm using Walsh method (%d/%d)... ', idx3, Nkz);
    cal_im = zpad(reshape(cal_data(:,:,idx3,:), [cal_shape(1:2) Nc]), [Nk Nky Nc]);
    for dim = 1:2
        cal_im = 1 / sqrt(size(cal_im,dim)) * fftshift(fft(ifftshift(cal_im, dim), [], dim), dim);
    end
    csm(:,:,idx3,:) = reshape(ismrm_estimate_csm_walsh(cal_im), [Nk Nky 1 Nc]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Perform FFT reconstruction
tstart = tic; fprintf('Performing FFT reconstruction... ');
%--------------------------------------------------------------------------
% FFT to image-space (k-space <=> image-space)
%--------------------------------------------------------------------------
imc_fft = kspace_all;
for dim = 1:nr_dims
    imc_fft = 1 / sqrt(size(imc_fft,dim)) * fftshift(fft(ifftshift(imc_fft, dim), [], dim), dim); % N1 x N2 x N3
end

%--------------------------------------------------------------------------
% Perform adaptive coil combination
%--------------------------------------------------------------------------
im_fft = complex(zeros(N1, N2, N3, Necho,Nshift, 'double'));
for n_shift = 1:Nshift
    for n_se = 1:Necho
        im = complex(zeros(N1, N2, N3, 'double'));
        for c = 1:Nc
            im = im + bsxfun(@times, conj(csm(:,:,:,c)), imc_fft(:,:,:,n_se,n_shift,c)); % N1 x N2 x N3
        end
        im_fft(:,:,:,n_se,n_shift) = im ./ sqrt(sum(abs(csm).^2,4));
    end
end
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Re-order temporal dimension (8x5 --> 40)
nt = 0;
for n_se = 1: Necho
    for n_shift = 1:Nshift
        nt = nt+1;
        im_fft_signal(:,:,:,nt)  = im_fft(:,:,:,n_se,n_shift);
    end
end

%% Smoothing images (gauss)
for k = 1: size(im_fft_signal,4)
    img_real_sm(:,:,:,k) = imgaussfilt(real(im_fft_signal(:,:,:,k)),1);
    img_imag_sm(:,:,:,k) = imgaussfilt(imag(im_fft_signal(:,:,:,k)),1);
end
img_sm = img_real_sm + 1i*img_imag_sm;

%% Draw mask
% figure,
% imagesc(abs(img_sm(:,:,1,8)),[0 25]);
% axis off image
% colormap(gray)
% 
% for n_lung = 1:2
%     hf = drawfreehand;
%     bw_lung(:,:,n_lung) = hf.createMask();
% end
% bw_draw= logical(sum(bw_lung,3));
% threshold_bw = (abs(img_sm(:,:,1,8)) < 10);
% bw_roi = bw_draw & threshold_bw ;

%% Peform fitting
tstart_fitting = tic; fprintf('Start fitting... ');
data_idx = [6:40];

im_fitting  =   img_sm(:,:,:,data_idx);
t_shift     =   repmat([-4,-2,0,2,4].',[Necho,1]) ;
t_shift     =   t_shift(data_idx);
t_SE        =   reshape(repmat(TE*1e3,[1,Nshift]).',[],1);
t_SE        =   t_SE(data_idx);

clear x_out_sm signal_TE  resnorm residual exitflag output
for k = 1: size(im_fitting,4)
    temp = im_fitting(:,:,:,k);
    signal_TE(:,k) = temp(:); % each columne for one TE; each row for one voxel in ROI
end
signal_TE = double(signal_TE);

tic
parfor p = 1:size(signal_TE,1)
    p
    y = [real(signal_TE(p,:)).'; imag(signal_TE(p,:)).'];
    x0 = [real(signal_TE(p,3)), imag(signal_TE(p,3)), 1/60, 1/10, 0].'; % inital values
    E = @(x) mixed_rho_r2_r2prime_B0_fitting_esmcse(x, y, t_shift, t_SE);
    options = optimoptions('lsqnonlin');
    [x_out_sm(p,:),resnorm(p),residual(:,p),exitflag,output] = lsqnonlin(E,x0,[-Inf,-Inf,0,1/100,-Inf],[Inf,Inf,1,1/4,Inf],options);
end
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart_fitting), toc(start_time));

%% Get the results
R2map_sm = zeros(size(im_fft_signal,[1 2]));
R2primemap_sm = zeros(size(im_fft_signal,[1 2]));
B0map_sm = zeros(size(im_fft_signal,[1 2]));
S0map_sm = zeros(size(im_fft_signal,[1 2]));

R2map_sm  = reshape(x_out_sm(:,3),[128,128])*1000;    % [Hz]
R2primemap_sm  = reshape(x_out_sm(:,4),[128,128])*1000;    % [Hz]
B0map_sm = reshape(x_out_sm(:,5),[128,128]);

mean(R2map_sm(bw_roi))
mean(R2primemap_sm(bw_roi))
mean(R2primemap_sm(bw_roi) + R2map_sm(bw_roi))

std(R2map_sm(bw_roi))
std(R2primemap_sm(bw_roi))
std(R2primemap_sm(bw_roi) + R2map_sm(bw_roi))