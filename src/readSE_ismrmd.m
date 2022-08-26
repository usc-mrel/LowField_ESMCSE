function [c_img, kdata, noise, info] = readSE_ismrmd(flags)
%   Read multiple SE dataset from ISMRMD .h5 and noise.h5
%   Input:
%       - flags: all info for reading

%   Output:
%   Saved as cell format, each cell contains one TE SE k-space or
%   images
%       - c_img: cell of coil images (Nx x Ny x Ncoil x Nslice)
%       - kdata: cell of coil k-space (Nx x Ny x Nslice x Ncoil)
%       - noise covariance (Ncoil x Ncoil)
%       - info: basic informatino of raw data

%   Author: Bochao Li
%   Email: bochaoli@usc.edu

%% Define constants
gamma = 2.67522212e+8; % gyromagnetic ratio for 1H [rad/sec/T]

%% Read an ismrmrd file
h5_fullpath = dir(flags.h5_fileList);
noise_fullpath = dir(flags.noise_fileList);
nF = length(h5_fullpath);
start_time = tic;

for nF = 1:length(h5_fullpath)
    h5_FileName = h5_fullpath(nF).name;
    h5_fullFileName = fullfile(h5_fullpath(nF).folder, h5_FileName);
    noise_FileNama = noise_fullpath(nF).name;
    noise_fullFileName = fullfile(noise_fullpath(nF).folder, noise_FileNama);
    
    tic; fprintf('Reading an ismrmrd file: %s... ', h5_fullFileName);
    if exist(h5_fullFileName, 'file')
        dset = ismrmrd.Dataset(h5_fullFileName, 'dataset');
        fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    else
        error('File %s does not exist.  Please generate it.' , cartesian_fullpath);
    end
    
    %% ---- Get imaging parameters from the XML header ----
    header = ismrmrd.xml.deserialize(dset.readxml);
    
    % ---- Sequence parameters ----
    info{nF}.TR         = header.sequenceParameters.TR;     % [msec]
    info{nF}.TE         = header.sequenceParameters.TE;     % [msec]
%     info{nF}.ESP        = header.sequenceParameters.echo_spacing; % [msec]
    info{nF}.flip_angle = header.sequenceParameters.flipAngle_deg; % [degrees]
    
    % ---- Experimental conditions ----
    info{nF}.B0 = header.experimentalConditions.H1resonanceFrequency_Hz * (2 * pi / gamma); % [Hz] * [2pi rad/cycle] / [rad/sec/T] => [T]
    
    % ---- Ecoding info ----
    info{nF}.field_of_view_mm(1) = header.encoding.reconSpace.fieldOfView_mm.x; % RO
    info{nF}.field_of_view_mm(2) = header.encoding.reconSpace.fieldOfView_mm.y; % PE
    info{nF}.field_of_view_mm(3) = header.encoding.reconSpace.fieldOfView_mm.z; % SS
    
    info{nF}.Nkx = header.encoding.encodedSpace.matrixSize.x; % number of samples in k-space (RO,row)
    info{nF}.Nky = header.encoding.encodedSpace.matrixSize.y; % number of samples in k-space (PE,col)
    info{nF}.Nkz = header.encoding.encodedSpace.matrixSize.z; % number of samples in k-space (SS,slice)
    info{nF}.Nx  = header.encoding.reconSpace.matrixSize.x;   % number of samples in image-space (RO,row)
    info{nF}.Ny  = header.encoding.reconSpace.matrixSize.y;   % number of samples in image-space (PE,col)
    info{nF}.Nz  = header.encoding.reconSpace.matrixSize.z;   % number of samples in image-space (SS,slice)
    
    %% ---- Parse the ISMRMRD header ----
    tic; fprintf('Parsing the ISMRMRD header... ');
    temp_raw_data = dset.readAcquisition(); % read all the acquisitions
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    
    info{nF}.rawData.number_of_samples  = single(max(temp_raw_data.head.number_of_samples));
    info{nF}.rawData.discard_pre        = single(max(temp_raw_data.head.discard_pre));
    info{nF}.rawData.discard_post       = single(max(temp_raw_data.head.discard_post));
    info{nF}.rawData.nr_samples         = info{nF}.rawData.number_of_samples - info{nF}.rawData.discard_pre  - info{nF}.rawData.discard_post;
    info{nF}.rawData.nr_channels        = single(max(temp_raw_data.head.active_channels));
    info{nF}.rawData.nr_phase_encodings = single(max(temp_raw_data.head.idx.kspace_encode_step_1)) + 1; % nr_interleaves for spiral imaging
    info{nF}.rawData.nr_slice_encodings = single(max(temp_raw_data.head.idx.kspace_encode_step_2)) + 1;
    info{nF}.rawData.nr_averages        = single(max(temp_raw_data.head.idx.average)) + 1;
    info{nF}.rawData.nr_slices          = single(max(temp_raw_data.head.idx.slice)) + 1; % imaging slice number
    info{nF}.rawData.nr_contrasts       = single(max(temp_raw_data.head.idx.contrast)) + 1; % echo number in multi-echo
    info{nF}.rawData.nr_phases          = single(max(temp_raw_data.head.idx.phase)) + 1; % cardiac phase
    info{nF}.rawData. nr_repetitions    = single(max(temp_raw_data.head.idx.repetition)) + 1; % dynamic number for dynamic scanning
    info{nF}.rawData.nr_sets            = single(max(temp_raw_data.head.idx.set)) + 1; % flow encoding set
    info{nF}.rawData.nr_segments        = single(max(temp_raw_data.head.idx.segment)) + 1;
    info{nF}.rawData.trajectory_dimensions = single(max(temp_raw_data.head.trajectory_dimensions)); % Get the dimensionality of the trajectory vector (0 means no trajectory)
    info{nF}.rawData.dt = single(max(temp_raw_data.head.sample_time_us)) * 1e-3; %  dwell time in [sec][usec] * [msec/1e-3 usec] => [msec]
    
    %% ---- Get all nominal k-space trajectories and k-space data ----
    tic; fprintf('Sorting out all k-space data... ');
    
    % ---- Initialize k-space: [samples, phase, slice_encodings, averages, slices, contrasts, cardiac_phases, repetitions, sets, segments, coils ]
    data = complex(zeros(info{nF}.rawData.nr_samples, info{nF}.rawData.nr_phase_encodings, info{nF}.rawData.nr_slice_encodings,...
        info{nF}.rawData.nr_averages, info{nF}.rawData.nr_slices, info{nF}.rawData.nr_contrasts, info{nF}.rawData.nr_phases,...
        info{nF}.rawData.nr_repetitions, info{nF}.rawData.nr_sets, info{nF}.rawData.nr_segments, info{nF}.rawData.nr_channels, 'double'));
    
    
    info{nF}.rawData.nr_profiles = length(temp_raw_data.data);
    for idx = 1:info{nF}.rawData.nr_profiles
        % ---- Get information about the current profile ----
        phase_encoding_nr = single(temp_raw_data.head.idx.kspace_encode_step_1(idx)) + 1;
        slice_encoding_nr = single(temp_raw_data.head.idx.kspace_encode_step_2(idx)) + 1;
        average_nr        = single(temp_raw_data.head.idx.average(idx)) + 1;
        slice_nr          = single(temp_raw_data.head.idx.slice(idx)) + 1;
        contrast_nr       = single(temp_raw_data.head.idx.contrast(idx)) + 1;
        phase_nr          = single(temp_raw_data.head.idx.phase(idx)) + 1;
        repetition_nr     = single(temp_raw_data.head.idx.repetition(idx)) + 1;
        set_nr            = single(temp_raw_data.head.idx.set(idx)) + 1;
        segment_nr        = single(temp_raw_data.head.idx.segment(idx)) + 1;
        
        % ---- Set the index range of a profile ----
        number_of_samples = single(temp_raw_data.head.number_of_samples(idx));
        discard_pre       = single(temp_raw_data.head.discard_pre(idx));
        discard_post      = single(temp_raw_data.head.discard_post(idx));
        start_index       = discard_pre + 1;
        end_index         = number_of_samples - discard_post;
        index_range       = (start_index:end_index).';
        
        % ---- Get k-space data ----
        full_profile = temp_raw_data.data{idx}; % number_of_samples x nr_channels
        data(index_range, phase_encoding_nr, slice_encoding_nr, average_nr, slice_nr,...
            contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr, :) = full_profile(index_range,:); % nr_samples x nr_channels
    end
    
    %% ---- Concatnate segmenets of Phase Encoding ----
%     if flags.IgnoreSeg
%         effectTE_idx = round(info{nF}.TE / info{nF}.ESP);
%         seg_idx = 1:info{nF}.rawData.nr_segments;
%         seg_idx = seg_idx(seg_idx~=effectTE_idx);
%         temp_kdata = sum(data(:, :, :, :, :, :, :, :, :,seg_idx, :), 10);
%         temp_kdata(:, temp_raw_data.head.idx.kspace_encode_step_1(1) + 1,:,:,:,:,:,:,:,:,:) = 0;
%         temp_kdata = temp_kdata + data(:, :, :, :, :, :, :, :, :,effectTE_idx, :);
%         kdata{nF} = temp_kdata;
%     else
%         kdata{nF} = data;
%     end
    
    %% ---- Crop k-space (Phse Encoding direction) same as recon space
    kdata{nF}  = data;
    if flags.CropPhaseEncoding
        remain_start = header.encoding.encodingLimits.kspace_encoding_step_1.center - header.encoding.reconSpace.matrixSize.y/2 + 1;
        remain_end   = header.encoding.encodingLimits.kspace_encoding_step_1.center + (header.encoding.reconSpace.matrixSize.y - 1 - header.encoding.reconSpace.matrixSize.y/2) + 1;
        kdata{nF} = kdata{nF}(:, [remain_start : remain_end], :, :, :, :, :, :, :, :, :);
    end
    
    %% ---- Squeeze kspace ----
    if flags.Squeeze
        kdata{nF} = squeeze(kdata{nF});
        kdata{nF} = reshape(kdata{nF},[info{nF}.rawData.nr_samples, info{nF}.rawData.nr_phase_encodings, info{nF}.rawData.nr_slice_encodings, info{nF}.rawData.nr_channels ]);
        % ---- Remove oversampling ----
        if flags.RemoveOS
            Nkx = size(kdata{nF},1);
            Nx = Nkx/2;
            imc_ = fftshift(fft(ifftshift(kdata{nF}, 1), [], 1), 1); % Nkx x Nky x Nslice x Ncoil
            idx1_range = (-floor(Nx/2):ceil(Nx/2)-1).' + floor(Nkx/2) + 1;
            imc_ = imc_(idx1_range,:,:,:);  % Nx x Nky x Nslice x Nc
            kdata{nF} = fftshift(ifft(ifftshift(imc_, 1), [], 1), 1); % Nkx x Nky x Nslice x Nc
        end
    end
    
    %% Process noise only ismrmrd data
    [noise{nF}] = process_noise_ismrmrd_data(noise_fullFileName);
    
    %% FFT to Image (k-space <=> image-space)
    c_img{nF} = fft2c(kdata{nF});
end

end

