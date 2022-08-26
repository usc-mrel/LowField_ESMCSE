function [kspace, sampling_pattern] = siemens_read_tse_shuffled_kspace(siemens_dat_fullpath, ismrmrd_data_fullpath, ismrmrd_noise_fullpath, user_opts)
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 06/01/2021, Last modified: 06/01/2022

%% Read k-space data (ISMRMRD format)
start_time = tic;
tstart = tic; fprintf('Reading an ISMRMRD file: %s... ', ismrmrd_data_fullpath);
if exist(ismrmrd_data_fullpath, 'file')
    dset = ismrmrd.Dataset(ismrmrd_data_fullpath, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_data_fullpath);
end

%% Get imaging parameters from the XML header
header = ismrmrd.xml.deserialize(dset.readxml);

%--------------------------------------------------------------------------
% measurement information
%--------------------------------------------------------------------------
patient_position = header.measurementInformation.patientPosition;

%--------------------------------------------------------------------------
% Sequence parameters
%--------------------------------------------------------------------------
TR         = header.sequenceParameters.TR * 1e-3;     % [msec] * [sec/1e3msec] => [sec]
TE         = header.sequenceParameters.TE * 1e-3;     % [msec] * [sec/1e3msec] => [sec]
flip_angle = header.sequenceParameters.flipAngle_deg; % [degrees]

if strcmp(header.sequenceParameters.sequence_type, 'TurboSpinEcho')
    echo_spacing = header.sequenceParameters.echo_spacing * 1e-3; % [msec] * [sec/1e3msec] => [sec]
end

%--------------------------------------------------------------------------
% Encoding
%--------------------------------------------------------------------------
encoded_fov(1) = header.encoding.encodedSpace.fieldOfView_mm.x; % RO
encoded_fov(2) = header.encoding.encodedSpace.fieldOfView_mm.y; % PE
encoded_fov(3) = header.encoding.encodedSpace.fieldOfView_mm.z; % SL

recon_fov(1) = header.encoding.reconSpace.fieldOfView_mm.x; % RO
recon_fov(2) = header.encoding.reconSpace.fieldOfView_mm.y; % PE
recon_fov(3) = header.encoding.reconSpace.fieldOfView_mm.z; % SL

Nkx = header.encoding.encodedSpace.matrixSize.x; % number of readout samples in k-space
Nky = header.encoding.encodedSpace.matrixSize.y; % number of phase encodes in k-space
Nkz = header.encoding.encodedSpace.matrixSize.z; % number of slice encodes in k-space
Nx  = header.encoding.reconSpace.matrixSize.x;   % number of samples in image-space (RO)
Ny  = header.encoding.reconSpace.matrixSize.y;   % number of samples in image-space (PE)
Nz  = header.encoding.reconSpace.matrixSize.z;   % number of samples in image-space (SL)
Nc  = header.acquisitionSystemInformation.receiverChannels;

encoded_resolution = encoded_fov ./ [Nkx Nky Nkz]; % [mm]
recon_resolution = recon_fov ./ [Nx Ny Nz]; % [mm]

if isfield(user_opts, 'remove_oversampling')
    remove_oversampling = user_opts.remove_oversampling;
else
    remove_oversampling = 1;
end

if remove_oversampling
    osf = encoded_fov(1) / recon_fov(1); % oversampling factor in the x direction
else
    osf = 1;
end

zpad_matrix_size(1) = floor(encoded_fov(1) / osf / recon_resolution(1));
zpad_matrix_size(2) = round(encoded_fov(2) / recon_resolution(2));
zpad_matrix_size(3) = round(encoded_fov(3) / recon_resolution(3));

%--------------------------------------------------------------------------
% Calculate reconstruction parameters
%--------------------------------------------------------------------------
Nk = Nkx / osf; % number of readout samples
N1 = Nk;
N2 = Nky;
N3 = Nkz;
N = N1 * N2 * N3;

%% Parse the ISMRMRD header
tic; fprintf('Parsing the ISMRMRD header... ');
raw_data = dset.readAcquisition(); % read all the acquisitions
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%--------------------------------------------------------------------------
% Get data type
%--------------------------------------------------------------------------
is_noise                            = raw_data.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
is_parallel_calibration             = raw_data.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION');
is_parallel_calibration_and_imaging = raw_data.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING');
is_reverse                          = raw_data.head.flagIsSet('ACQ_IS_REVERSE');
is_navigation                       = raw_data.head.flagIsSet('ACQ_IS_NAVIGATION_DATA');
is_phasecorr                        = raw_data.head.flagIsSet('ACQ_IS_PHASECORR_DATA');
is_hpfeedback                       = raw_data.head.flagIsSet('ACQ_IS_HPFEEDBACK_DATA');
is_dummyscan                        = raw_data.head.flagIsSet('ACQ_IS_DUMMYSCAN_DATA');
is_rtfeedback                       = raw_data.head.flagIsSet('ACQ_IS_RTFEEDBACK_DATA');
is_surfacecoilcorrectionscan        = raw_data.head.flagIsSet('ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA');

%--------------------------------------------------------------------------
% ISMRMRD header
%--------------------------------------------------------------------------
% uint16_t version;                                    /**< First unsigned int indicates the version */
% uint64_t flags;                                      /**< bit field with flags */
% uint32_t measurement_uid;                            /**< Unique ID for the measurement */
% uint32_t scan_counter;                               /**< Current acquisition number in the measurement */
% uint32_t acquisition_time_stamp;                     /**< Acquisition clock */
% uint32_t physiology_time_stamp[ISMRMRD_PHYS_STAMPS]; /**< Physiology time stamps, e.g. ecg, breating, etc. */
% uint16_t number_of_samples;                          /**< Number of samples acquired */
% uint16_t available_channels;                         /**< Available coils */
% uint16_t active_channels;                            /**< Active coils on current acquisiton */
% uint64_t channel_mask[ISMRMRD_CHANNEL_MASKS];        /**< Mask to indicate which channels are active. Support for 1024 channels */
% uint16_t discard_pre;                                /**< Samples to be discarded at the beginning of  acquisition */
% uint16_t discard_post;                               /**< Samples to be discarded at the end of acquisition */
% uint16_t center_sample;                              /**< Sample at the center of k-space */
% uint16_t encoding_space_ref;                         /**< Reference to an encoding space, typically only one per acquisition */
% uint16_t trajectory_dimensions;                      /**< Indicates the dimensionality of the trajectory vector (0 means no trajectory) */
% float sample_time_us;                                /**< Time between samples in micro seconds, sampling BW */
% float position[3];                                   /**< Three-dimensional spatial offsets from isocenter */
% float read_dir[3];                                   /**< Directional cosines of the readout/frequency encoding */
% float phase_dir[3];                                  /**< Directional cosines of the phase */
% float slice_dir[3];                                  /**< Directional cosines of the slice direction */
% float patient_table_position[3];                     /**< Patient table off-center */
% ISMRMRD_EncodingCounters idx;                        /**< Encoding loop counters, see above */
% int32_t user_int[ISMRMRD_USER_INTS];                 /**< Free user parameters */
% float user_float[ISMRMRD_USER_FLOATS];               /**< Free user parameters */
%--------------------------------------------------------------------------
% Where EncodingCounters are defined as:
% uint16_t kspace_encode_step_1;    /**< e.g. phase encoding line number */
% uint16_t kspace_encode_step_2;    /**< e.g. partition encoding number */
% uint16_t average;                 /**< e.g. signal average number */
% uint16_t slice;                   /**< e.g. imaging slice number */
% uint16_t contrast;                /**< e.g. echo number in multi-echo */
% uint16_t phase;                   /**< e.g. cardiac phase number */
% uint16_t repetition;              /**< e.g. dynamic number for dynamic scanning */
% uint16_t set;                     /**< e.g. flow encoding set */
% uint16_t segment;                 /**< e.g. segment number for segmented acquisition */
% uint16_t user[ISMRMRD_USER_INTS]; /**< Free user parameters */
%--------------------------------------------------------------------------
number_of_samples  = double(max(raw_data.head.number_of_samples(~is_noise)));
discard_pre        = double(max(raw_data.head.discard_pre));
discard_post       = double(max(raw_data.head.discard_post));
center_sample      = double(max(raw_data.head.center_sample));
nr_channels        = double(max(raw_data.head.active_channels));
nr_phase_encodings = double(max(raw_data.head.idx.kspace_encode_step_1)) + 1; % nr_interleaves for spiral imaging
nr_slice_encodings = double(max(raw_data.head.idx.kspace_encode_step_2)) + 1;
nr_averages        = double(max(raw_data.head.idx.average)) + 1;
nr_slices          = double(max(raw_data.head.idx.slice)) + 1;
nr_contrasts       = double(max(raw_data.head.idx.contrast)) + 1;
nr_phases          = double(max(raw_data.head.idx.phase)) + 1;
nr_repetitions     = double(max(raw_data.head.idx.repetition)) + 1;
nr_sets            = double(max(raw_data.head.idx.set)) + 1;
nr_segments        = double(max(raw_data.head.idx.segment)) + 1; % echo train length or turbo factor
nr_samples         = number_of_samples - discard_pre - discard_post;
nr_shots           = nr_phase_encodings / nr_segments; % number of shots

%--------------------------------------------------------------------------
% Get the dimensionality of the trajectory vector (0 means no trajectory)
%--------------------------------------------------------------------------
trajectory_dimensions = double(max(raw_data.head.trajectory_dimensions));

%--------------------------------------------------------------------------
% Get the dwell time in [sec]
%--------------------------------------------------------------------------
dt = double(max(raw_data.head.sample_time_us)) * 1e-6; % [usec] * [sec/1e-6 usec] => [sec]

%--------------------------------------------------------------------------
% Calculate the readout duration [sec]
%--------------------------------------------------------------------------
T = nr_samples * dt; % readout duration [sec]

%% Display ISMRMRD header
fprintf('========================= ISMRMRD header ========================\n');
fprintf('encoded_fov        = %8.4f %8.4f %8.4f\n', encoded_fov(1), encoded_fov(2), encoded_fov(3));
fprintf('Nkx Nky Nkz        = %d      %d        %d\n', Nkx, Nky, Nkz);
fprintf('encoded_resolution = %8.4f %8.4f %8.4f\n', encoded_resolution(1), encoded_resolution(2), encoded_resolution(3));
fprintf('-----------------------------------------------------------------\n');
fprintf('recon_fov          = %8.4f %8.4f %8.4f\n', recon_fov(1), recon_fov(2), recon_fov(3));
fprintf('Nx Ny Nz           = %d      %d        %d\n', Nx, Ny, Nz);
fprintf('recon_resolution   = %8.4f %8.4f %8.4f\n', recon_resolution(1), recon_resolution(2), recon_resolution(3));
fprintf('-----------------------------------------------------------------\n');
fprintf('trajectory         = %s\n', header.encoding.trajectory);
fprintf('number_of_samples  = %d\n', number_of_samples);
fprintf('discard_pre        = %d\n', discard_pre);
fprintf('discard_post       = %d\n', discard_post);
fprintf('center_sample      = %d\n', center_sample);
fprintf('nr_channels        = %d\n', nr_channels);
fprintf('nr_phase_encodings = %d\n', nr_phase_encodings);
fprintf('nr_slice_encodings = %d\n', nr_slice_encodings);
fprintf('nr_averages        = %d\n', nr_averages);
fprintf('nr_slices          = %d\n', nr_slices);
fprintf('nr_contrasts       = %d\n', nr_contrasts);
fprintf('nr_phases          = %d\n', nr_phases);
fprintf('nr_repetitions     = %d\n', nr_repetitions);
fprintf('nr_sets            = %d\n', nr_sets);
fprintf('nr_segments        = %d\n', nr_segments); % echo train length!
fprintf('nr_shots           = %d\n', nr_shots);
fprintf('dt                 = %5.2f [usec]\n', dt * 1e6);
fprintf('readout duration   = %5.2f [msec]\n', T * 1e3);
fprintf('=================================================================\n');

%% Calculate the receiver noise matrix
[Psi,inv_L] = calculate_noise_decorrelation_matrix(ismrmrd_noise_fullpath);

%% Read a Siemens .dat file
fprintf('Reading a Siemens .dat file: %s\n', siemens_dat_fullpath);
twix = mapVBVD(siemens_dat_fullpath);
if length(twix) > 1
    twix = twix{end};
end

%% Reconstruct images per slice
nr_recons = nr_slices * nr_contrasts * nr_phases * nr_repetitions * nr_sets;

for idx = 1:nr_recons
    %% Get information about the current slice
    [slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr] = ind2sub([nr_slices nr_contrasts nr_phases nr_repetitions nr_sets], idx);
    fprintf('(%2d/%2d): Reconstructing slice (slice = %2d, contrast = %2d, phase = %2d, repetition = %2d, set = %2d)\n', idx, nr_recons, slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr);

    %% Get a list of profiles for all segments
    profile_list = find(~is_noise & ~is_phasecorr & ...
                        (raw_data.head.idx.slice      == (slice_nr - 1))      & ...
                        (raw_data.head.idx.contrast   == (contrast_nr - 1))   & ...
                        (raw_data.head.idx.phase      == (phase_nr - 1))      & ...
                        (raw_data.head.idx.repetition == (repetition_nr - 1)) & ...
                        (raw_data.head.idx.set        == (set_nr - 1)));

    %% Calculate the actual slice number for Siemens interleaved multislice imaging
    if isfield(twix.hdr.MeasYaps.sSliceArray, 'alSliceAcqOrder')
        alSliceAcqOrder = twix.hdr.MeasYaps.sSliceArray.alSliceAcqOrder{slice_nr};
        if isempty(alSliceAcqOrder), alSliceAcqOrder = 0; end
        actual_slice_nr = nr_slices - alSliceAcqOrder;
    else
        actual_slice_nr = slice_nr;
    end

    %% Get a slice normal vector from Siemens TWIX format
    %----------------------------------------------------------------------
    % dNormalSag: Sagittal component of a slice normal vector (in PCS)
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal, 'dSag')
        dNormalSag = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal.dSag;
    else
        dNormalSag = 0;
    end

    %----------------------------------------------------------------------
    % dNormalCor: Coronal component of a slice normal vector (in PCS)
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal, 'dCor')
        dNormalCor = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal.dCor;
    else
        dNormalCor = 0;
    end

    %----------------------------------------------------------------------
    % dNormalTra: Transverse component of a slice normal vector (in PCS)
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal, 'dTra')
        dNormalTra = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal.dTra;
    else
        dNormalTra = 0;
    end

    %----------------------------------------------------------------------
    % dRotAngle: Slice rotation angle ("swap Fre/Pha")
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}, 'dInPlaneRot')
        dRotAngle = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.dInPlaneRot; % [rad]
    else
        dRotAngle = 0; % [rad]
    end

    %% Determine the main orientation of an imaging stack
    main_orientation = fGSLClassOri(dNormalSag, dNormalCor, dNormalTra);

    %% Get a slice offset in the PCS from Siemens TWIX format
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}, 'sPosition')
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition, 'dSag')
            sag_offset = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition.dSag; % [mm]
        else
            sag_offset = 0; % [mm]
        end
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition, 'dCor')
            cor_offset = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition.dCor; % [mm]
        else
            cor_offset = 0; % [mm]
        end
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition, 'dTra')
            tra_offset = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition.dTra; % [mm]
        else
            tra_offset = 0; % [mm]
        end
    else
        sag_offset = 0; % [mm]
        cor_offset = 0; % [mm]
        tra_offset = 0; % [mm]
    end
    pcs_offset = [sag_offset; cor_offset; tra_offset] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Get a slice offset in the PCS from ISMRMRD format
    %----------------------------------------------------------------------
    % Get a list of profiles in the current slice
    %----------------------------------------------------------------------
    sag_offset_ismrmrd = double(raw_data.head.position(1,profile_list(1))); % [mm]
    cor_offset_ismrmrd = double(raw_data.head.position(2,profile_list(1))); % [mm]
    tra_offset_ismrmrd = double(raw_data.head.position(3,profile_list(1))); % [mm]
    pcs_offset_ismrmrd = [sag_offset_ismrmrd; cor_offset_ismrmrd; tra_offset_ismrmrd] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Get a rotation matrix from the GCS to the PCS (ISMRMRD format)
    phase_dir = double(raw_data.head.phase_dir(:,profile_list(1)));
    read_dir  = double(raw_data.head.read_dir(:,profile_list(1)));
    slice_dir = double(raw_data.head.slice_dir(:,profile_list(1)));
    R_gcs2pcs_ismrmrd = [phase_dir read_dir slice_dir];

    %% Calculate coordinate transform matrices
    %----------------------------------------------------------------------
    % Calculate a scaling matrix
    %----------------------------------------------------------------------
    scaling_matrix = diag(encoded_resolution) * 1e-3; % [mm] * [m/1e3mm] => [m]

    %----------------------------------------------------------------------
    % Calculate a transformation matrix from the RCS to the GCS [r,c,s] <=> [PE,RO,SL]
    %----------------------------------------------------------------------
    R_rcs2gcs = [0    1    0 ; % [PE]   [0 1 0] * [r]
                 1    0    0 ; % [RO] = [1 0 0] * [c]
                 0    0    1]; % [SL]   [0 0 1] * [s]

    %----------------------------------------------------------------------
    % Calculate a rotation matrix from the GCS to the PCS
    %----------------------------------------------------------------------
    [R_gcs2pcs,phase_sign,read_sign] = siemens_calculate_matrix_gcs_to_pcs(dNormalSag, dNormalCor, dNormalTra, dRotAngle);

    %----------------------------------------------------------------------
    % Calculate a rotation matrix from the PCS to the DCS
    %----------------------------------------------------------------------
    R_pcs2dcs = siemens_calculate_matrix_pcs_to_dcs(patient_position);

    %----------------------------------------------------------------------
    % Calculate a rotation matrix from the GCS to the DCS
    %----------------------------------------------------------------------
    R_gcs2dcs = R_pcs2dcs * R_gcs2pcs;

    %----------------------------------------------------------------------
    % Calculate a slice offset in the DCS [m]
    %----------------------------------------------------------------------
    dcs_offset = R_pcs2dcs * pcs_offset; % 3 x 1

    %% Display slice information
    fprintf('======================= SLICE INFORMATION =======================\n');
    fprintf('main_orientation = %d (SAGITTAL/CORONAL/TRANSVERSE = 0/1/2)\n', main_orientation);
    fprintf('actual_slice_nr = %d, slice_nr = %d\n', actual_slice_nr, slice_nr);
    fprintf('dNormalSag = %+g \ndNormalCor = %+g \ndNormalTra = %+g \ndRotAngle = %g [rad]\n', dNormalSag, dNormalCor, dNormalTra, dRotAngle);
    fprintf('phase_sign = %+g, read_sign = %+g\n', phase_sign, read_sign);
    fprintf('---------------------- From Siemens TWIX format ------------------\n');
    fprintf('                   [sag]   %10.5f [mm]\n', sag_offset);
    fprintf('slice offset     : [cor] = %10.5f [mm]\n', cor_offset);
    fprintf('                   [tra]   %10.5f [mm]\n', tra_offset);
    fprintf('---------------------- From ISMRMRD format -----------------------\n');
    fprintf('                   [sag]   %10.5f [mm]\n', sag_offset_ismrmrd);
    fprintf('slice offset     : [cor] = %10.5f [mm]\n', cor_offset_ismrmrd);
    fprintf('                   [tra]   %10.5f [mm]\n', tra_offset_ismrmrd);
    fprintf('---------------------- From Siemens TWIX format ------------------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2pcs(1,1), R_gcs2pcs(1,2), R_gcs2pcs(1,3));
    fprintf('R_gcs2pcs        : [cor] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2pcs(2,1), R_gcs2pcs(2,2), R_gcs2pcs(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2pcs(3,1), R_gcs2pcs(3,2), R_gcs2pcs(3,3));
    fprintf('---------------------- From ISMRMRD format (incorrect!)-----------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2pcs_ismrmrd(1,1), R_gcs2pcs_ismrmrd(1,2), R_gcs2pcs_ismrmrd(1,3));
    fprintf('R_gcs2pcs        : [cor] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2pcs_ismrmrd(2,1), R_gcs2pcs_ismrmrd(2,2), R_gcs2pcs_ismrmrd(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2pcs_ismrmrd(3,1), R_gcs2pcs_ismrmrd(3,2), R_gcs2pcs_ismrmrd(3,3));
    fprintf('------------------------------------------------------------------\n');
    fprintf('                   [ x ]   [%10.5f %10.5f %10.5f][sag]\n', R_pcs2dcs(1,1), R_pcs2dcs(1,2), R_pcs2dcs(1,3));
    fprintf('R_pcs2dcs        : [ y ] = [%10.5f %10.5f %10.5f][cor]\n', R_pcs2dcs(2,1), R_pcs2dcs(2,2), R_pcs2dcs(2,3));
    fprintf('                   [ z ]   [%10.5f %10.5f %10.5f][tra]\n', R_pcs2dcs(3,1), R_pcs2dcs(3,2), R_pcs2dcs(3,3));
    fprintf('------------------------------------------------------------------\n');
    fprintf('                   [ x ]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2dcs(1,1), R_gcs2dcs(1,2), R_gcs2dcs(1,3));
    fprintf('R_gcs2dcs        : [ y ] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2dcs(2,1), R_gcs2dcs(2,2), R_gcs2dcs(2,3));
    fprintf('                   [ z ]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2dcs(3,1), R_gcs2dcs(3,2), R_gcs2dcs(3,3));
    fprintf('=================================================================\n');

    %% Get information about k-space encoding
    kspace_encoding_step_1_center = header.encoding.encodingLimits.kspace_encoding_step_1.center;
    kspace_encoding_step_2_center = header.encoding.encodingLimits.kspace_encoding_step_2.center;
    kspace_encoding_step_1_maximum = header.encoding.encodingLimits.kspace_encoding_step_1.maximum;
    kspace_encoding_step_2_maximum = header.encoding.encodingLimits.kspace_encoding_step_2.maximum;

    %% Calculate the phase encode view ordering (Nkx x nr_segments x nr_shots)
    %----------------------------------------------------------------------
    % Calculate the index range (1st index) of a profile
    %----------------------------------------------------------------------
    k1_range = (discard_pre+1:number_of_samples-discard_post).' - (center_sample + 1) + floor(Nkx/2) + 1;

    scan_counter = double(raw_data.head.scan_counter(profile_list));
    scan_counter_min = min(scan_counter);

    view_ordering = zeros(Nky, nr_segments, nr_shots, 'double');
    for idx2 = 1:nr_shots
        scan_counter_offset = scan_counter_min + (idx2 - 1) * nr_segments;
        for idx1 = 1:nr_segments % ETL or turbo factor
            profile_nr = scan_counter_offset + (idx1 - 1);
            tstart = tic; fprintf('Calculating the view ordering (echo=%2d/%2d): shot=%2d, profile=%3d... ', idx1, nr_segments, idx2, profile_nr);

            %--------------------------------------------------------------
            % Calculate the (2nd,3rd) matrix index of a profile
            %--------------------------------------------------------------
            kspace_encode_step_1 = double(raw_data.head.idx.kspace_encode_step_1(profile_nr));
            kspace_encode_step_2 = double(raw_data.head.idx.kspace_encode_step_2(profile_nr));
            k2_index = kspace_encode_step_1 - (kspace_encoding_step_1_maximum - kspace_encoding_step_1_center + 1) + floor(Nky/2) + 1;
            k3_index = kspace_encode_step_2 - (kspace_encoding_step_2_maximum - kspace_encoding_step_2_center + 1) + floor(Nkz/2) + 1;
            if k3_index == 0, k3_index = k3_index + 1; end % For 2D imaging
            view_ordering(k2_index,idx1,idx2) = 1;
            fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
        end
    end

    %% Get raw k-space data per echo time (Nkx x Nky x Nkz x Nc x nr_segments)
    sampling_pattern = zeros(Nky, Nkz, nr_segments, 'double');
    kspace = complex(zeros(Nkx, Nky, Nkz, Nc, nr_segments, 'double'));
    for segment_nr = 1:nr_segments
        tstart = tic; fprintf('Reading k-space data (%2d/%2d)... ', segment_nr, nr_segments);
        %------------------------------------------------------------------
        % Get a list of profiles per segment
        %------------------------------------------------------------------
        profile_list_segment = find(~is_noise & ~is_phasecorr & ...
                                    (raw_data.head.idx.slice      == (slice_nr - 1))      & ...
                                    (raw_data.head.idx.contrast   == (contrast_nr - 1))   & ...
                                    (raw_data.head.idx.phase      == (phase_nr - 1))      & ...
                                    (raw_data.head.idx.repetition == (repetition_nr - 1)) & ...
                                    (raw_data.head.idx.set        == (set_nr - 1))        & ...
                                    (raw_data.head.idx.segment    == (segment_nr - 1)));

        %------------------------------------------------------------------
        % Calculate a sampling pattern (Nky x Nkz)
        %------------------------------------------------------------------
        for idx1 = 1:nr_shots
            kspace_encode_step_1 = double(raw_data.head.idx.kspace_encode_step_1(profile_list_segment(idx1)));
            kspace_encode_step_2 = double(raw_data.head.idx.kspace_encode_step_2(profile_list_segment(idx1)));
            k2_index = kspace_encode_step_1 - (kspace_encoding_step_1_maximum - kspace_encoding_step_1_center + 1) + floor(Nky/2) + 1;
            k3_index = kspace_encode_step_2 - (kspace_encoding_step_2_maximum - kspace_encoding_step_2_center + 1) + floor(Nkz/2) + 1;
            if k3_index == 0, k3_index = k3_index + 1; end % For 2D imaging
            sampling_pattern(k2_index,k3_index,segment_nr) = sampling_pattern(k2_index,k3_index,segment_nr) + 1;

            %--------------------------------------------------------------
            % Prewhiten k-space data
            %--------------------------------------------------------------
            profile = raw_data.data{profile_list_segment(idx1)}; % nr_samples x Nc
            profile = (inv_L * profile.').';

            %--------------------------------------------------------------
            % Accumulate k-space
            %--------------------------------------------------------------
            kspace(k1_range,k2_index,k3_index,:,segment_nr) = kspace(k1_range,k2_index,k3_index,:,segment_nr) + reshape(profile, [nr_samples 1 1 Nc]);
        end
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Flip k-space
    tstart = tic; fprintf('Flipping k-space... ');
    if read_sign == -1
        kspace = flip(kspace,1);
    end
    if phase_sign == -1
        kspace = flip(kspace,2);
        sampling_pattern = flip(sampling_pattern,1);
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Remove oversampling: (Nkx x Nky x Nkz x Nc x nr_segments) => (Nk x Nky x Nkz x Nc x nr_segments)
    %----------------------------------------------------------------------
    % Reconstruct in x (k-space <=> image-space)
    %----------------------------------------------------------------------
    if osf > 1
        tstart = tic; fprintf('Removing oversampling along the readout direction... ');
        projection = 1 / sqrt(Nkx) * fftshift(fft(ifftshift(kspace, 1), [], 1), 1);
        idx1_range = (-floor(Nk/2):ceil(Nk/2)-1).' + floor(Nkx/2) + 1;
        kspace = sqrt(Nk / osf) * fftshift(ifft(ifftshift(projection(idx1_range,:,:,:,:), 1), [], 1), 1);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
        clear projection;
    end
end
