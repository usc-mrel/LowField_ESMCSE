% demo_estse_dictionary_simulation_R2.m
% Written by Bochao Li
% Email: bochaoli@usc.edu
% Started: 06/08/2022, Last modified: 08/26/2022

%% Clean slate
%% Define parameters
echo_spacing = 20 * 1e-3; % [sec]
ETL          = 8;         % echo train length
sim_times    = 10000;      % number of simulation for each R2 and R2' case
SNR          = 18 ;       % SNR in lung region (from SNR recon)

R2p_array = (60:5:110);    % [Hz]
R2_array  = (14:1:20);     % [Hz]

%R2_array    = (20);     %   [Hz]
%R2p_array   = (60);     %   [Hz]


df_array    = (-30:5:30);       %   [Hz]
%df_array = 10;
%% Calculate a continuous signal evolution
nr_shifts_continous = 1201;
t_shift_continous = (-floor(nr_shifts_continous/2):ceil(nr_shifts_continous/2)-1).' * echo_spacing / nr_shifts_continous; % [sec]
t_se_continous = (1:ETL).' * echo_spacing;

t_shift = [-4; -2; 0; 2; 4] * 1e-3; % [sec]
t_se = (1+(1:7)).' * echo_spacing;

neg_shift_index_continous = find(t_shift_continous < 0);
pos_shift_index_continous = find(t_shift_continous >= 0);

neg_shift_continous = t_shift_continous(neg_shift_index_continous);
pos_shift_continous = t_shift_continous(pos_shift_index_continous);

acquired_samples_noisy  = zeros(length(t_shift) * length(t_se), sim_times , length(R2p_array), length(R2_array),length(df_array), 'double');
acquired_samples   = zeros(length(t_shift) * length(t_se), length(R2p_array), length(R2_array),length(df_array), 'double');
signal_evolution   = zeros(nr_shifts_continous * ETL, length(R2p_array), length(R2_array),length(df_array), 'double');

%% Loop for all R2 and R2p
for ndf = 1:length(df_array)
    df = df_array(ndf);

    for nR2 = 1:length(R2_array)
        R2 = R2_array(nR2);

        for nR2p = 1:length(R2p_array)
            R2p = R2p_array(nR2p);

            %% Generate noiseless signal (continous)
            for idx = 1:ETL
                %----------------------------------------------------------------------
                % t_shift < 0
                %----------------------------------------------------------------------
                index_range1 = neg_shift_index_continous + nr_shifts_continous * (idx - 1);
                signal_evolution(index_range1,nR2p,nR2,ndf) = exp(-1j * 2 * pi * df * neg_shift_continous) .* exp(-R2 * (t_se_continous(idx) + neg_shift_continous)) .* exp(R2p * neg_shift_continous);
                %----------------------------------------------------------------------
                % t_shift > 0
                %----------------------------------------------------------------------
                index_range2 = pos_shift_index_continous + nr_shifts_continous * (idx - 1);
                signal_evolution(index_range2,nR2p,nR2,ndf) = exp(-1j * 2 * pi * df * pos_shift_continous) .* exp(-R2 * (t_se_continous(idx) + pos_shift_continous)) .* exp(-R2p * pos_shift_continous);
            end

            %% Calculate a time vector
            t_continuous = reshape(bsxfun(@plus, t_shift_continous, t_se_continous.'), [nr_shifts_continous * ETL 1]);

            %% Calculate actual acquired sample points
            nr_shifts = length(t_shift);
            neg_shift_index = find(t_shift < 0);
            pos_shift_index = find(t_shift >= 0);

            neg_shift = t_shift(neg_shift_index);
            pos_shift = t_shift(pos_shift_index);

            nr_echoes = length(t_se);

            for idx = 1:nr_echoes
                %----------------------------------------------------------------------
                % t_shift < 0
                %----------------------------------------------------------------------
                index_range1 = neg_shift_index + nr_shifts * (idx - 1);
                acquired_samples(index_range1,nR2p,nR2,ndf) = exp(-1j * 2 * pi * df * neg_shift) .* exp(-R2 * (t_se(idx) + neg_shift)) .* exp(R2p * neg_shift);

                %----------------------------------------------------------------------
                % t_shift > 0
                %----------------------------------------------------------------------
                index_range2 = pos_shift_index + nr_shifts * (idx - 1);
                acquired_samples(index_range2,nR2p,nR2,ndf) = exp(-1j * 2 * pi * df * pos_shift) .* exp(-R2 * (t_se(idx) + pos_shift)) .* exp(-R2p * pos_shift);
            end
        end
    end
end
%% Calculate noisy sample
S0      = max(acquired_samples(:));
sigma   = S0/SNR;
for ndf = 1:length(df_array)
    for nR2 = 1:length(R2_array)
        for nR2p = 1:length(R2p_array)
            for n_time = 1:sim_times
                R_rand = sigma * randn(35,1)/sqrt(2);
                I_rand = sigma * randn(35,1)/sqrt(2);
                %acquired_samples_noisy_real(:,n_time,nT2p,nT2) = abs(acquired_samples_real(:,nT2p,nT2) +...
                % (R_rand/sqrt(2) + 1j*I_rand/sqrt(2)));
                acquired_samples_noisy(:,n_time,nR2p,nR2,ndf) = acquired_samples(:,nR2p,nR2,ndf) + R_rand+1j*I_rand;
            end
        end
    end
end
%% Calculate a time vector
t_sampled = reshape(bsxfun(@plus, t_shift, t_se.'), [nr_shifts * nr_echoes 1]);

%% fitting the acquired_samples for paratmeter estimation
%Set source directories
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
    src_directory = '/server/home/nlee/GNL_lowrank_recon_3d_gpu';
    thirdparty_directory = 'E:\gradient_nonlinearity\thirdparty';
    ismrmrd_directory = '/server/home/nlee/ismrmrd';
end
%Add source directories to search path
addpath(genpath(src_directory));
addpath(genpath(thirdparty_directory));
addpath(genpath(ismrmrd_directory));

fprintf('Start noiseless fitting... ');
Necho  = 7;
Nshift = 5;
t_shift     =   repmat([-4,-2,0,2,4].',[Necho,1]);        % [ms]
t_SE        =   reshape(repmat(t_se*1e3,[1,Nshift]).',[],1); % [ms]
NR2p    = size(acquired_samples_noisy,3);
NR2     = size(acquired_samples_noisy,4);
Ndf     = size(acquired_samples_noisy,5);

tic
parfor ndf = 1:Ndf
    df = df_array(ndf)
    for nR2 = 1:NR2
        nR2
        for nR2p = 1:NR2p
            for ntime = 1:sim_times
                %y = [acquired_samples_noisy_real(:,ntime,nT2p,nT2);acquired_samples_noisy_img(:,ntime,nT2p,nT2) ];
                %x0 = [acquired_samples_noisy_real(3,ntime,nT2p,nT2),acquired_samples_noisy_img(3,ntime,nT2p,nT2),50, 10].'; % inital values [a.u., kHz, kHz,Hz]
                %E = @(x) dictionary_complex_mixed_rho_t2_t2prime_B0_fitting_esmcse(x, y, t_shift, t_SE);
                y =  [real(acquired_samples_noisy(:,ntime,nR2p,nR2,ndf));imag(acquired_samples_noisy(:,ntime,nR2p,nR2,ndf))];
                x0 = [real(acquired_samples_noisy(3,ntime,nR2p,nR2)),imag(acquired_samples_noisy(3,ntime,nR2p,nR2)),16/1000, 10/1000, 0].';
                %E = @(x) dictionary_mixed_rho_r2_r2prime_B0_fitting_esmcse(x, y, t_shift, t_SE);
                E = @(x) dictionary_seperate_mixed_rho_t2_t2prime_B0_fitting_esmcse(x, y, t_shift, t_SE);
                options = optimoptions('lsqnonlin','Display','none');
                [x_out_dictionary(:,ntime,nR2p,nR2,ndf),resnorm(ntime,nR2p,nR2,ndf),residual(:,ntime,nR2p,nR2,ndf),exitflag,output] = lsqnonlin(E,x0,[-Inf,-Inf,-Inf,-Inf,-Inf],[Inf,Inf,Inf,Inf,Inf],options);
            end
        end
    end
end
toc


%% analyze the error
R2_est          = squeeze(x_out_dictionary_sp(3,:,:,:,:)) * 1e3;  %   [Hz]
R2p_est         = squeeze(x_out_dictionary_sp(4,:,:,:,:)) * 1e3;  %   [Hz]
df_est          = squeeze(x_out_dictionary_sp(5,:,:,:,:)) * 1e3;  %   [Hz]

error_R2        = R2_est - permute(repmat(R2_array.',1,sim_times,NR2p,Ndf),[2,3,1,4]);
error_R2p       = R2p_est - permute(repmat(R2p_array.',1,sim_times,NR2,Ndf),[2,1,3,4]);
error_df        = df_est - permute(repmat(df_array.',1,sim_times,NR2p,NR2),[2,3,4,1]);


bias_R2   = squeeze(mean(error_R2,1));
bias_R2p  = squeeze(mean(error_R2p,1));
bias_df   = squeeze(mean(error_df,1));

std_R2  = squeeze(std(R2_est,[],1));
std_R2p = squeeze(std(R2p_est,[],1));
std_df  = squeeze(std(df_est,[],1));

%% Calculate fitted signal
signal_evolution_fitting = zeros(nr_shifts_continous * ETL,1, 'double');
nR2     = 7;
nR2p    = 11;
nSim    = 50;
ndf     = 11;
R2val   = R2_est(nSim,nR2p,nR2,ndf)
R2pval  = R2p_est(nSim,nR2p,nR2,ndf)
dfval   = df_est(nSim,nR2p,nR2,ndf)
for idx = 1:ETL
    %----------------------------------------------------------------------
    % t_shift < 0
    %----------------------------------------------------------------------
    index_range1 = neg_shift_index_continous + nr_shifts_continous * (idx - 1);
    signal_evolution_fitting(index_range1) = exp(-1j * 2 * pi * dfval * neg_shift_continous) .* exp(-R2val * (t_se_continous(idx) + neg_shift_continous)) .* exp(R2pval * neg_shift_continous);

    %----------------------------------------------------------------------
    % t_shift > 0
    %----------------------------------------------------------------------
    index_range2 = pos_shift_index_continous + nr_shifts_continous * (idx - 1);
    signal_evolution_fitting(index_range2) = exp(-1j * 2 * pi * dfval * pos_shift_continous) .* exp(-R2val * (t_se_continous(idx) + pos_shift_continous)) .* exp(-R2pval * pos_shift_continous);
end

%% Display a signal evolution
R2  = R2_array(nR2);
R2p = R2p_array(nR2p);
df  =df_array(ndf);
FontSize = 20;
MarkerSize = 25;
LineWidth = 2;
figure('Color', 'w', 'Position', [9 388 1071 410]);
hold on;
grid on;
color_order = get(gca, 'colororder');
plot(t_continuous*1e3, signal_evolution(:,nR2p,nR2,ndf), '-', 'LineWidth', LineWidth);
plot(t_continuous*1e3, real(signal_evolution_fitting), 'r--', 'LineWidth', LineWidth);
plot(t_sampled(1:5:end)*1e3, real(acquired_samples_noisy(1:5:end,nSim,nR2p,nR2,ndf)), '.', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'Color', color_order(2,:));
plot(t_sampled(2:5:end)*1e3, real(acquired_samples_noisy(2:5:end,nSim,nR2p,nR2,ndf)), '.', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'Color', color_order(3,:));
plot(t_sampled(3:5:end)*1e3, real(acquired_samples_noisy(3:5:end,nSim,nR2p,nR2,ndf)), '.', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'Color', color_order(4,:));
plot(t_sampled(4:5:end)*1e3, real(acquired_samples_noisy(4:5:end,nSim,nR2p,nR2,ndf)), '.', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'Color', color_order(5,:));
plot(t_sampled(5:5:end)*1e3, real(acquired_samples_noisy(5:5:end,nSim,nR2p,nR2,ndf)), '.', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'Color', color_order(7,:));
legend({'Simulated signal','fitted signal', sprintf('$t_{\\mathrm{shift}}$=%4.1f ms', t_shift(1)*1e3), ...
    sprintf('$t_{\\mathrm{shift}}$=%4.1f ms', t_shift(2)*1e3), ...
    sprintf('$t_{\\mathrm{shift}}$=%4.1f ms', t_shift(3)*1e3), ...
    sprintf('$t_{\\mathrm{shift}}$=%4.1f ms', t_shift(4)*1e3), ...
    sprintf('$t_{\\mathrm{shift}}$=%4.1f ms', t_shift(5)*1e3)}, 'Interpreter', 'latex', 'FontSize', FontSize);

text(t_se(1)*1e3, real(acquired_samples_noisy(3+5*0,nSim,nR2p,nR2,ndf))+0.05, sprintf('$t_{\\mathrm{SE}}$=%2.0f ms', t_se(1)*1e3), 'Interpreter', 'latex', 'FontSize', FontSize);
text(t_se(2)*1e3, real(acquired_samples_noisy(3+5*1,nSim,nR2p,nR2,ndf))+0.05, sprintf('%2.0f ms', t_se(2)*1e3), 'Interpreter', 'latex', 'FontSize', FontSize);
text(t_se(3)*1e3, real(acquired_samples_noisy(3+5*2,nSim,nR2p,nR2,ndf))+0.05, sprintf('%2.0f ms', t_se(3)*1e3), 'Interpreter', 'latex', 'FontSize', FontSize);
text(t_se(4)*1e3, real(acquired_samples_noisy(3+5*3,nSim,nR2p,nR2,ndf))+0.05, sprintf('%2.0f ms', t_se(4)*1e3), 'Interpreter', 'latex', 'FontSize', FontSize);
text(t_se(5)*1e3, real(acquired_samples_noisy(3+5*4,nSim,nR2p,nR2,ndf))+0.05, sprintf('%2.0f ms', t_se(5)*1e3), 'Interpreter', 'latex', 'FontSize', FontSize);
text(t_se(6)*1e3, real(acquired_samples_noisy(3+5*5,nSim,nR2p,nR2,ndf))+0.05, sprintf('%2.0f ms', t_se(6)*1e3), 'Interpreter', 'latex', 'FontSize', FontSize);
text(t_se(7)*1e3, real(acquired_samples_noisy(3+5*6,nSim,nR2p,nR2,ndf))+0.05, sprintf('%2.0f ms', t_se(7)*1e3), 'Interpreter', 'latex', 'FontSize', FontSize);
set(gca, 'Box', 'On', 'FontSize', FontSize);
xlabel('Time [ms]', 'FontSize', FontSize);
ylabel('S/S_0', 'FontSize', FontSize);
title({sprintf('ES-MCSE signal evolution, $\\mathrm{R}_2/\\mathrm{R}_2''$ = %4.2f/%4.2f [Hz]', R2, R2p)}, 'Interpreter', 'latex');
ylim([0 0.7]);
xlim([30 180])
set(gca,'FontSize',20)
%export_fig(sprintf('estse_signal_evolution_R2_%2.0f_R2p_%2.0f_subject%d', R2, R2p, subject_nr), '-r300', '-tif');

%% Show histograms
figure,histogram(bias_R2);title('Bias of R_2');set(gca,'FontSize',16)
figure,histogram(bias_R2p);title('Bias of R_2''');set(gca,'FontSize',16)
figure,histogram(bias_df);title('Bias of \Deltaf');set(gca,'FontSize',16)

%% Show R2 Bias and Std maps
ha = tight_subplot(3,6,[.1 .03],[.1 .05],[.03 .01]);
for k = 1:17
    axes(ha(k));
    imagesc(bias_R2(:,:,k),[-0.03 0.03]);
    title("\Deltaf = " + df_array(k) + "Hz");
    y = colorbar;title(y,'Hz'); 
    ylabel('R_2'' [Hz]'); xlabel('R_2 [Hz]');
    xticks([1:2:length(R2_array)]);xticklabels(num2str(R2_array(1:2:end).'));
    yticks([1:2:length(R2p_array)]);yticklabels(num2str(R2p_array(1:2:end).'));
    set(gca,'FontSize',16)
end

ha = tight_subplot(3,6,[.1 .03],[.1 .05],[.03 .01]);
for k = 1:17
    axes(ha(k));
    imagesc(std_R2(:,:,k),[0.4 1]);
    title("\Deltaf = " + df_array(k) + "Hz");
    y = colorbar;title(y,'Hz'); 
    ylabel('R_2'' [Hz]'); xlabel('R_2 [Hz]');
    xticks([1:2:length(R2_array)]);xticklabels(num2str(R2_array(1:2:end).'));
    yticks([1:2:length(R2p_array)]);yticklabels(num2str(R2p_array(1:2:end).'));
    set(gca,'FontSize',16)
end

%% Show R2prime Bias and Std maps
ha = tight_subplot(3,6,[.1 .03],[.1 .05],[.03 .01]);
for k = 1:17
    axes(ha(k));
    imagesc(bias_R2p(:,:,k),[-0.5 0.5]);
    title("\Deltaf = " + df_array(k) + "Hz");
    y = colorbar;title(y,'Hz'); 
    ylabel('R_2'' [Hz]'); xlabel('R_2 [Hz]');
    xticks([1:2:length(R2_array)]);xticklabels(num2str(R2_array(1:2:end).'));
    yticks([1:2:length(R2p_array)]);yticklabels(num2str(R2p_array(1:2:end).'));
    set(gca,'FontSize',16)
end

ha = tight_subplot(3,6,[.1 .03],[.1 .05],[.03 .01]);
for k = 1:17
    axes(ha(k));
    imagesc(std_R2p(:,:,k),[0 15]);
    title("\Deltaf = " + df_array(k) + "Hz");
    y = colorbar;title(y,'Hz'); 
    ylabel('R_2'' [Hz]'); xlabel('R_2 [Hz]');
    xticks([1:2:length(R2_array)]);xticklabels(num2str(R2_array(1:2:end).'));
    yticks([1:2:length(R2p_array)]);yticklabels(num2str(R2p_array(1:2:end).'));
    set(gca,'FontSize',16)
end

%%
for k = 5:17
    R2_est_30(10000*(k-5)+1:10000*(k-5)+10000,:,:) =  R2_est(:,:,:,k);
    R2p_est_30(10000*(k-5)+1:10000*(k-5)+10000,:,:) =  R2p_est(:,:,:,k);
end

std_R2_est_30 = squeeze(std(R2_est_30,[],1));
std_R2p_est_30 = squeeze(std(R2p_est_30,[],1));

figure,imagesc(std_R2_est_30,[0,1]);y = colorbar;title(y,'Hz','FontSize',18);
xticks([1:2:length(R2_array)]);xticklabels(num2str(R2_array(1:2:end).'));
yticks([1:2:length(R2p_array)]);yticklabels(num2str(R2p_array(1:2:end).'));
ylabel('R_2'' [Hz]');xlabel('R_2 [Hz]');set(gca,'FontSize',18);
colormap(turbo);axis image

figure,imagesc(std_R2p_est_30,[0 14]);y = colorbar;title(y,'Hz','FontSize',18);
xticks([1:2:length(R2_array)]);xticklabels(num2str(R2_array(1:2:end).'));
yticks([1:2:length(R2p_array)]);yticklabels(num2str(R2p_array(1:2:end).'));
ylabel('R_2'' [Hz]');xlabel('R_2 [Hz]');set(gca,'FontSize',18);
colormap(turbo);axis image
