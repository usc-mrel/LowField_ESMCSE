video_file = VideoWriter('/Users/bli/MREL/LungT2_Paper/Figures/S6/Mag.mp4', 'MPEG-4');
video_file.FrameRate = 10; % need to adjust this frame rate
video_file.Quality = 100;
open(video_file);


echo_spacing = 20 * 1e-3; % [sec]
t_se = (1+(1:8)).' * echo_spacing;
t_shift     =   repmat([-4,-2,0,2,4].',[8,1]);        % [ms]
t_SE        =   reshape(repmat(t_se*1e3,[1,5]).',[],1); % [ms]
figure,
for idx = 6:40
    tseVal      = t_SE(idx);
    tshiftVal   = t_shift(idx);
    imagesc( abs(squeeze(im_fft_signal(:,:,:,idx))) ./ max(abs(squeeze(im_fft_signal(:,:,:,8))),[],"all"),[0 0.5]); % Display an image
    colormap("gray")
    axis image off
    title(sprintf('t_{SE} = %4.1f [ms], t_{shift} = %4.1f [ms]',tseVal,tshiftVal))
    set(gca,'FontSize',18)
    writeVideo(video_file, getframe(gcf));
end
close(video_file);