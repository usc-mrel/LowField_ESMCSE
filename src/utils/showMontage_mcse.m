function [outputArg1,outputArg2] = showMontage_mcse(im,nr_rows,nr_cols,str_title,colorrange)

[N1,N2,N3,M] = size(im);
im_montage = complex(zeros(N2*nr_rows, N1*nr_cols, 'double'));
for m = 1:M
    % 1   2  3  4  5
    % 6   7  8  9 10
    % 11 12 13 14 15
    col = mod(m-1,nr_cols) + 1;
    row = floor((m - 1) / nr_cols) + 1;
    idx1_range = (1:N2).' + (row - 1) * N2;
    idx2_range = (1:N1).' + (col - 1) * N1;
    im_montage(idx1_range,idx2_range) = im(:,:,m);
end


figure;
imagesc(abs(im_montage),colorrange);
axis image;
colormap(gca, gray(256));
%colorbar;
title(str_title)

figure('Color', 'w', 'Position', [8 110 1580 706]);
imagesc(angle(im_montage)*180/pi);
axis image;
colormap(gca, hsv(256));
colorbar;
title(str_title)
end



