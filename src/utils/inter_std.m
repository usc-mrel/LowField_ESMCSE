function [sigma_b] = inter_std(data1,data2,data3,K1,K2,K3)

if (numel(data1) ~= K1) && (numel(data1) ~= K3) && (numel(data3) ~= K3)
    exit
end
mean_all = mean([data1,data2,data3]);

bms      = (K1*(mean(data1)-mean_all)^2 + K2*(mean(data2)-mean_all)^2 + K3*(mean(data3)-mean_all)^2)/3;

wms      = (sum((data1-mean(data1)).^2) + sum((data2-mean(data2)).^2) + sum((data3-mean(data3)).^2)) / (K1+K2+K3-3);

sigma_b  = sqrt((bms-wms)/min([K1,K2,K3]));
end