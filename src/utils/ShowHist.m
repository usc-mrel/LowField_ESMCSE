function ShowHist(input1,bcounts,BLm,xText,rm)
Fsize = 20;
figure,
histogram(input1,bcounts,'BinLimits',BLm);
hold on;
if rm == 1
    histogram(rmoutliers(input1,'percentiles',[25 75]),bcounts,'BinLimits',BLm);
    xlabel(xText,'FontSize',Fsize)
end
set(gca,'FontSize',Fsize)
end

