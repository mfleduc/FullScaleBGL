clear variables;close all;
%% We are investigating the temporal correlation in the data 
data = load('/home/male7736/Desktop/Research/FullScaleBGL/data/6 hour time resolution neutral temp.mat');
[LON,LAT] = ndgrid(data.lon,data.lat);
SLT = LON/15;
detrended = data.data-mean(data.data);
ndcs = randi( numel(SLT), 16, 1 );
figure;
for ii=1:16
    [acfs(ii,:) , lag] =  xcov(detrended(:,ndcs(ii)),'coeff');
    subplot(4,4,ii)
    plot(lag*6,acfs(ii,:))
    xlabel("Lag, hours");ylabel("Autocorrelation")
    grid on
    title(sprintf("Latitude %.1f deg, SLT %.1f hrs",LAT(ndcs(ii)),SLT(ndcs(ii))))
    ylim([-0.5,1])
end
sgtitle("Sample temporal autocorrelation functions at different locations")