close all
clear all
clc

%% Parameters
plot                 = 0; 
FSR                  = 0.0165; % Free Spectral Range 0.0165 cm-1
MPH                  = 0.03;
sel                  = 0.03;

%% Main
raw_etalon_data      = csvread('EtalonData.csv');

etalon_data          = raw_etalon_data(944:4691); 
midpoint             = round((4691-944)/2);

upscan_data          = etalon_data(midpoint:length(etalon_data));
downscan_data        = etalon_data(1:midpoint-1);

[I1,y1]              = peakfinder(upscan_data, sel, MPH, 1);

for i=1:length(I1)
    peaks(i,1)       = I1(i,1);
    peaks(i,2)       = y1(i,1);
    peaks(i,3)       = -FSR*i;
end

C_upscan             = polyfit(peaks(1:end,1),peaks(1:end,3),6);
 
etalon_fit_upscan    = polyval(C_upscan,peaks(:,1));

[I2,y2]              = peakfinder(downscan_data, sel, MPH, 1);
I2                   = I2(3:(length(I2)-1)); 
I2(5)                = [];
I2(3)                = [];% remove bad point
y2                   = y2(3:(length(y2)-1)); 
y2(5)                = [];
y2(3)                = []; % remove bad point

if plot
plot(downscan_data);hold on;plot(I2,y2,'r*')
%plot(etalon_data);hold on;plot(I1+midpoint,y1,'b*');plot(I2,y2,'r*')
end

for i=1:length(I2)
    peaks1(i,1)      = I2(i,1);
    peaks1(i,2)      = y2(i,1);
    peaks1(i,3)      = -0.891+FSR*i;
end

C_downscan           = polyfit(peaks1(1:end,1),peaks1(1:end,3),6); 
etalon_fit_downscan  = polyval(C_downscan,peaks1(:,1));

%%
if plot 
figure(1)
plot(upscan_data,'b');hold on;
plot(I1,y1,'r*')
xlabel('Index Number')
ylabel('Detector Voltage, V')

figure(2)
subplot(2,1,1)
plot(upscan_data,'b');hold on;
plot(peaks(:,1),peaks(:,2),'r*')
title('Etalon Upscan and Relative Wavenumber')
hold off
xlabel('Index Number')
ylabel('Detector Voltage, V')
subplot(2,1,2)
plot(peaks(:,1),peaks(:,3),'b*');hold on;
plot(peaks(:,1),etalon_fit_upscan,'r')
xlabel('Index Number')
ylabel('Wavenumber')

figure(3)
plot(downscan_data,'b');hold on;
plot(I2,y2,'r*')
xlabel('Index Number')
ylabel('Detector Voltage, V')

figure(4)
subplot(2,1,1)
plot(downscan_data,'b');hold on;
plot(peaks1(:,1),peaks1(:,2),'r*')
hold off
title('Etalon Downscan and Relative Wavenumber')
xlabel('Index Number')
ylabel('Detector Voltage, V')
subplot(2,1,2)
plot(peaks1(:,1),peaks1(:,3),'b*');hold on;
plot(peaks1(:,1),etalon_fit_downscan,'r')
hold off
xlabel('Index Number')
ylabel('Wavenumber')
end
