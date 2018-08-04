% Laser temperature 7.75C
% Laser current 117 mA
% Triangle Wave Voltage 0-600 mV
% Triangle Wave Frequency Voltage 500 Hz
% Sample Rate for detector signal 1.875 MHz
close all
clear all
clc

global  vd Absorbance_Fit

%% Read data file
% Read steady flame data or unsteady flame data

raw_data         = csvread('SteadyFlameData.csv');
% raw_data         = csvread('UnsteadyFlameData.csv');

%% Parameters
SR                      = 1.875e6; % Sample rate
Scan                    = 500;     % Scan frequency 500 Hz
Tguess                  = 1500; 
vo                      = 2060;
MW                      = 28;      % Molecular weight of carbon monoxide
points_per_scan         = 1/Scan*SR;
L                       = 1.25;    % 1.25 cm path length
save                    = 0;
load                    = 0;
plot                    = 0;

%% Etalon result
% Polynomial coefficient for up scan
C_polyfit_upscan        = [-1.386212101804218e-20,9.105727297641090e-17,-2.472392336240578e-13,...
    3.633386562263428e-10,-4.005959930756314e-07,-1.852629667917438e-04,0.012603424011227];

% Polynomial coefficient for down scan
C_polyfit_downscan      = [3.101554445027156e-20,-1.983690922113016e-16,5.109167693079562e-13,...
    -6.945628359599489e-10,5.473487595375121e-07,2.649765224741536e-04,-0.885175102651763];

%% Main
[I,y]                   = peakfinder(raw_data,1.2,1.2,1); % find the peaks for dataset
num_scan                = length(I)-1;                    % number of scans 
A                       = zeros(num_scan*2,2);            % Initialize a vector to store Integrated Absorbance 

% down scan clipped region
Region1                 = [90:270]'; 
Region2                 = [980:1180]';
Region3                 = [1660:1730]';

% up scan clipped region
Region4                 = [2000:2200]'; 
Region5                 = [2650:2850]';
Region6                 = [3520:3640]';

% background clipped region
Region7                 = [1840:1890]';

for ii = 1:num_scan
    
    index               = (I(ii):I(ii+1))'; % ii scan region
    midpoint            = round((I(ii)+I(ii+1))/2); 
    
    clipped_data        = raw_data(index);
    offset              = mean(clipped_data(Region7)); % background noise
    clipped_data        = clipped_data - offset;
    
    %% DownScan
    index_downscan      = (Region1(1,1):Region3(end,1))';
    X_downscan          = [Region1;Region2;Region3];
    Y_downscan          = [clipped_data(Region1);clipped_data(Region2);clipped_data(Region3)];
    
    C_absorb            = polyfit(X_downscan,Y_downscan,3);
    I_0                 = polyval(C_absorb,index_downscan); % Fit 
    I_t                 = clipped_data(index_downscan); % Data
    Absorbance_Exp      = -log(I_t./I_0);
    v_downscan          = polyval(C_polyfit_downscan,index_downscan); % wavenumber
    
    % Determine guess values for spectroscopic parameters needed for
    % least-squares fitting
    % Use nlinfit to find best-fit absorbance spectra
    vd                  = vo*(7.1623e-7)*(Tguess/MW)^(1/2);         % Doppler width (cm^-1)
    vo_guesses          = [ -0.6  -0.16];                           % Linecenter frequency of each transition
    peak_guesses        = [ 0.46   0.18];                           % Peak absorbance, used for determing guess of A
    vc_guesses          = [ 0.04   0.04];                           % Collisional width guess
    Aint1_guess         = get_IntArea_guess(vd,vc_guesses(1,1),peak_guesses(1,1));  % Convert peak absorbance guess to integrated absorbance guess
    Aint2_guess         = get_IntArea_guess(vd,vc_guesses(1,2),peak_guesses(1,2));
    Free_Parameters     = [vo_guesses(1,1) vc_guesses(1,1) Aint1_guess vo_guesses(1,2) vc_guesses(1,2)  Aint2_guess]; % Bundle free-parameters

    options             = statset('MaxIter',200);
    estimates           = nlinfit(v_downscan,Absorbance_Exp,@Voigt_Approx_McLean_Vectorized_Fit,Free_Parameters, options);
    Linecenters         = [estimates(1,1) estimates(1,4)];
    Collisional_Widths  = [estimates(1,2) estimates(1,5)];
    Integrated_Areas    = [estimates(1,3) estimates(1,6)];
    A(2*ii-1,1:2) = Integrated_Areas; % store downscan integrated areas
    
    if plot
        
    figure(1)
    plot(v_downscan,Absorbance_Exp,'k');hold on;
    plot(v_downscan,Absorbance_Fit,'r--')
    xlabel('Relative Frequency, cm^{-1}')
    ylabel('Absorbance')
    title('Unsteady Flame Downscan Fitting')
    legend('Data','Best Fit')
    hold off    
%     figure(2)
%     plot(v_downscan,Absorbance_Exp);hold on;
    end
    
    %% UpScan
    index_upscan        = (Region4(1,1):Region6(end,1))';
    X_upscan            = [Region4;Region5;Region6];
    Y_upscan            = [clipped_data(Region4);clipped_data(Region5);clipped_data(Region6)];
    
    C_absorb            = polyfit(X_upscan,Y_upscan,3);
    I_0                 = polyval(C_absorb,index_upscan);
    I_t                 = clipped_data(index_upscan);
    Absorbance_Exp      = -log(I_t./I_0);
    v_upscan            = polyval(C_polyfit_upscan,(index_upscan-1875+1));
    
    % Determine guess values for spectroscopic parameters needed for
    % least-squares fitting
    % Use nlinfit to find best-fit absorbance spectra 
    vd                  = vo*(7.1623e-7)*(Tguess/MW)^(1/2);         % Doppler width (cm^-1)
    vo_guesses          = [ -0.6  -0.184];                          % Linecenter frequency of each transition
    peak_guesses        = [ 0.465  0.16];                           % Peak absorbance, used for determing guess of A
    vc_guesses          = [ 0.04   0.04];                           % Collisional width guess
    Aint1_guess         = get_IntArea_guess(vd,vc_guesses(1,1),peak_guesses(1,1));  % Convert peak absorbance guess to integrated absorbance guess
    Aint2_guess         = get_IntArea_guess(vd,vc_guesses(1,2),peak_guesses(1,2));
    Free_Parameters     = [vo_guesses(1,1) vc_guesses(1,1) Aint1_guess vo_guesses(1,2) vc_guesses(1,2)  Aint2_guess]; % Bundle free-parameters 
    
    options             = statset('MaxIter',200);
    estimates           = nlinfit(v_upscan,Absorbance_Exp,@Voigt_Approx_McLean_Vectorized_Fit,Free_Parameters, options);
    Linecenters         = [estimates(1,1) estimates(1,4)];
    Collisional_Widths  = [estimates(1,2) estimates(1,5)];
    Integrated_Areas    = [estimates(1,3) estimates(1,6)];
    A(2*ii,1:2) = Integrated_Areas;
    
    if plot
    figure(2)
    plot(v_upscan,Absorbance_Exp,'b');hold on;
    plot(v_upscan,Absorbance_Fit,'c--');
    xlabel('Relative Frequency, cm^{-1}')
    ylabel('Absorbance')
    title('Unsteady Flame Upscan Fitting')
    legend('Data','Best Fit')
    hold off
%      plot(v_upscan,Absorbance_Exp)
    end
end

%% Save result
if save 
     save('absorbance_steady.mat','A');
end
%% Post-procesing
if load
     load('absorbance_steady.mat','A');
end

% Solve Temperature
R    =  A(:,1)./A(:,2);
h    =  6.626e-34;          % Prank constant
c    =  3e10;               % Light speed, cm/s
k    =  1.38e-23;           % Boltzmann constant
T0   =  296;                % T0 temperature, K
S1   =  3.53e-20*2.4797e19; % Linestrenth(T0) for P(0,20) transition, cm^-1/atm
S2   =  1.06e-23*2.4797e19; % Linestrenth(T0) for P(1,14) transition, cm^-1/atm
E1   =  806.38;             % Lower state energy for P(0,20), cm^-1
E2   =  2543.05;            % Lower state energy for P(1,14), cm^-1
T    =  h*c/k*(E2-E1)./(log(R)+log(S2/S1)+h*c*(E2-E1)/(k*296)); % Temperature, K
L    =  1.25;                

if plot 
    
figure(3)
plot(T,'r');
xlim([0 1066])
xlabel('Time ms');
ylabel('Temperature K')
title('Temperature Plot')

end 

% Solve CO Concentration
v0_1 =  2059.91;            % P(0,20) transition, cm^-1
Q0   =  112.77498;          % partition function T0
P    =  1;                  % Pressure, atm
S    =  S1*(Q0./Q_CO(T)).*(T0./T).*exp(-h*c*E1/k.*(1./T-1/T0)).*(1-exp(-h*c.*v0_1./(k*T)))*(1-exp(-h*c*v0_1/(k*T0)))^-1; % Linestrength
X_CO =  A(:,1)./(S*P*L);    % Mole Fraction

if plot
figure(4)
plot(X_CO,'b')
xlim([0 1066])
xlabel('Time ms');
ylabel('CO Concentration')
title('CO Concentration Plot')

% j = 533;
% aa = zeros(533,1);
% bb = zeros(533,1);
% cc = zeros(533,1);
% dd = zeros(533,1);
% 
% for j = 1:533
%  aa(j) = T(2*j-1);
%  bb(j) = T(2*j);
%  cc(j) = X_CO(2*j-1);
%  dd(j) = X_CO(2*j);
% end

figure(5);
plot(aa,'r');
hold on;
plot(bb,'b')
xlabel('Time ms');
ylabel('Temperature K')
title('Downscan vs. Upscan Temperature')
legend('downscan','upscan')

figure(6);
plot(cc,'r');
hold on;
plot(dd,'b')
xlabel('Time ms');
ylabel('CO Concentration')
title('Downscan vs. Upscan CO Concentration')
legend('downscan','upscan')
end
