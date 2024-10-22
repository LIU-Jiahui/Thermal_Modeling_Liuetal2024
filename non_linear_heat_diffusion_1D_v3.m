clear; close all;
% 1D code to solve for heat diffusion equation taking into account
% radiogenic heat production and non-linear conductivity
% Here we use an implicit solver in order to easily absorbe the
% non-linearity and to prescibe large timesteps

%% Load interpolated phase diagrams
load('phase_diagrams/volcanics_1687MTK2_mb.mat','volcanics_1687MTK2_mb')
load('phase_diagrams/volcanics_1687MTK2_dryer_mb.mat','volcanics_1687MTK2_dryer_mb')
load('phase_diagrams/volcanics_12WT66_1_mb.mat','volcanics_12WT66_1_mb')
load('phase_diagrams/volcanics_12WT66_1_dryer_mb.mat','volcanics_12WT66_1_dryer_mb')
load('phase_diagrams/TTG_AveTTG2_mp.mat','TTG_AveTTG2_mp')
load('phase_diagrams/TTG_AveTTG2_0_5H2O_mp.mat','TTG_AveTTG2_0_5H2O_mp')
load('phase_diagrams/TTG_AveTTG1_mp.mat','TTG_AveTTG1_mp')
load('phase_diagrams/pelite_L17_mp.mat','pelite_L17_mp')
load('phase_diagrams/granite_2715_mp.mat','granite_2715_mp')
load('phase_diagrams/KLB1_mantle.mat','KLB1_mantle')



%% Define parameters
h       = 280e3;    % model height
hc      = 60e3;     % crust thickness
dk      = 30e3;     % thickness of high radiogenic heat production

tol     = 1e-6;     % non linear relative tolerance
max_ite = 64;       % maximum number of non-linear iteration
nz      = 281;      % number of grid points
max_time= 60e6;     % years
timestep= 1e6;      % years
Tsurface= 25;

s_in_yr = 365*24*3600;
dt      = s_in_yr*timestep;  
max_ts  = max_time/timestep;     %
dz      = h/(nz-1);
z       = 0:dz:h;% define depth array

% material properties
mm_peridotite   = 140.7;    % molar mass g/mol
mm_crust        = 250;      % molar mass g/mol

rho_c   = 2800;     % crust density
rho_m   = 3300;     % mantle density

A_uc    = 1.2e-6;   % radiogenic heating uW/m-3 - upper crust
A_mc    = 2.3e-6;     % radiogenic heating uW/m-3 - middle crust
A_lc    = 1.1e-6;   % radiogenic heating uW/m-3 - lower crust 
A_M     = 0.006e-6;      % radiogenic heating uW/m-3 - upper mantle

K       = 273;
Ttop    = K;
Tbot    = 1550+K;     

% conductivity param
b1       = 5.3;
c1       = 0.0015;
d0      = 1.753e-2;
d1      = -1.0365e-4;
d2      = 2.2451e-7;
d3      = -3.4071e-11;


%% memory allocation
T0      = zeros(nz,1); T0(1) = Ttop; T0(nz) = Tbot;

dT      = (Tbot-Ttop)/(nz-1);
T0(:)   = Ttop:dT:Tbot;
T1      = zeros(nz,1);      % °C
Tini    = zeros(nz,1);      % °C
time    = zeros(max_ts,1); 
S       = zeros(nz,1);      % Source term

c       = zeros(nz-1,1);    % conductivity 
rho     = zeros(nz-1,1);    % density km/m3
cp      = zeros(nz-1,1);    % heat capacity 
kappa   = zeros(nz-1,1);    % thermal diffusivity

A       = zeros(nz);
b       = zeros(nz,1);
alpha   = dt/(dz*dz);

% define density
for i = 1:nz-1
    if (z(i)<=hc)
        rho(i) = 2800;
    else
        rho(i) = 3300;
    end
end
P_kbar = cumsum(rho.*9.81.*dz)./1e8;
% define radiogenic heat production (constant so outside calculation loop)
for i=1:nz
    if (z(i)>= 0 && z(i)<= 2/3*dk)           % 0 km <= TTG or others <=20 km
        S(i)=A_uc;    
    elseif (z(i)> 2/3*dk && z(i)<= 5/3*dk)  % 20 km <= High Ar <= 50 km
        S(i)=A_mc;    
    elseif (z(i)> 5/3*dk && z(i)<= hc)      % 50 km <= TTG <=60 km
        S(i)=A_lc;    
    elseif (z(i)> hc)                      % 60 km <= mantle <=280 km
        S(i)=A_M;   
    end
end

%load initial temperature
data1=load("Dataset 1.csv"); % Extracted from Reimink and Smye, Nature, 2024
t_reimink   = data1(:,1);
d_reimink   = data1(:,2);

% Interpolation
Tstart      = interp1(d_reimink,t_reimink,z./1000,"linear",'extrap') + K + Tsurface;
T0          = Tstart;

Tmat        = zeros(nz,max_ts+1);
Tmat(:,1)   = Tstart;

%% non-linear implicit solver
for ts=1:max_ts
   time(ts) = ts * dt / (s_in_yr*1e6);
    
   conv = 1.0;
   ite  = 1;
   while conv > tol && ite < max_ite

        for i=1:nz-1
            if (z(i)>= 0 && z(i)<= 2/3*dk)           % 0 km <= TTG or others <=20 km
                cp(i) = 0.20*(pelite_L17_mp(P_kbar(i),T0(i)-K))+0.20*(volcanics_1687MTK2_dryer_mb(P_kbar(i),T0(i)-K))+0.60*(TTG_AveTTG1_mp(P_kbar(i),T0(i)-K));
            elseif (z(i)> 2/3*dk && z(i)<= 5/3*dk)  % 20 km <= High Ar <= 50 km
                cp(i) = 0.40*(pelite_L17_mp(P_kbar(i),T0(i)-K))+0.3*(volcanics_12WT66_1_dryer_mb(P_kbar(i),T0(i)-K))+0.3*(granite_2715_mp(P_kbar(i),T0(i)-K));
            elseif (z(i)> 5/3*dk && z(i)<= hc)      % 50 km <= TTG <=60 km
                cp(i) = TTG_AveTTG2_mp(P_kbar(i),T0(i)-K);
            elseif (z(i)> hc && z(i)<= h)                      % 60 km <= mantle <=280 km
                cp(i) = KLB1_mantle(P_kbar(i),T0(i)-K);
            end
        end

        for i=1:nz-1
            if (z(i) <= dk)
                kappa(i)= (0.603+1.902*exp(-1*T0(i)/286.18))*1e-6;
                c(i)= rho(i)*cp(i)*kappa(i);
            elseif (z(i)>dk && z(i)<=hc)
                kappa(i) = (0.558+7.862*exp(-1*T0(i)/176.54))*1e-6;
                c(i)= rho(i)*cp(i)*kappa(i);
            else
                %                 c(i)=4.08*(298/T0(i)).^0.406;       % McKenzie 2005, equation 6
                c(i)=b1/(1+c1*(T0(i)-273))+d0*(T0(i)-273)^0+d1*(T0(i)-273)^1+d2*(T0(i)-273)^2+d3*(T0(i)-273)^3; % McKenzie 2005, equation 4
            end
        end


        for i=1:nz
            if i == 1 || i == nz
                b(i) = T0(i);
                A(i,i) = 1.0;
            else
                b(i) = T0(i)+ dt*S(i)/(rho(i)*cp(i));
                
                A(i,i)    = 1.0 + alpha*(c(i)/(rho(i)*cp(i))) + alpha*(c(i-1)/(rho(i-1)*cp(i-1)));
                A(i,i+1)  = - alpha*(c(i)/(rho(i)*cp(i)));
                A(i,i-1)  = - alpha*(c(i-1)/(rho(i-1)*cp(i-1)));
            end
        end
        Tini    = T1;
        T1      = A\b;
        conv    = norm(T1-Tini);
        ite     = ite + 1;
        if ite == max_ite
           disp('maximum number of non-linear iteration has been reached...')
        end
    end

    T0      = T1;
    Tmat(:,ts+1) = T1;
    figure(1)
    plot(T1-K,-z./1000)
    xlim([Ttop-K,Tbot-K])
    ylim([-80 0])
    xlabel('Temperature [°C]');
    ylabel('Depth [km]');
    cap   = sprintf("Temperature evolution - %.2f Myrs",time(ts));
    title(cap);
end


%% Metamorphic geothermos
Filename2 = 'WTFP_PT_SameAge.xlsx';
Filecontent2 = readtable(Filename2);
data2 = Filecontent2{:,["T_C","P_GPa"]};
Ts2 = data2(:,1);
Ps2 = data2(:,2);
Ts2 = Ts2 + 273;
Ps2 = Ps2.*33;

Filename3 = 'WTFP_PT_all.xlsx';
Filecontent3 = readtable(Filename3);
data3 = Filecontent3{:,["T_C","P_GPa"]};
Ts3 = data3(:,1);
Ps3 = data3(:,2);
Ts3 = Ts3 + 273;
Ps3 = Ps3.*33;



%% Visualize results

figure(2);
% Plotting implicit solution
plot(T1-K,z./1000,'r',Tstart-K,z./1000,'k-.');
axis ij;

title(['Time = ',num2str(max_time/1e6),' Myr']);
xlabel('Temperature (C)');
ylabel('Depth (km)');
hold on
errorbar(Ts2-273,Ps2,3.3,3.3,50,50,'bo','MarkerSize',9,'MarkerFaceColor','b')
plot(Ts2-273,Ps2,'bo','MarkerSize',9,'MarkerFaceColor','b');
hold on
plot(Ts3-273,Ps3,'ko');
hold on
xlim([0 1000])
ylim([0 80])
set(gca,'YDir','reverse');

figure(3)
for i = 1:5:length(time)
    if i < length(time)
        plot(Tmat(:,i)-K,z/1000,'r')
        axis ij;
        hold on
    end
end
ylim([0 80])
xlim([0 1000])

