%~~~~~Program to numerically solve SBS optical fields~~~~~~%
%                                                          %
% Author:       Sam Weller                                 %
% Institution:  Royal Holloway University of London        %
% Department:   Electronic Engineering                     %
% Email:        samuel.weller.2017@live.rhul.ac.uk         %
%                                                          %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Variables and Constants
L = 1000; % Waveguide length (m)
Aeff = 85e-12; % Effective area (m^2)

c = 3e8; % Speed of light in a vacuum (m/s)
n = 1.46; % Core index
v = c/n; % Velocity in waveguide

alpha = 0.2; % Wavegiuide attenuation (dB/km)
alpha = alpha/(1000*10*log10(exp(1))); % Recalculate in /m

eta = 0.2; % Polarisation coefficient
g0 = 4.4e-11 / Aeff; % Brillouin gain coefficient scaled by the effective area
g0 = g0 * eta; % Recalculated with polarisation coefficient

a0_input_P = 70e-3; % Pump input power (W)
a1_input_P = 251e-6; % Stokes input power (W)

a0_input = @(t) (a0_input_P) ./(1+exp(-0.0025*t*v)); % Pump input wave (CW)
a1_input = @(t) (a1_input_P) ./(1+exp(-0.0025*t*v)); % Stokes input wave (CW)

% Simulation Parameters
simulation_time = 10*L/v; % Simulation time
N_z = 501; % Number of z increments
zz = linspace(0,L,N_z); % z space
d_z = zz(2) - zz(1); % calc d_z
d_t = d_z / v; % calc d_t
N_t = ceil(simulation_time/(d_t)); % calc num t increments
tt = linspace(0, simulation_time, N_t); % t space
d_t = tt(2)-tt(1); % recalculate d_t

[TT, ZZ] = meshgrid(tt, zz); % Make meshgrid for z and t

% Initialise 
a_0 = a0_input(TT-ZZ/v);
a_1 = a1_input(TT-(L-ZZ)/v);

a_0(1,:) = a0_input(tt);
a_1(N_z,:) = a1_input(tt);

% Solve Loop
tic
for i = 1:N_t-1
      
  a_0(2:N_z, i+1) = a_0(1:N_z-1, i);
  a_1(1:N_z-1, i+1) = a_1(2:N_z, i);
  
  a_0(1, i+1) = a0_input(tt(i+1));
  a_1(end, i+1) = a1_input(tt(i+1));
  
  a0_temp = a_0(:, i+1);
  a1_temp = a_1(:, i+1);

  a_0(:, i+1) = a0_temp.*(1-alpha*v*d_t-g0*v*d_t*a1_temp);
  a_1(:, i+1) = a1_temp.*(1-alpha*v*d_t+g0*v*d_t*a0_temp);
end
toc

P_pump = abs(a_0);
P_stokes = abs(a_1);

%% Plot

steps = 2.*ceil(N_t/100);
steps2 = 2.*ceil(N_z/100);
[ZZ1,TT1] = meshgrid(zz(:,1:steps2:N_z),tt(:,1:steps:N_t));
yLIMS = [0,max(tt.*1e6)];
xLIMS = [0,L];
xticks_v = linspace(0,L,6);
map = parula;

az_angle1 = 30;
el_angle1 = 50;

fig = figure;

subplot(1,2,1)
s = mesh(ZZ1,TT1.*1e6,(P_pump(1:steps2:N_z,1:steps:N_t).')*1e3);
xlim(xLIMS);
ylim(yLIMS);
xticks(xticks_v);
view(az_angle1,el_angle1);
colormap(map);
xlabel('z (m)')
ylabel('t (\mus)')
zlabel('P_1(z,t) (mW)')
title('a)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax.FontWeight = "normal";
grid off;

subplot(1,2,2)
s = mesh(ZZ1,TT1.*1e6,(P_stokes(1:steps2:N_z,1:steps:N_t).')*1e3);
xlim(xLIMS);
ylim(yLIMS);
xticks(xticks_v);
view(az_angle1,el_angle1);
colormap(map);
xlabel('z (m)')
ylabel('t (\mus)')
zlabel('P_2(z,t) (mW)')
title('b)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ax.FontWeight = "normal";
grid off