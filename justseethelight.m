%% Simulating the photoacoustic effect using k-Wave: kwavetest.m
% This example demonstrates simulation of a pressure field generated through 
% the absorption of an externally introduced light pulse.
% The light propagation is simulated using ValoMC and the propagation and
% detection of pressure wavefield is simulated using k-Wave, see 
% http://www.k-wave.org/documentation/k-wave_initial_value_problems.php.
% The example also shows how the computation grid of k-Wave and mesh of
% ValoMC can be made compatible.
% Note that k-Wave uses SI units (e.g. [m]) and ValoMC works in millimetre-scale 
% (e.g. [mm]).

%% k-Wave initialisation
% The initialisation is done as normal in k-Wave.
% Care must be taken at the initialization ValoMC, to make a 
% matching computational simulation area for (see ValoMC initialization)
% the photon propagation simulation.
% if ValoMC doesnot work complile these three line codes.
% mex cpp/2d/MC2Dmex.cpp COMPFLAGS='$COMPFLAGS -O3'
% mex cpp/3d/MC3Dmex.cpp COMPFLAGS='$COMPFLAGS -O3'
% mex cpp/3d/createBH3mex.cpp COMPFLAGS='$COMPFLAGS -O3'
%% 
% photon_count default=1e6
clear all;

% Create the k-Wave grid
PML_size = 20;              % size of the PML in grid points
Nx = 300;           % number of grid points in the x (row) direction
Ny = 300;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

wavelength = 990; 

% Create two internal structures using makeDisk
% discs = makeDisc(Nx, Ny, 55, 55, 5) + makeDisc(Nx, Ny, 75, 85, 10);
% discs = loadImage('/Users/liutieming/Desktop/mice_skull.png');
discs = load("digimouse_brian_300.mat");
discs = discs.V_brain_300;
% discs(1:280,:) = discs(21:300,:);
% discs(281:300,:) = discs(1:20,:);
% Define the acoustic properties 
disc_background = find(discs==0);
disc_scalp = find(discs==1);
disc_skull = find(discs==2);
disc_braintissue = find(discs>2);
% disc_blood = makeDisc(Nx, Ny, 100, 145, 5);
% disc_blood_indices = find(disc_blood==1);
% discs = discs+ disc_blood; % simulate the blood vessel

% medium.sound_speed = 1500*ones(Nx, Ny);    % [m/s]
% medium.sound_speed(disc_indices) = 1800;   % [m/s]
% 
% medium.density = 1000*ones(Nx, Ny);        % [kg/m^3]
% medium.density(:,Ny/2:end) = 1400;
% figure
% imagesc(discs);
% axis square
% set(gca,'xtick',[],'xticklabel',[]);
% set(gca,'ytick',[],'yticklabel',[]);
% export_fig test.png -transparent

%% mouse skull

% convert to density and speed of sound
density = ones([Nx,Ny]); 
sound_speed = ones([Nx,Ny]); 

density(disc_background) = 1000;                           % Minimum density equal to water
density(disc_scalp) = 955;
density(disc_skull) = 1658;
density(disc_braintissue) = 1046;
% density(disc_blood_indices) = 930;

sound_speed(disc_background) = 1500;    % 37C: Speed of sound at body temperature
sound_speed(disc_scalp) = 1450;
sound_speed(disc_skull) = 2850;
sound_speed(disc_braintissue) = 1545;
% sound_speed(disc_blood_indices) = 1555;

% display CT image
% imagesc(CT_resized)
% colormap gray
% colorbar
% title('CT data', FontSize=18)
% xlabel('x [grid points]', FontSize=18)
% ylabel('y [grid points]', FontSize=18)
%% define the properties of the propagation medium

% define the properties of the propagation medium from CT data
medium.sound_speed = sound_speed;   % [m/s]
medium.density     = density;       % [kg/m^-3]
medium.alpha_power = 1.1;
 
% create alpha coefficient mask
alpha_coeff = 2.4e-3  * ones(Nx, Ny);      % attenuation coefficient in water [dB/(MHz^y cm)]
alpha_coeff(disc_scalp) = 5e-2;
alpha_coeff(disc_skull) = 8.83;
alpha_coeff(disc_braintissue) = 6.9e-2;
% alpha_coeff(disc_blood_indices) = 0.05;
% pix = CT_resized .* (8.83/ max(CT_resized(:)));           % scale pixel values between 0 and 8.83/13.3 (attenuation through mice/human skull)
% pix=round(pix);                           % round scales values
% medium.alpha_coeff = alpha_coeff + pix;   % add scaled values to the attenuation of water
medium.alpha_coeff = alpha_coeff;
% display sound speed map
% figure
% imagesc(density)
% colormap winter
% xlabel('x [grid points]', FontSize=18)
% ylabel('y [grid points]', FontSize=18)
% title('Density map', FontSize=18)
% colorbar
% c = colorbar;
% c.Label.String = 'Density [kg/m^3]';
% c.Label.FontSize = 18;

% create the time array
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);
kgrid.t_array = (kgrid.t_array(1):dt:100e-6); % Total simulation time (10ms)
%% Create a ValoMC mesh
% ValoMC uses triangles and tetrahedrons as the basis shape, whereas
% in k-Wave pixels and voxels are used. createGridMesh can be used to
% create a straightforward mapping between the triangles and pixels.
% *Note that K-Wave uses matrices in the format matrix(X,Y), whereas* 
% *MATLAB t(c.f. meshgrid) and ValoMC uses matrix(Y,X)*
% Therefore x and y should be swapped when moving between ValoMC
% arrays and k-Wave arrays

vmcmesh = createGridMesh(kgrid.y_vec*1e3, kgrid.x_vec*1e3); % [m to mm]

%% Define optical coefficients
% For users accustomed to k-Wave, the optical coefficients can be set 
% in similar fashion as in k-Wave, i.e. using multidimensional arrays.
% If one of the arrays defining the medium is given as multidimensional 
% array to ValoMC, the code will assume that the mesh was created
% using 'createGridMesh' and the output fluence will also given as 
% a two dimensional array in solution.grid_fluence. See the example 
% 'Working with pixel and  voxel data' on how to achieve similar 
% assignments using one dimensional indexing.

load("caculate_us_blood.mat");
load("caculate_us_brain.mat");
load("caculate_us_scalp.mat");
load("caculate_us_skull.mat");
load("createFit_blood.mat");
load("createFit_brain.mat");

vmcmedium.scattering_coefficient = 0.01*ones(Nx, Ny);
vmcmedium.absorption_coefficient = 0.001*ones(Nx, Ny);



%% Create a light source 

% Set a light source with a width of 2 mm and cosinic directional profile 
% in -x direction
% the light width is set to be 10mm (LDs or LEDs)
% the angle of the light is set to be 45 deg (pi/4) or 60 deg (pi/6)

% angle = pi/4.8; % 52.5 deg
% angle = pi/4.24; % 47.5 deg
% angle = pi/4; % 45 deg
angle = 0; % 60 deg
% angle = pi/4.5; % 50 deg
% angle = pi/5.14; % 55 deg
% angle = pi/5.54; % 57.5 deg
% angle = pi/6.55; % 62.5 deg


    lightsource1 = findBoundaries(vmcmesh, 'direction', ...
        [0 0], ...
        [-12 -15], ...
        5);

        % lightsource1 = findBoundaries(vmcmesh, 'direction', ...
        % [0 12], ...
        % [0 0], ...
        % 7);

    lightsource2 = findBoundaries(vmcmesh, 'direction', ...
        [0 0], ...
        [12 -15], ...
        5);
    % vmcboundary.lightsource(lightsource1) = {'cosinic'};
    % vmcboundary.lightsource(lightsource2) = {'cosinic'};

    vmcboundary.lightsource(lightsource1) = {'gaussian'};
    vmcboundary.lightsource(lightsource2) = {'gaussian'};

    vmcboundary.lightsource_gaussian_sigma(lightsource1) = 0.2;
    vmcboundary.lightsource_direction(lightsource1,1) = sin(-angle);
    vmcboundary.lightsource_direction(lightsource1,2) = cos(-angle);

    vmcboundary.lightsource_gaussian_sigma(lightsource2) = 0.2;
    vmcboundary.lightsource_direction(lightsource2,1) = sin(angle);
    vmcboundary.lightsource_direction(lightsource2,2) = cos(angle);
    vmcboundary.lightsource_direction_type(lightsource1) = {'relative'};
    vmcboundary.lightsource_direction_type(lightsource2) = {'relative'};
    %% Run the Monte Carlo simulation
    solution = ValoMC(vmcmesh, vmcmedium, vmcboundary);

    %% Compute the initial pressure from the photon fluence
    %
    % Note that since the medium was defined as two dimensional arrays,
    % the output is also given as a two-dimensional array.
    %
    % <html>
    % <font color="red">Corrected explanation</font>
    % </html>

    % Compute the absorbed optical energy density.
    % multiply by
    % 1e6   to convert [J/mm^2] to [J/m^2]
    % 1e-3  to set the total energy in the pulse to 1 mJ
    %
    vmcmedium.absorbed_energy = vmcmedium.absorption_coefficient .* solution.grid_fluence*1e3; % [J/m3]

    % Compute the initial pressure distribution
    % source.p0 = vmcmedium.gruneisen_parameter .* vmcmedium.absorbed_energy;  % [Pa]
    %%

    log_data = vmcmedium.absorbed_energy;
    min_data = min(min(log_data(:)));
    max_data = max(max(log_data(:)));
    norm_data = (log_data-min_data)/(max_data-min_data);

    filename = ['absorbenergy' '_' num2str(angle) '.mat'];
    absorbed_energy = norm_data;
    save(['/Users/liutieming/Documents/MATLAB/ValoMC/result_angle/',filename],'absorbed_energy');

    mean_absorbedenergy = mean(norm_data(disc_braintissue),"all");
    std_absorbedenergy = std2(norm_data(disc_braintissue));
    
    COV(i,1) = std_absorbedenergy/mean_absorbedenergy;
   
    i = i+1;


%% Define the k-Wave sensor mask

% % Define a circular sensor
% 
% % define a binary line sensor
% % sensor.mask = zeros(Nx, Ny);
% % sensor.mask(400, 373:2:627) = 1; % 128 elements
% sensor.mask(20, 22:2:278) = 1; % 128 elements
% 
% % sensor_radius = 4e-3;       % [m]
% % num_sensor_points = 50;     % number of sensor points
% % 
% % sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);
% 
% % set the input arguements: force the PML to be outside the computational
% % grid; switch off p0 smoothing within kspaceFirstOrder2D
% input_args = {'PMLInside', false, 'PMLSize', PML_size, 'Smooth', false, 'PlotPML', false,'RecordMovie',false,'MovieName','PA_skull sim','PlotScale',[-2,2],'PlotFreq',5,'MovieProfile','MPEG-4','MovieArgs',{'FrameRate',40}};
% %% Move the perfectly matched layer (PML) outside of the computation domain and run the acoustic simulation
% % The PML is a layer that absorbs waves for simulating free regions and 
% % is normally contained within the computation  region of k-Wave. 
% % For a more straightforward mapping between the
% % computation region of k-Wave and ValoMC, the PML is moved outside
% % of the computation region.
% 
% sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLInside', false);
% 
% %% 
% % assign the time reversal data
% sensor.time_reversal_boundary_data = sensor_data;
% 
% % run the time reversal reconstruction
% p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
% 
% % add first order compensation for only recording over a half plane
% p0_recon = 2 * p0_recon;
% 
% % repeat the FFT reconstruction for comparison
% % p_xy = kspaceLineRecon(sensor_data.', dy, kgrid.dt, medium.sound_speed, ...
% %     'PosCond', true, 'Interp', '*linear');
% p_xy = kspaceLineRecon(sensor_data.', dy, kgrid.dt, 1500, ...
%     'PosCond', true, 'Interp', '*linear');
% % define a second k-space grid using the dimensions of p_xy
% [Nx_recon, Ny_recon] = size(p_xy);
% kgrid_recon = kWaveGrid(Nx_recon, kgrid.dt * 1500, Ny_recon, dy);
% 
% % resample p_xy to be the same size as source.p0
% p_xy_rs = interp2(kgrid_recon.y, kgrid_recon.x - min(kgrid_recon.x(:)), p_xy, kgrid.y, kgrid.x - min(kgrid.x(:)));

%% Plot the solution Absorbed energy [norm.]

% figure('rend','painters','pos',[10 10 1200 400])

% plot the initial pressure and sensor distribution

figure;
log_data = log10(vmcmedium.absorbed_energy);
min_data = min(min(log_data(:)));
max_data = max(max(log_data(:)));
norm_data = (log_data-min_data)/(max_data-min_data);



% We have to swap x and y again 
% imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, vmcmedium.absorbed_energy);
imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, norm_data);

xticks = ([-15 -10 -5 0 5 10 15]);
yticks = ([-15 -10 -5 0 5 10 15]);
xticklabels = ({'-15','-10','-5','0','5','10','15'});
yticklabels = ({'-15','-10','-5','0','5','10','15'});
xlabel('x-position [mm]');
ylabel('y-position [mm]');
set(gca,'xtick',xticks,'xticklabel',xticklabels)
colormap(getColorMap);

c = colorbar;  % create a colorbar
axis image;
title('Absorbed energy [norm.]');
set(gca,'FontSize',18);
%% plot the absorbed_energy [J/m3]

figure;
log_data = log10(vmcmedium.absorbed_energy);
min_data = min(min(log_data(:)));
max_data = max(max(log_data(:)));
norm_data = (log_data-min_data)/(max_data-min_data);



% We have to swap x and y again 
% imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, vmcmedium.absorbed_energy);
imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, vmcmedium.absorbed_energy);
xlabel('x-position [mm]');
ylabel('y-position [mm]');
colormap(getColorMap);

c = colorbar;  % create a colorbar
% colorbar('ylim',[0,4])
axis image;
title('Absorbed energy [J/m3]');
set(gca,'FontSize',18);
% save AE_700_53.mat , vmcmedium.absorbed_energy;
%% Initial pressure [Pa]
    
figure;

log_data_source = log10(source.p0);
min_data_source = min(min(log_data_source(:)));
max_data_source = max(max(log_data_source(:)));
norm_data_source = (log_data_source-min_data_source)/(max_data_source-min_data_source);

% imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, transpose(norm_data_source));
imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, source.p0);
% imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, transpose(source.p0), [min(source.p0(:)) ...
%         max(source.p0(:))]);

% imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, transpose(norm_data_source));
   
colormap(getColorMap);
xlabel('x-position [mm]');
ylabel('y-position [mm]');
c = colorbar;  % create a colorbar
axis image;
title('Initial pressure [Pa]');
set(gca,'FontSize',18);

%% plot the simulated sensor data
% figure;
% imagesc(sensor_data,  [min(sensor_data(:)) max(sensor_data(:))]);
% colormap(getColorMap);
% ylabel('Sensor Position');
% xlabel('Time Step');
% c = colorbar;  % create a colorbar
% colorbar;
% title('Sensor data');
% set(gca,'FontSize',18);
%% absorbed_energy varies with depth[J/m3]
    
% figure;
% load AE_depth_1064.mat
% AE_depth_1064 = AE_depth;
% load AE_depth_700.mat
% AE_depth_700 = AE_depth;

% load AE_1064_53.mat
% AE_depth_1064 = AE_1064_53;
% load AE_700_53.mat
% AE_depth_700 = AE_700_53;
% 
% plot(kgrid.y_vec*1e3,AE_depth_1064,LineWidth=2,Color='b');
% hold on;
% plot(kgrid.y_vec*1e3,AE_depth_700,LineWidth=2,Color='k',LineStyle=':');
% grid on
% xlim([-15 15])
% xlabel('y-position [mm]');
% ylabel('Absorbed energy [J/mm3]');
% legend('1064nm','700nm')                    
% title('Absorbed energy varies with depth');
% set(gca,'FontSize',18);

%% 
% plot the reconstructed initial pressure with positivity condition
% figure;
% imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0_recon+sensor.mask);
% c = colorbar;  % create a colorbar
% colorbar;
% ylabel('x-position [mm]');
% xlabel('y-position [mm]');
% axis image;
% colorbar;
%% Normalized absorbed energy
figure;
log_data = vmcmedium.absorbed_energy;
min_data = min(min(log_data(:)));
max_data = max(max(log_data(:)));
norm_data = (log_data-min_data)/(max_data-min_data);


light_mask = zeros(Nx,Ny);
lightbegin1_column = min(lightsource1)-(fix(min(lightsource1)./Nx))*Nx;
lightbegin1_row = fix(min(lightsource1)./Nx);
lightend1_column = max(lightsource1)-(fix(max(lightsource1)./Nx))*Nx;
lightend1_row = fix(max(lightsource1)./Nx);
lightbegin2_column = min(lightsource2)-(fix(min(lightsource2)./Nx))*Nx;
lightbegin2_row = fix(min(lightsource1)./Nx);
lightend2_column = max(lightsource2)-(fix(max(lightsource2)./Nx))*Nx;
lightend2_row = fix(max(lightsource2)./Nx);
light_mask(lightbegin1_row:lightend1_row,lightbegin1_column:lightend1_column) = 0.5;
light_mask(lightbegin2_row:lightend2_row,lightbegin2_column:lightend2_column) = 0.5;
% We have to swap x and y again 
% imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, vmcmedium.absorbed_energy);
imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, norm_data+light_mask);
xlabel('x-position [mm]');
ylabel('y-position [mm]');
colormap(getColorMap);
c = colorbar;  % create a colorbar
axis image;
title('Normalized absorbed energy');
set(gca,'FontSize',18);
%% 
% us 计散射系数计算
% lam = 700 nm, 810 nm, 850 nm, 980 nm and 1064 nm  
syms a b;
eq1 = a*(700/500)^(-b)==85.7;
eq2 = a*(1064/500)^(-b)==59.6;
[a,b] = solve([eq1,eq2],[a,b])
lam = 980;
us = a*(lam/500)^(-b)

%% calculate the Coefficient of Variation COV (inside the brain)

mean_absorbedenergy = mean(reshape(norm_data(disc_braintissue)))





