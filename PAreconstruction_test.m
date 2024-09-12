% Using An Ultrasound Transducer As A Sensor Example
%
% This example shows how an ultrasound transducer can be used as a detector
% by substituting a transducer object for the normal sensor input
% structure. It builds on the Defining An Ultrasound Transducer and
% Simulating Ultrasound Beam Patterns examples.
%
% author: Bradley Treeby
% date: 1st August 2011
% last update: 4th June 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011-2017 Bradley Treeby

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

clearvars;

% simulation settings
DATA_CAST = 'single';

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 20;            % [grid points]
PML_Y_SIZE = 10;            % [grid points]
PML_Z_SIZE = 10;            % [grid points]

% set total number of grid points not including the PML
Nx = 128 - 2*PML_X_SIZE;    % [grid points]
Ny = 128 - 2*PML_Y_SIZE;    % [grid points]
Nz = 64 - 2*PML_Z_SIZE;     % [grid points]

% set desired grid size in the x-direction not including the PML
x = 40e-3;                  % [m]

% calculate the spacing between the grid points
dx = x/Nx;                  % [m]
dy = dx;                    % [m]
dz = dx;                    % [m]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
medium.sound_speed = 1500;      % [m/s]
medium.density = 1000;          % [kg/m^3]
medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
medium.BonA = 6;

% create the time array
t_end = 40e-6;                  % [s]
kgrid.makeTime(medium.sound_speed, [], t_end);

% =========================================================================
% DEFINE THE SOURCE
% =========================================================================

% create source mask
source.p_mask = makeBall(Nx, Ny, Nz, round(25e-3/dx), Ny/2, Nz/2, 3)...
    + makeBall(Nx, Ny, Nz, round(8e-3/dx), Ny/4, Nz/2, 3);

% define properties of the input signal
source_strength = 10;          % [Pa]
tone_burst_freq = 0.5e6;        % [Hz]
tone_burst_cycles = 5;

% create the input signal using toneBurst 
source.p = source_strength .* toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
% p0 = source_strength * source.p_mask;
% source.p = smooth(p0, true);



% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% physical properties of the transducer
transducer.number_elements = 72;    % total number of transducer elements
transducer.element_width = 1;       % width of each element [grid points/voxels]
transducer.element_length = 12;     % length of each element [grid points/voxels]
transducer.element_spacing = 0;     % spacing (kerf  width) between the elements [grid points/voxels]
transducer.radius = inf;            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = 1540;                  % sound speed [m/s]
transducer.focus_distance = 25e-3;              % focus distance [m]
transducer.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
transducer.steering_angle = 0;                  % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Rectangular';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
transducer.active_elements = zeros(transducer.number_elements,1,1);
transducer.active_elements(21:52,:,:) = 1;

% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties;

voxelPlot(double(source.p_mask | transducer.all_elements_mask) , 'AxisTight', false, 'Color', [1 0 0], 'Transparency', 0.5);

% transducer.plot

% =========================================================================
% RUN THE SIMULATION
% =========================================================================

% set the input settings
input_args = {'DisplayMask', transducer.active_elements_mask, ...
    'PMLInside', false, 'PlotPML', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
    'DataCast', DATA_CAST, 'PlotScale', [-1/4, 1/4] * source_strength};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, transducer, input_args{:});

% extract a single scan line from the sensor data using the current
% beamforming settings
scan_line = transducer.scan_line(sensor_data);
[~, as_1] = spect(sensor_data(1, :), 1/kgrid.dt);
[~, as_2] = spect(sensor_data(2, :), 1/kgrid.dt);
[f, as_3] = spect(sensor_data(3, :), 1/kgrid.dt);


%% =========================================================================
% VISUALISATION
% =========================================================================

% plot the recorded time series
figure;
stackedPlot(kgrid.t_array * 1e6, sensor_data(1,:));
xlabel('Time [\mus]');
ylabel('Transducer Element');
title('Recorded Pressure');
scaleFig(1, 1.5);

% plot the scan line
figure;
plot(kgrid.t_array * 1e6, scan_line * 1e-6, 'k-');
xlabel('Time [\mus]');
ylabel('Pressure [MPa]');
title('Scan Line After Beamforming');
%% 
% plot the recorded time series
figure;
stackedPlot(kgrid.t_array * 1e6, {'Sensor Position 1', 'Sensor Position 2', 'Sensor Position 3'}, sensor_data);
xlabel('Time [\mus]');

% plot the corresponding amplitude spectrums
figure;
plot(f .* 1e-6, as_1 ./ max(as_1(:)), 'k-', ...
     f .* 1e-6, as_2 ./ max(as_1(:)), 'b-', ...
     f .* 1e-6, as_3 ./ max(as_1(:)), 'r-');
legend('Sensor Position 1', 'Sensor Position 2', 'Sensor Position 3');
xlabel('Frequency [MHz]');
ylabel('Normalised Amplitude Spectrum [au]');
f_max = medium.sound_speed / (2 * dx);
set(gca, 'XLim', [0, f_max .* 1e-6]);

%% PA reconstruction

sensor.mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
sensor.mask(1, :, :) = 1;

% create the time array
% kgrid.makeTime(medium.sound_speed);

% set the input arguements
input_args = {'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], 'PMLInside', false, ...
    'PlotPML', false, 'Smooth', false, 'DataCast', 'single'};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% reset the initial pressure
source.p0 = 0;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

% run the time-reversal reconstruction
p0_recon = kspaceFirstOrder3D(kgrid, medium, source, transducer, input_args{:});

% add first order compensation for only recording over a half plane
p0_recon = 2 * p0_recon;

% apply a positivity condition
p0_recon(p0_recon < 0) = 0;
%% 
% plot the reconstructed initial pressure
plot_scale = [-10, 10];
figure;
subplot(2, 2, 1);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, squeeze(p0_recon(:, :, 1)), plot_scale);
title('x-y plane');
axis image;

subplot(2, 2, 2);
imagesc(kgrid.z_vec * 1e3, kgrid.x_vec * 1e3, squeeze(p0_recon(:, Ny/2, :)), plot_scale);
title('x-z plane');
axis image;
xlabel('(All axes in mm)');

subplot(2, 2, 3);
imagesc(kgrid.z_vec * 1e3, kgrid.y_vec * 1e3, squeeze(p0_recon(32, :, :)), plot_scale);
title('y-z plane');
axis image;
colormap(getColorMap);

