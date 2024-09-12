%% Creating an inhomogeneous medium: inhomogeneous.m
% This example demonstrates how define inhomogeneous optical properties
clear all;

% Create the k-Wave grid
PML_size = 20;              % size of the PML in grid points
Nx = 712;           % number of grid points in the x (row) direction
Ny = 512;           % number of grid points in the y (column) direction
dx = 0.01e-3;        % grid point spacing in the x direction [m]
dy = 0.01e-3;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);


vmcmesh = createGridMesh(kgrid.y_vec*1e3, kgrid.x_vec*1e3); 
%% calculate absorption coefficient (ua mm-1) and scattering coefficient (us mm-1)
% paper: Toward whole-brain in vivo optoacoustic angiography of rodents: modeling and experimental observations
% absorption coefficient of whole blood: us_blood
% absorption coefficient of brain tissue: us_brain
ua_blood_700 = 0.18;
ua_blood_750 = 0.27;
ua_blood_800 = 0.39;
ua_blood_850 = 0.5;
ua_blood_900 = 0.55;
ua_blood_950 = 0.64;
ua_blood_1000 = 0.64;
ua_blood_1050 = 0.48;
ua_blood_1100 = 0.34;
ua = [700:50:1100];
ua_blood = [0.18 0.27 0.39 0.5 0.55 0.64 0.64 0.48 0.34];

ua_brain_700 = 0.16;
ua_brain_750 = 0.15;
ua_brain_800 = 0.145;
ua_brain_850 = 0.14;
ua_brain_900 = 0.138;
ua_brain_950 = 0.165;
ua_brain_1000 = 0.195;
ua_brain_1050 = 0.155;
ua_brain_1100 = 0.12;
ua_brain = [0.16 0.15 0.145 0.14 0.138 0.165 0.195 0.155 0.12];
ua_brain = ua_brain./4;
% ua fitting functions
fitresult_blood = fit(ua', ua_blood', 'poly7');
createFit_blood = @(x) fitresult_blood.p1*x.^7 + fitresult_blood.p2*x.^6 + fitresult_blood.p3*x.^5 + fitresult_blood.p4*x.^4 + fitresult_blood.p5*x.^3 + fitresult_blood.p6*x.^2 + fitresult_blood.p7*x + fitresult_blood.p8;

fitresult_brain = fit(ua', ua_brain', 'poly8');
createFit_brain = @(x) fitresult_brain.p1*x.^8 + fitresult_brain.p2*x.^7 + fitresult_brain.p3*x.^6 + fitresult_brain.p4*x.^5 + fitresult_brain.p5*x.^4 + fitresult_brain.p6*x.^3 + fitresult_brain.p7*x.^2 + fitresult_brain.p8*x + fitresult_brain.p9;
% ua of scalp and skull = 0.0186 [mm-1]

figure(1);
ua = [700:10:1100];
absorption_blood = createFit_blood(ua);
absorption_brain = createFit_brain(ua);
grid on
plot([700 1100],[0.0186 0.0186],LineWidth=4,Color='k'); % scalp
hold on
plot([700 1100],[0.0136 0.0136],LineWidth=4,Color='g'); % skull
hold on
plot(ua,absorption_brain,LineWidth=4,Color='b');
hold on
plot(ua,absorption_blood,LineWidth=4,Color='r');


ylim([0 0.8])
xlabel('Wavelength [nm]');
ylabel('Absorption coefficients \mu_a [mm^{-1}]');
legend('\mu_a of scalp','\mu_a of skull','\mu_a of brain tissue','\mu_a of blood','Location','northwest') ;
set(gca,'FontSize',18);
%% calculate absorption coefficient (ua mm-1) and scattering coefficient (us mm-1)

% us fitting functions
caculate_us_brain = inline('1245*x.^-1.05 + 2*10.^10*x.^-4');
caculate_us_scalp = inline('16.3*x.^-0.28 + 3.4*10.^10*x.^-4');
caculate_us_skull = inline('22.7*x.^-0.23');
caculate_us_blood = inline('0.01*x + 9');
us = [700:10:1100];
us_brain = caculate_us_brain(us);
us_scalp = caculate_us_scalp(us);
us_skull = caculate_us_skull(us);
us_blood = caculate_us_blood(us) - 10*createFit_blood(us);
figure(2);
plot(us,us_scalp,LineWidth=4,Color='k');
hold on
plot(us,us_skull,LineWidth=4,Color='g');
grid on
plot(us,us_brain,LineWidth=4,Color='b');
hold on
plot(us,us_blood,LineWidth=4,Color='r');
grid on
ylim([0 25])
xlabel('Wavelength [nm]');
ylabel('Scattering coefficients \mu_s [mm^{-1}]');
legend('\mu_s of scalp','\mu_s of skull','\mu_s of brain tissue','\mu_s of blood','Location','northwest') ;
set(gca,'FontSize',18);

%% Give varying optical parameters
% paper:Optical windows for head tissues in near-infrared and short-wave infrared regions: Approaching transcranial light applications 


% for depth = 0.5:0.1:5
%     for diameter = 0.2:0.1:0.5
%         for wavelength = 700:10:1100

% depth = 2;
% diameter = 0.3;
% wavelength = 940;
load('/Users/liutieming/Documents/MATLAB/ValoMC/minivessel2D.mat');

discs = ones(712,512);
discs(1:50,:) = 4; % background
discs(51:70,:) = 2; % scalp
discs(71:100,:) = 3; % skull
discs(101:150,:) = 0; % brain
discs(151:662,:) = minivessel(:,:); % blood vessel % depth in the brain = 0.5-5 mm
discs(663:end,:) = 0; % brain
% Define the acoustic properties 
disc_background = find(discs==4);
disc_scalp = find(discs==2);
disc_skull = find(discs==3);
disc_braintissue = find(discs==0);
disc_blood = find(discs==1);


% for wavelength = 700:10:1100

% Set constant background coefficients
vmcmedium.scattering_coefficient = 0.01*ones(Nx, Ny);
vmcmedium.absorption_coefficient = 0.001*ones(Nx, Ny);

vmcmedium.absorption_coefficient(disc_scalp) = 0.0186;
vmcmedium.absorption_coefficient(disc_skull) = 0.0186; 
vmcmedium.absorption_coefficient(disc_braintissue) = createFit_brain(wavelength);
vmcmedium.absorption_coefficient(disc_blood) = createFit_blood(wavelength);

vmcmedium.scattering_anisotropy(disc_scalp) = caculate_us_scalp(wavelength);
vmcmedium.scattering_anisotropy(disc_skull) = caculate_us_skull(wavelength);
vmcmedium.scattering_anisotropy(disc_braintissue) = caculate_us_brain(wavelength);
vmcmedium.scattering_anisotropy(disc_blood) = caculate_us_blood(wavelength) - 10*createFit_blood(wavelength);

vmcmedium.refractive_index = 1.4*ones(Nx, Ny);

% Define the Gruneisen parameter describing photoacoustic efficiency
vmcmedium.gruneisen_parameter = 0.02*ones(Nx, Ny);      % [unitless]
vmcmedium.gruneisen_parameter(disc_skull) = 0.9;
vmcmedium.gruneisen_parameter(disc_scalp) = 0.9;
vmcmedium.gruneisen_parameter(disc_braintissue) = 0.9;
vmcmedium.gruneisen_parameter(disc_blood) = 0.9;

% Resize the fields in vmcmedium so that they match the number of elements in the mesh
% vmcmedium = createMedium(vmcmesh,vmcmedium);


%% Plot the parameter distribution
% The resulting absorption coefficient distribution is shown in the figure
% below. Shown is also the circle that was used to select the elements.
%
% Note that some triangles intersect the circle, resulting in a poor
% representation of the circular shape. The representation can be improved
% by a better mesh or by incresing the discretisation size.
%
% <<circle.png>>

%% Create a light source
photon_count = 1e7;
line_start = [0 -20];
line_end = [0 0];
line_width = 5; % a light source with a width of 5 mm
lightsource = findBoundaries(vmcmesh,'direction',line_start,line_end,line_width);
vmcboundary.lightsource(lightsource) = {'gaussian'};
vmcboundary.lightsource_gaussian_sigma(lightsource) = 0.2;
% vmcboundary.lightsource_direction_type(lightsource) = {'relative'};
%% Run the Monte Carlo simulation

solution = ValoMC(vmcmesh, vmcmedium, vmcboundary);

%% 
vmcmedium.absorbed_energy = vmcmedium.absorption_coefficient .* solution.grid_fluence*1e3; % [J/m3]

% Compute the initial pressure distribution
source.p0 = vmcmedium.gruneisen_parameter .* vmcmedium.absorbed_energy;  % [Pa]
%%  save data absorbed energy
% filename = ['absorbenergy' '_' num2str(wavelength) '_' num2str(depth) '_' num2str(diameter) '.mat'];
% absorbed_energy = vmcmedium.absorbed_energy;
% save(['/Users/liutieming/Documents/MATLAB/ValoMC/result/modified_absorbedenergy/',filename],'absorbed_energy');

%%  save data absorbed energy
% filename = ['absorbenergy' '_' num2str(wavelength) .mat'];
% absorbed_energy = vmcmedium.absorbed_energy;
% save(['/Users/liutieming/Documents/MATLAB/ValoMC/result/absorbedenergy_wavelength/',filename],'absorbed_energy');
% end

%%  save data Initial pressure [Pa]
% filename = ['initialpressure' '_' num2str(wavelength) '_' num2str(depth) '_' num2str(diameter)  '.mat'];
% initialpressure = source.p0;
% save(['/Users/liutieming/Documents/MATLAB/ValoMC/result/initialpressure/',filename],'initialpressure');
%         end
%     end
% end

%% Plot all the result of wavelengths(700 to 1100nm)
clear clc

for wavelength = 700:10:1100
    filename = ['/Users/liutieming/Documents/MATLAB/ValoMC/result/absorbedenergy_wavelength/absorbenergy' '_' num2str(wavelength)  '.mat'];
    load(filename);
    eval(['absorbenergy','_',num2str(wavelength),'=','absorbed_energy',';']);
end

%% 

PML_size = 20;              % size of the PML in grid points
Nx = 712;           % number of grid points in the x (row) direction
Ny = 512;           % number of grid points in the y (column) direction
dx = 0.01e-3;        % grid point spacing in the x direction [m]
dy = 0.01e-3;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);
% tiledlayout(5,8);

for i = 1:8:40
    % nexttile
    figure
    wavelength=(i-1)*10+700;
    eval(['maxenergy(i)=max(max(absorbenergy_',num2str(wavelength),'));']);
    eval(['minenergy(i)=min(min(absorbenergy_',num2str(wavelength),'));']);
    maxmaxenergy = max(maxenergy);
    minminenergy = min(minenergy);
    eval(['norm_',num2str(wavelength),'=','(absorbenergy_',num2str(wavelength),'-','minminenergy)','/','(maxmaxenergy-minminenergy)',';']);
    eval(['imagesc(kgrid.y_vec*1e3, (kgrid.x_vec*1e3+3.56), norm_',num2str(wavelength),')',';']);
    colormap(getColorMap);
    caxis([0,1]);
    hold on
    plot([-2.5 2.5],[0.5 0.5],'--k',[-2.5 2.5],[0.7 0.7],'--k',[-2.5 2.5],[1 1],'--k');
    xticks = ([-4 -3 -2 -1 0 1 2 3 4]);
    yticks = ([-4 -3 -2 -1 0 1 2 3 4]);
    xticklabels = ({'-4','-3','-2','-1','0','1','2','3','4'});
    yticklabels = ({'-3','-2','-1','0','1','2','3'});
    % xlabel('x-position [mm]');
    ylabel('y-position [mm]');
    % set(gca,'xtick',xticks,'xticklabel',xticklabels)
    set(gca,'xtick',[],'xticklabel',[])
    % set(gca,'ytick',[],'yticklabel',[])
    c = colorbar;  % create a colorbar
    axis image;
    set(gca,'FontSize',30);
end

%colormap(getColorMap);

%% Plot Absorbed energy [norm.]

 figure;
% log_data = log10(vmcmedium.absorbed_energy);
log_data = vmcmedium.absorbed_energy;
min_data = min(min(log_data(:)));
max_data = max(max(log_data(:)));
norm_data = (log_data-min_data)/(max_data-min_data);



% We have to swap x and y again 
% imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, vmcmedium.absorbed_energy);
imagesc(kgrid.y_vec*1e3, (kgrid.x_vec*1e3+3.56), norm_data);
hold on
plot([-2.5 2.5],[0.5 0.5],'--k',[-2.5 2.5],[0.7 0.7],'--k',[-2.5 2.5],[1 1],'--k');
xticks = ([-4 -3 -2 -1 0 1 2 3 4]);
yticks = ([-4 -3 -2 -1 0 1 2 3 4]);
xticklabels = ({'-4','-3','-2','-1','0','1','2','3','4'});
yticklabels = ({'-3','-2','-1','0','1','2','3'});
xlabel('x-position [mm]');
ylabel('y-position [mm]');
set(gca,'xtick',xticks,'xticklabel',xticklabels)
colormap(getColorMap);

c = colorbar;  % create a colorbar
axis image;
title('Absorbed energy [norm.]');
set(gca,'FontSize',18);

%% Plot absorbed energy [J/m3]

figure;
imagesc(kgrid.y_vec*1e3, (kgrid.x_vec*1e3+3.56), vmcmedium.absorbed_energy);
hold on
plot([-2.5 2.5],[0.5 0.5],'--k',[-2.5 2.5],[0.7 0.7],'--k',[-2.5 2.5],[1 1],'--k');
xticks = ([-4 -3 -2 -1 0 1 2 3 4]);
yticks = ([-4 -3 -2 -1 0 1 2 3 4]);
xticklabels = ({'-4','-3','-2','-1','0','1','2','3','4'});
yticklabels = ({'-3','-2','-1','0','1','2','3'});
xlabel('x-position [mm]');
ylabel('y-position [mm]');
set(gca,'xtick',xticks,'xticklabel',xticklabels)
colormap(getColorMap);

c = colorbar;  % create a colorbar
axis image;
% title('Absorbed energy [J/m3]');
set(gca,'FontSize',18);
caxis([0,136]);
ylabel(c,'Absorbed energy [J/m^3]')

%% Plot Initial pressure [Pa]

figure;
imagesc(kgrid.y_vec*1e3, (kgrid.x_vec*1e3+3.56), source.p0);
hold on
plot([-2.5 2.5],[0.5 0.5],'--k',[-2.5 2.5],[0.7 0.7],'--k',[-2.5 2.5],[1 1],'--k');
xticks = ([-4 -3 -2 -1 0 1 2 3 4]);
yticks = ([-4 -3 -2 -1 0 1 2 3 4]);
xticklabels = ({'-4','-3','-2','-1','0','1','2','3','4'});
yticklabels = ({'-3','-2','-1','0','1','2','3'});
xlabel('x-position [mm]');
ylabel('y-position [mm]');
set(gca,'xtick',xticks,'xticklabel',xticklabels)
colormap(getColorMap);

c = colorbar;  % create a colorbar
axis image;
title('Initial pressure [Pa]');
set(gca,'FontSize',18);

%% Absorbed energy [norm.]

figure;
% log_data = log10(vmcmedium.absorbed_energy);
% min_data = min(min(log_data(:)));
% max_data = max(max(log_data(:)));
% norm_data = (log_data-min_data)/(max_data-min_data);

log_data = vmcmedium.absorbed_energy;
min_data = min(min(log_data(:)));
max_data = max(max(log_data(:)));
norm_data = (log_data-min_data)/(max_data-min_data);


% We have to swap x and y again 
% imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, vmcmedium.absorbed_energy);
contourf(kgrid.y_vec*1e3, (kgrid.x_vec*1e3+3.56), rot90(norm_data,2));
hold on
plot([-2.5 2.5],[0.5 0.5],'--k',[-2.5 2.5],[0.7 0.7],'--k',[-2.5 2.5],[1 1],'--k');
xticks = ([-4 -3 -2 -1 0 1 2 3 4]);
yticks = ([-4 -3 -2 -1 0 1 2 3 4]);
xticklabels = ({'-4','-3','-2','-1','0','1','2','3','4'});
yticklabels = ({'-3','-2','-1','0','1','2','3'});
xlabel('x-position [mm]');
ylabel('y-position [mm]');
set(gca,'xtick',xticks,'xticklabel',xticklabels)
colormap(getColorMap);

c = colorbar;  % create a colorbar
axis image;
title('Absorbed energy [norm.]');
set(gca,'FontSize',18);

%% Optical fluence [W/mm^2]
figure;
% imagesc(kgrid.y_vec*1e3, (kgrid.x_vec*1e3+15), rot90(solution.grid_fluence*1e3,2));
imagesc(kgrid.y_vec*1e3, (kgrid.x_vec*1e3+3.56), solution.grid_fluence*1e3);
hold on
plot([-2.5 2.5],[0.5 0.5],'--k',[-2.5 2.5],[0.7 0.7],'--k',[-2.5 2.5],[1 1],'--k');
xticks = ([-4 -3 -2 -1 0 1 2 3 4]);
yticks = ([-4 -3 -2 -1 0 1 2 3 4]);
xticklabels = ({'-4','-3','-2','-1','0','1','2','3','4'});
yticklabels = ({'-3','-2','-1','0','1','2','3'});
xlabel('x-position [mm]');
ylabel('y-position [mm]');
set(gca,'xtick',xticks,'xticklabel',xticklabels)
colormap(getColorMap);

c = colorbar;  % create a colorbar
axis image;
title('Optical fluence [W/mm^2]');
set(gca,'FontSize',18);
%% Optical fluence [W/mm^2] 2D Contour line

figure;
% imagesc(kgrid.y_vec*1e3, (kgrid.x_vec*1e3+15), rot90(solution.grid_fluence*1e3,2));
% hold on
% contourf(kgrid.y_vec*1e3, (kgrid.x_vec*1e3+15),rot90(solution.grid_fluence*1e3,2));
imagesc(kgrid.y_vec*1e3, (kgrid.x_vec*1e3+15), solution.grid_fluence*1e3);
hold on
contourf(kgrid.y_vec*1e3, (kgrid.x_vec*1e3+15),solution.grid_fluence*1e3);
hold on
plot([-15 15],[5 5],'--k',[-15 15],[6 6],'--k',[-15 15],[7 7],'--k',[-15 15],[7+depth 7+depth],'-r',[-15 15],[7+depth+diameter,7+depth+diameter],'-r');

xticks = ([-15 -10 -5 0 5 10 15]);
yticks = ([0 5 10 15 20 25 30]);
xticklabels = ({'-15','-10','-5','0','5','10','15'});
yticklabels = ({'0','5','10','15','20','25','30'});
% set(gca,'XTick',-15:5:15);
% set(gca,'YTick',0:5:30);
% set(gca,'xticklabel',{'-15','-10','-5','0','5','10','15'});
% set(gca,'yticklabel',{'0','5','10','15','20','25','30'});

xlabel('x-position [mm]');
ylabel('y-position [mm]');
set(gca,'xtick',xticks,'xticklabel',xticklabels)
colormap(getColorMap);

c = colorbar;  % create a colorbar
axis image;
title('Optical fluence [W/mm^2]');
set(gca,'FontSize',18);

%% compare different wavelength's abesorbed energy

max_absorption = zeros(41,1);
for wavelength = 700:10:1100
    filename = ['absorbenergy' '_' num2str(wavelength)  '.mat'];
    load(filename);
    eval(['absorbenergy','_',num2str(wavelength),'=','absorbed_energy',';']);
    eval(['a_midline','_',num2str(wavelength),'=','absorbenergy_',num2str(wavelength),'(:,150)',';']);
    eval(['max_absorption((',num2str(wavelength),'-700)./10+1',',1)','=','max(','a_midline','_',num2str(wavelength),'(70:end,1)',')',';']);
end

plot(700:10:1100,max_absorption);
figure;
plot(kgrid.y_vec*1e3+15,a_midline_700,LineWidth=2,Color='b');
hold on;
plot(kgrid.y_vec*1e3+15,a_midline_800,LineWidth=2,Color='k');
hold on;
plot(kgrid.y_vec*1e3+15,a_midline_900,LineWidth=2,Color='y');
hold on;
plot(kgrid.y_vec*1e3+15,a_midline_1000,LineWidth=2,Color='g');
hold on;
plot(kgrid.y_vec*1e3+15,a_midline_1100,LineWidth=2,Color='r');

grid on
xlim([0 30])
xlabel('y-position [mm]');
ylabel('Absorbed energy [J/mm3]');
legend('700nm','800nm','900nm','1000nm','1100nm')                    
title('Absorbed energy varies with depth');
set(gca,'FontSize',18);

%% compare different wavelength's initial pressure

load("initialpressure_700.mat");
initialpressure_700 = initialpressure;

load("initialpressure_800.mat");
initialpressure_800 = initialpressure;

load("initialpressure_900.mat");
initialpressure_900 = initialpressure;

load("initialpressure_1000.mat");
initialpressure_1000 = initialpressure;

load("initialpressure_1100.mat");
initialpressure_1100 = initialpressure;

midline_700 = initialpressure_700(:,150);
midline_800 = initialpressure_800(:,150);
midline_900 = initialpressure_900(:,150);
midline_1000 = initialpressure_1000(:,150);
midline_1100 = initialpressure_1100(:,150);

figure;
plot(kgrid.y_vec*1e3+15,midline_700,LineWidth=2,Color='b');
hold on;
plot(kgrid.y_vec*1e3+15,midline_800,LineWidth=2,Color='k');
hold on;
plot(kgrid.y_vec*1e3+15,midline_900,LineWidth=2,Color='y');
hold on;
plot(kgrid.y_vec*1e3+15,midline_1000,LineWidth=2,Color='g');
hold on;
plot(kgrid.y_vec*1e3+15,midline_1100,LineWidth=2,Color='r');
xticks = ([0 5 10 15 20 25 30]);
xticklabels = ({'0','5','10','15','20','25','30'});
grid on
xlim([0 30])
xlabel('y-position [mm]');
ylabel('Initial pressure [Pa]');
legend('700nm','800nm','900nm','1000nm','1100nm')                    
title('Initial pressure varies with depth');
set(gca,'FontSize',18);