clc;clear;
datadir = '../EXAMPLES/2D_inplane';

% read fault data from file **_sem2d.dat
data_fault = sem2d_read_fault([datadir '/Flt01']);

% Read seismogram data from SEM2DPACK
data_seis = sem2d_read_seis();

% read spectral grid data from file grid_sem2d.dat
g = sem2d_read_specgrid(datadir);


% grid distribution
figure(1);
h=sem2d_plot_grid(g);
title('Spectral elements distribution')

% rupture_front
figure(2)
Veps = 1e-3;
Dc=0.4;
plot_fronts(data_fault.v,Veps,data_fault.d,Dc,data_fault.x,data_fault.dt)
title('Rupture front arrival time')

% % snapshot of faut
figure(3);
sem2d_snapshot_movie('vy',g,'DATADIR',datadir,'TPAUSE',1 )

% SEM seismograms
figure(4)
% Plot all x-component traces together, offset by distance to first station
d = sqrt((data_seis.x-data_seis.x(1)).^2 +(data_seis.z-data_seis.z(1)).^2);
plot_seis(d,data_seis.dt,data_seis.ux)


% % Selective uncomment to plot
% % Try plot several plots, be careful the units (stress - Pa, slip - m, slip rate - m/s)
% data_fault.t=linspace(0,data_fault.nt * data_fault.dt, data_fault.nt);
% plot(data_fault.st0) % initial shear stress distribution
% % v(nx,nt);
% plot(data_fault.t,data_fault.v(400,:)) % slip rate evolution at location x = data_fault.x(400)
% plot(data_fault.x,data_fault.v(:,data_fault.nt/2)) % slip rate snapshot at t = data_fault.t(data_fault.nt/2)
% % d(nx,nt)
% plot(data_fault.t,data_fault.d(400,:)) % same as v(nx,nt)
% plot(data_fault.x,data_fault.d(:,data_fault.nt/2)) 
