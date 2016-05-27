clear all,clc,close all;

% plot snapshots 
%addpath (genpath('/home/yuan/Desktop/GPX/Geophysics_3.0'));

NX=400;
NZ=200;


ntstep = 1500;
nstep= nrecv*ntstep;

filename_snap = '/gpfs/scratch/yuan/gim_2D/sou_dipole/test_dipole/2d_start.vs';

 
snap_file = fopen (filename_snap);
data_snap = fread(snap_file, 'real*4');
fclose('all');


data_snap=reshape(data_snap,NX,NZ);

figure
imagesc(data_snap');
% colorbar
% caxis([1300 1800])


% mesh(data_snap)
% pcolor(flipud(data_snap'))
% shading interp
% 

