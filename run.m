clear, clc, close all;
N=128; 
T=0.9;
dt=0.001;
dtn=0.1; % The time interval for saving the snapshots
NSBFS_ETDRK3G(N, T, dt, dtn)
draw_figure(N, T, dt, dtn)
