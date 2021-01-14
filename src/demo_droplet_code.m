% demo_droplet_code.m
% =================================================
% Demo master code for calling dropTimescales.m
%
% Last updated 01/14/2020 by T. Mandel
% For more information, see:
%   Mandel, T.L., Zhou, D.Z., Waldrop, L., Theillard, M., Kleckner, D., &
%       Khatri, S. (2020). Retention of rising droplets in density
%       stratification. Physical Review Fluids 12(5): 124803.
%       https://link.aps.org/doi/10.1103/PhysRevFluids.5.124803

clear all; close all; clc

load('PRF2020_droplet_data.mat')

% Use case #15 as an example for droplet retention and entrainment time
% processing.
exp = 25;

t = drops(exp).t;       % time (seconds)
y = drops(exp).y;       % vertical position (cm)

% Can use figure to approximate y_ll and y_ul (see below)
%   (drops(exp).h_lo is also a good estimate for y_ul)
figure
plot(t,y,'b-','linewidth',1.5)
set(gca,'fontsize',14)
xlabel('Time (sec)'); ylabel('Droplet position (cm)')

%% Example code: How to use code to compute retention and entrainment times
% Other inputs to dropTimescales function
d = drops(exp).d;       % droplet diameter (cm)
y_ll = -2;              % approx. highest extent of lower-layer terminal velocity (cm)
%                           OR use pre-determined drops(exp).y_ll
y_ul = 5.6;             % approx. lowest extent of upper-layer terminal velocity (cm)
%                           OR use pre-determinedd drops(exp).y_ul
y_r = 10;               % upper treshold point to compute retention time (cm)
thresh_frac = 0.05;     % fraction of drop diameter used in determining when
                        % the droplet has converged to its upper-layer terminal velocity
yint_offset_guess = 2;  % initial guess for offset of y-intercept of tangent line to droplet position

[t_e,t_r,u_lower,u_upper]= dropTimescales(t,y,y_ll,y_ul,y_r,d,thresh_frac,yint_offset_guess);

%% Check calculations vs. processed data from Mandel et al. (2020)
figure

subplot(1,2,1)
plot(0:.1:5,0:.1:5,'k-','linewidth',1.2); hold on
plot(t_e,drops(exp).t_e,'ko','markersize',10,'linewidth',1.5)
plot(t_r,drops(exp).t_r,'rs','markersize',10,'linewidth',1.5)
set(gca,'fontsize',14)
xlabel('Time computed (sec)'); ylabel('Time from processed data (sec)')
legend('1:1 relationship','Entrainment time','Retention time')

subplot(1,2,2)
plot(0:.1:20,0:.1:20,'k-','linewidth',1.2); hold on
plot(u_upper,drops(exp).u_upper,'ko','markersize',10,'linewidth',1.5)
plot(u_lower,drops(exp).u_lower,'rs','markersize',10,'linewidth',1.5)
set(gca,'fontsize',14)
xlabel('Time computed (sec)'); ylabel('Time from processed data (sec)')
legend('1:1 relationship','U_{upper}','U_{lower}')