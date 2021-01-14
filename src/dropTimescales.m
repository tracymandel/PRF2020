function [t_e,t_r,u_lower,u_upper]= dropTimescales(t_drop,y_drop,y_ll,y_ul,y_r,d,thresh_frac,yint_offset_guess)
% [t_e,t_r,u_lower,u_upper]= dropTimescales(t_drop,y_drop,y_ll,y_ul,y_r,d)
%
% dropTimescales computes the retention and entrainment times from
% Lagrangian tracking data of an object traversing a finite transition in
% fluid density. The terminal velocities in the homogeneous-density upper
% and lower layers are also calculated.
%
% INPUTS
% t_drop:   time vector
% y_drop:   vector of vertical drop position, where z=0 at the center of
%               the stratified layer
% y_ll:     approx. highest extent of lower-layer terminal velocity
% y_ul:     approx. lowest extent of upper-layer terminal velocity
% y_r:      uppermost cutoff point to compute retention times
% d:        drop diameter
% thresh_frac:
%           fraction of drop diameter used in determining when the
%           droplet has converged to its upper-layer velocity
%           (suggested value: 0.05 or 5%)
% yint_offset_ guess:
%           initial guess for offset of y-intercept of tangent line to
%           droplet position  (suggested value: 2 to 5 cm)
%    
% OUTPUTS
% t_e:      entrainment time, i.e. time for drop to recover from
%               interactions with density transition
% t_r:      retention time, i.e. degree to which drop is delayed in
%               traversing the FOV
% u_lower:  terminal drop speed in lower layer
% u_upper:  terminal drop speed in upper layer
%
%
% This code uses the function "intersections.m" by Douglas M. Schwartz (dmschwarz@ieee.org).
% https://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections
%
% Last updated 01/14/2020 by T. Mandel (tracy.mandel@unh.edu)
%
% For more information, see:
%   Mandel, T.L., Zhou, D.Z., Waldrop, L., Theillard, M., Kleckner, D., &
%       Khatri, S. (2020). Retention of rising droplets in density
%       stratification. Physical Review Fluids 12(5): 124803.
%       https://link.aps.org/doi/10.1103/PhysRevFluids.5.124803



% i.   Find highest extent of lower layer, and compute drop terminal velocity
ind1= find(y_drop > y_ll, 1);
p= polyfit(t_drop(1:ind1),y_drop(1:ind1),1);
u_lower= p(1);

% ii.  Find lowest extent of upper layer, and compute drop terminal velocity
ind2= find(y_drop > y_ul, 1);
p= polyfit(t_drop(ind2:end),y_drop(ind2:end),1);
u_upper= p(1);


% ========== Entrainment time ==========
% t_e: Find point where drop has asymptoted to terminal upper layer velocity

% Euclidean distance between y_drop and linear fit to upper layer terminal
% velocity
mag= sqrt((y_drop - (p(1)*(t_drop)+p(2))).^2); 
mag_filt= medfilt1(mag);    % filtered Euclidean distance

% Plot this value vs. time
figure(2); clf
plot(t_drop,mag_filt,'k-','linewidth',1.3)
set(gca,'fontsize',14)
latexlabel('Time',20,'Filtered Euclidean distance',20)

% Index where the dist. between y_drop and linear fit is < some % of drop diameter
thresh= thresh_frac*d;             % threshold
hold on; plot(linspace(min(t_drop),max(t_drop),10),thresh*ones(1,10),'--','color',[.6 .6 .6],'linewidth',1.5)

% find first peak: this is initial intersection; throw points before this out
[pks,locs]= findpeaks(-mag_filt);
i1= locs(1);
hold on; plot(t_drop(locs),-pks,'b+')
plot(t_drop(i1),-pks(1),'go')

% subset of data following initial intersection
mag_sub= mag_filt(i1+1:end);
t_sub= t_drop(i1+1:end);
y_sub= y_drop(i1+1:end);

% find second local peak: throw these points out
[~,idx]= max(mag_sub);
mag_sub= mag_sub(idx:end);
t_sub= t_sub(idx:end);
y_sub= y_sub(idx:end);

% find index where this subset of the magnitude function first becomes less
% than the threshold
i_asym= find(mag_sub < thresh,1);
% if accidentally get the last point --- remove last point
if i_asym==length(mag_sub)
    mag_sub(end)= []; t_sub(end)= [];
    i_asym= find(mag_sub > thresh,1,'last');
end
hold on; plot(t_sub(i_asym),mag_sub(i_asym),'ro')

t_asym= t_sub(i_asym);  % Time at which drop reaches this asymptotic speed

% Find point in transition region where dz/dt = u_upper
%   Looking for a line with slope of u_upper that lies tangent to y_drop
%   Use fzero to solve for the y-intercept needed to have 0 intersection points
fun= @(yint) length(intersections(t_drop,u_upper*t_drop+yint,t_drop,y_drop)) - 0; % can change this to -1 to find just 1 intersection point (but does not work as well for some reason)
guess= p(2)+yint_offset_guess;       % !! adjust this input if get an error
y_int= fzero(fun,guess);

% In case the output of intersections.m is empty, gradually decrease y_int
% in increments of 0.001 until the correct y-intercept is found
[X0,Y0] = intersections(t_drop,u_upper*t_drop+y_int,t_drop,y_drop);
const= 1;
while isempty(X0)
    [X0,Y0] = intersections(t_drop,u_upper*t_drop+y_int-const*.001,t_drop,y_drop);
    const=const+1;
end
t_tan= nanmean(X0);     % Time at which drop moving at speed of u_upper
y_int= y_int-(const-1)*.001;


t_e= t_asym - t_tan;


% ========== Retention time ==========
% t_r: Find time at which drop crosses upper threshold

t= 0:.001:ceil(max(t_drop))+1;
y0= p(1)*t+y_int;   % tangent line in transition region (dot-dashed line)
ya= p(1)*t+p(2);    % line fit to upper-layer terminal behavior (dashed line)

% look at when these two lines cross an upper threshold
y_thresh_r= y_r;   
%OR
% y_thresh_r = min(max(y0),max(ya));  % whichever line ends sooner...
i0= find(y0 >= y_thresh_r,1);
ia= find(ya >= y_thresh_r,1);

t0= t(i0);
ta= t(ia);


t_r= ta-t0;


% !!!! Figures !!!!

% ENTRAINMENT TIME
figure
subplot(1,2,1)

plot(0:.1:ceil(max(t_drop))+1,p(1)*(0:.1:ceil(max(t_drop))+1)+p(2),'--','color',[.6 .6 .6],'linewidth',1.3); hold on
plot(t_tan,nanmean(Y0),'b^','markersize',7,'linewidth',1.2)
plot(t_sub(i_asym),y_sub(i_asym),'b^','markersize',7,'linewidth',1.2)
plot((0:.1:ceil(max(t_drop))+1),u_upper*(0:.1:ceil(max(t_drop))+1)+y_int,'-.','color',[.6 .6 .6],'linewidth',1.3)
plot(linspace(t_tan,t_asym,10),nanmean(Y0)*ones(1,10),'color',[.4 .4 .4],'linewidth',1.8)

plot(t_drop,y_drop,'k-','linewidth',1.6); hold on
set(gca,'fontsize',14)
xlabel('t (sec)'); ylabel('z (cm)')
title('Entrainment time')

box on

% RETENTION TIME

subplot(1,2,2)
plot(0:.1:ceil(max(t_drop))+1,p(1)*(0:.1:ceil(max(t_drop))+1)+p(2),'--','color',[.6 .6 .6],'linewidth',1.3); hold on
plot(t_tan,nanmean(Y0),'b^','markersize',7,'linewidth',1.2)
plot((0:.1:ceil(max(t_drop))+1),u_upper*(0:.1:ceil(max(t_drop))+1)+y_int,'-.','color',[.6 .6 .6],'linewidth',1.3)
plot(t0,y_thresh_r,'rs','markersize',7,'linewidth',1.2)
plot(ta,y_thresh_r,'rs','markersize',7,'linewidth',1.2)
plot(linspace(t0,ta,10),y_thresh_r*ones(1,10),'color',[.3 .3 .3],'linewidth',1.8)

plot(t_drop,y_drop,'k-','linewidth',1.6); hold on
set(gca,'fontsize',14)
xlabel('t (sec)'); ylabel('z (cm)')
title('Retention time')

box on
