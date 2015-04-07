function [dots_px,still_legal] = project_dot_world(dots_m, view_window, view_dist, scale_factor)
% dots_px = project_dot_world(dots_m, view_window, view_dist, scale_factor)
%
% Takes the 3D world of dots, and projects it into two dimensions,
% according to the observer's distance from the screen and the size of the
% viewing window.
% 
% dots_m is a 3 row vector containing the [X,Y,Z] positions of the dots, in 
%       m.
% 
% view_window is the size [X,Y] of the viewing window for the projected
%       dots, in dva.
% 
% view_dist is the observer's distance from the screen, in m.
%
% scale_factor is the conversion factor between visual angle and pixels, in
%       arcmin/px.
%
% dots_px is a two row vector of the [X;Y] of each visible dot, in px. Note
% that invisble dots (those that project outside of the viewing window, or 
% are closer than the window) are not included, so it will have fewer 
% columns than dots_m.
% 
% DMG
% Updated 12/03/14

dots_2d_m = [view_dist*dots_m(1,:)./dots_m(3,:);view_dist*dots_m(2,:)./dots_m(3,:)];
dots_deg = 2*atand((dots_2d_m/2)/view_dist); %convert coordinates from m to dva

still_legalX=abs(dots_deg(1,:))<=(view_window(1)/2);
still_legalY=abs(dots_deg(2,:))<=(view_window(2)/2);
still_legalW=dots_m(3,:)>=view_dist;
still_legal = find(still_legalX.*still_legalY.*still_legalW);

dots_m(:,dots_m(3,:)<view_dist)=[]; %clip out the dots that are closer than the window
dots_deg = 2*atand((dots_2d_m/2)/view_dist); %convert coordinates from m to dva

%dots_deg(:,abs(dots_deg(1,:))>(view_window(1)/2))=[]; %clip out dots that are outside the view window (x)
%dots_deg(:,abs(dots_deg(2,:))>(view_window(2)/2))=[]; %clip out dots that are outside the view window (y)

dots_px = dots_deg*60/scale_factor.*repmat([1;-1],1,size(dots_deg,2)); %convert coordinates from dva to px. Flipping Y row, because px# goes up as you go down the screen