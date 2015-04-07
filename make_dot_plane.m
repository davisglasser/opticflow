function [dots_m, dots_deg] = make_dot_plane(density, dist, dims, exclude)
% dots = make_dot_plane(density, dist, view_dist, dims, exclude)
%
% Make a plane of random dots.
%
% density tells you the dot density in dots/deg^2 for surfaces.
%        Density from Warren & Hannon (1990) works out to about .22
%        dots/deg^2.
%
% dist tells you how far the observer starts from a dot surface in m.
%
% dims tells you the [X,Y] dimensions of the dot field on the screen in
%       deg
%
% exclude is a rectangle [X1, Y1, X2, Y2] where no dots will be placed,
%       lets you do things like block out the FOE. Can be set to 0 to
%       display the whole world.
%
% dots_m is a three row vector containing the [X;Y;Z] positions of the
%       dots, in m
%
% dots_deg is a two row vector containing the [X;Y] positions of the dots
%       on the screen, in deg
%
% DMG
% Updated 3/19/15

if nargin==0
    density = .1;
    dist = 12.5;
    dims = [40 32];
    exclude = [-20 -3 20 3];
end

nDots = ceil(prod(dims)*density);
dots_deg = [(rand(1,nDots)-.5)*dims(1); (rand(1,nDots)-.5)*dims(2)];
legal = ones(1,nDots);
if exclude
    for i=1:nDots
        if IsInRect(dots_deg(1,i),dots_deg(2,i),exclude)
            legal(i) = 0;
        end
    end
end
dots_deg = dots_deg(:,legal==1);
dots_m(1:2,:) = dist*tand(dots_deg);
dots_m(3,:) = dist;