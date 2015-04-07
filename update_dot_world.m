function dots = update_dot_world(dots, trans, rot, dt, nSteps)
% dots = update_dot_world(dots, trans, rot, dt, nSteps)
%
% Takes the initial array of dots, and moves them to new positions in
% observer-centered coordinates. Right now the dots have infinite
% lifetimes, but eventually this would be the place to implement death and
% rebirth if we wanted. 
% 
% dots is a 3 row vector containing the [X,Y,Z] positions of the dots, in 
%       m.
% 
% trans is a vector of the observer's velocity [X,Y,Z] in m/s.
% 
% rot is a vector of the observer's rotation [pitch, yaw, roll] in deg/s 
%       (i.e. around the X, Y, and Z axes, respectively). Positive values
%       are counter-clockwise.
%
% dt is how long between refreshes of the world, in seconds. This could
%       theoretically be as fast as the refresh rate of the monitor, but
%       usually is less frequent.
%
% nSteps is how far into the movie we are. Each time we call this function,
%       we're starting with the initial world conditions to avoid any
%       unpleasantness from the observer rotating and translating at the
%       same time (the translation is fixed in world-centered coordinates,
%       but changes in viewer-centered, if the observer is rotating).
% 
% dots the output is the new [X,Y,Z] position (in m) of all the dots,
%       viewer-centered.
% 
% DMG
% Updated 10/6/14

% First, move all the dots because of observer translation. (subtract 
% because we're effectively moving the world instead of the observer)
movement = trans*dt*nSteps; % how far is the observer moving
dots = dots-repmat(movement',1,size(dots,2)); % update the positions, 


% Second, rotate all the dots as a result of observer rotation. We're going
% to do this axis by axis, because it makes our lives easier
rot = rot*dt*nSteps; % how far is the observer rotating?

%These are rotation matrices, which move [X;Y;Z] positions to the new
%coordinate system c.f.:
% http://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
rotate_x = [1 0 0; 0 cosd(rot(1)) -sind(rot(1)); 0 sind(rot(1)) cosd(rot(1))];
rotate_y = [cosd(rot(2)) 0 sind(rot(2)); 0 1 0; -sind(rot(2)) 0 cosd(rot(2))];
rotate_z = [cosd(rot(3)) -sind(rot(3)) 0; sind(rot(3)) cosd(rot(3)) 0; 0 0 1];
rotation = rotate_x*rotate_y*rotate_z;
dots = rotation*dots;
