function [velocity_field] = calculate_plane_flow(XYZ, centers_deg, translate, rotate, view_dist)
% velocity_field = calculate_flow_field(XYZ,translate, rotate, view_dist)
%
% The world is a plane, the observer is translating and rotating past it.
% The observer is always at the origin, and is looking down the Z axis. I'm
% not really checking the values you give me, so don't do stupid things
% like flying through the plane because weird stuff will probably happen.
%
% XYZ is the XYZ position (in m) of the locations where we are measuring
% the speed.
%
% translate tells you the observer's translation speed (Vx, Vy, Vz) in
%       m/sec
%
% rotate tells you the observer's rotation (wx, wy, wz) in rad/s.
%
% View dist is how far the observer is the from the screen in meters.
%
% velocity_field is a matrix (2x[number of locations]) reporting
% the x and y velocity in deg/s
%
% DMG
% Updated 11/24/14

if nargin==0
    [XYZ, centers_deg] = make_dot_plane(.1, 12.5, [40 32], [-20 -3 20 3]);
    translate = [0.14 0 1.9];  % the observer is walking forward at 1.9 m/s, a brisk walk
    rotate = [0 0 0];       % no rotation
    view_dist = .57;
end


velocity_field1 = zeros(size(centers_deg)); % number of patches, X and Y vel

%gonna convert the positions in deg to m, because the other measurements
%are already in m.
centers_m = view_dist*tand(centers_deg);

%solve for the inverse depth at the X,Y,Z position that projects to the
%desired x,y position, then use that to compute the flow vector at each
%location
for i=1:length(XYZ)
        inv_depth = 1/XYZ(3,i);
        
        % All these equations come from Heeger and Jepson's chapter
        x = centers_m(1,i);
        y = centers_m(2,i);
        f = view_dist;
        Pxy = inv_depth;
        
        A = [-f 0 x; 0 -f y];
        B = [(x*y)/f -(f+x^2/f) y; f+y^2/f -(x*y)/f -x];
        
        velocity_field1(1:2,i) = Pxy*A*translate'+B*rotate';
end
%gonna convert velocities in m/s to deg/sec, because thats what we actually
%care about.
velocity_field = atand(velocity_field1/view_dist);

% scatter3(XYZ(:,1),XYZ(:,3),XYZ(:,2));
% figure
% quiver(centers_deg(:,1),centers_deg(:,2),velocity_field(:,1),velocity_field(:,2));
end