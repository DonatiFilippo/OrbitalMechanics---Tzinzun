function Planet3d(plan, position, unid, R1)

%% FUNCTION DESCRIPTION
%
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function visualizes a selected planet or celestial body in 3D, 
% using texture mapping to represent its appearance accurately. It plots 
% the selected body at the specified position, scale, and orientation.
%
%--------------------------------------------------------------------------
% INPUTS:
%   plan      [1x1]  Identifier of the planet or celestial body (integer):
%                    0:  Earth
%                    1:  Moon
%                    2:  Mercury
%                    3:  Venus
%                    4:  Mars
%                    5:  Jupiter
%                    6:  Saturn
%                    7:  Uranus
%                    8:  Neptune
%                    9:  Pluto
%                    10: Sun
%                    11: Generic algorithm (used for custom textures)
%   position  [1x3]  Position of the planet in Cartesian coordinates [km].
%   unid      [1x1]  Unit scaling factor (e.g., 1 for km, 1000 for m, etc.).
%   R1        [1x1]  Radius of the planet or celestial body [km].
%
%--------------------------------------------------------------------------
% OUTPUTS:
%   []       [figure]   Opens a figure displaying the selected planet with
%                       the applied texture.
%
%--------------------------------------------------------------------------
% Group number : 27
%
% Created and maintained by : 
%   Azevedo Da Silva Esteban
%   Gavidia Pantoja Maria Paulina
%   Donati Filippo
%   Domenichelli Eleonora
%--------------------------------------------------------------------------

switch plan
    case 0 %Earth
        R = 6371.01;                                       % [km]
        image = 'https://www.h-schmidt.net/map/map.jpg';
    case 1 %Moon
        R = 1737.4;                                       % [km]
        %image = 'https://public-files.gumroad.com/j341xwrkvrduzmdlcsiymq8wz6ev';
        image = 'https://static.wikia.nocookie.net/planet-texture-maps/images/2/22/RS3_Moon.jpg/revision/latest?cb=20220815013406';
    case 2 %Mercury
        R = 2439.7;                                       % [km]
        image = 'https://static.wikia.nocookie.net/planet-texture-maps/images/7/78/Dh_mercury_texture.png/revision/latest?cb=20211010180321';
        %image = 'https://static.wikia.nocookie.net/planet-texture-maps/images/9/9d/Mercury_Texture_Map.png/revision/latest?cb=20240821203752';
    case 3 %Venus
        R =  6051.8;                                       % [km]
        image = 'https://static.wikia.nocookie.net/planet-texture-maps/images/f/fd/Dh_venus_surface_texture.png/revision/latest?cb=20211010212031';  
        % with clouds image = 'https://static.wikia.nocookie.net/planet-texture-maps/images/0/05/Dh_venus_cloud_texture.png/revision/latest?cb=20211010212201';
    case 4 %Mars
        R = 3396.2;                                       % [km]
        image = 'https://planet-texture-maps.fandom.com/wiki/Special:FilePath/Mars_Map.png';
    case 5 %Jupiter
        R = 71492;                                       % [km]
        image = 'https://www.solarsystemscope.com/textures/download/2k_jupiter.jpg';
    case 6 %Saturn
        R = 60268;                                       % [km]
        image = 'https://static.wikia.nocookie.net/planet-texture-maps/images/7/79/Saturn-0.png/revision/latest?cb=20161015093009';   
    case 7 %Uranus
        R = 25559;                                       % [km]
        image = 'https://static.wikia.nocookie.net/planet-texture-maps/images/c/c2/Dh_uranus_texture.png/revision/latest?cb=20211011001243';     
    case 8 %Neptune
        R = 24764;                                       % [km]
        image = 'https://www.solarsystemscope.com/textures/download/2k_neptune.jpg'; 
    case 9 %Pluto
        R = 1188.3;                                       % [km]
        image = 'https://static.wikia.nocookie.net/planet-texture-maps/images/6/64/Pluto_Made.png/revision/latest?cb=20190331055010';  
    case 10 %Sun
        R = 696340;                                       % [km]
        image = 'https://www.solarsystemscope.com/textures/download/2k_sun.jpg';
    case 11 %Asteroid
        R = 100.2;                                       % [km]
        image = 'https://img.freepik.com/free-photo/close-up-various-types-mold-nature_23-2149006048.jpg';
        

end     


if exist('R1', 'var') % Verifica si R1 existe como variable
    R = R1;
end

switch unid
    case "KM"
        R = R  /5000
    case "AU"
        R = R / 149597870.7;
    case "RE"
        R = R / 6371;
    case "M"
        R = R / 1000;
    case "~"
        R = R;
    otherwise
        error("Unid not defined, must be capital letters.");
end


%% Create Earth surface as a wireframe

% Define the number of panels to be used to model the sphere 
npanels = 180;  

% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x, y, z] = ellipsoid(0, 0, 0, R, R, R, npanels);

x = x + position(1);
y = y + position(2);
z = z + position(3);

% Create the globe with the surf function
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'none');

%% Texturemap the globe

% Set the transparency of the globe: 1 = opaque, 0 = invisible
alpha = 1; 

set(globe, 'FaceColor', 'texturemap', 'CData', imread(image), 'FaceAlpha', alpha, 'EdgeColor', 'none');

if plan == 6 %If saturn, rings
% Define the radius and texture of the rings
ring_inner_radius = 66900; ring_outer_radius = 136775; 

theta = linspace(0, 2*pi, 360);
r = linspace(ring_inner_radius, ring_outer_radius, 100);

[theta_grid, r_grid] = meshgrid(theta, r);
x_ring = r_grid .* cos(theta_grid);
y_ring = r_grid .* sin(theta_grid);
z_ring = zeros(size(x_ring));
ring_texture = imread('https://www.shutterstock.com/image-illustration/saturn-ring-background-texture-3d-260nw-1856278243.jpg'); % Cambia el nombre al archivo si es necesario

hold on;
surface(x_ring, y_ring, z_ring, 'FaceColor', 'texturemap', ...
    'CData', ring_texture, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
hold off;
axis equal;
end

end

