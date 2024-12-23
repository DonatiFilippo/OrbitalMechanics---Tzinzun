function Earth_3D(Rt)

% Function to load the Earth modelled as a sphere inside a figure.
%
% Input arguments
% -------------------------------------------------------
% Rt    [1x1]       Earth mean radius       [km]
%
% ------------------------------------------------------------------------

% Set the default value for the Earth radius in case of no inputs.
if nargin < 1
    Rt = 6371.01;                                       % [km]
end

%  Load the Earth image from a website
Earth_image = 'https://cff2.earth.com/uploads/2016/10/05180920/earth-continents-map-flat_1big_stock.jpg';

% Create the figure
figure(1);
clf;
hold on;
grid on;
axis equal;
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');

% Set initial view
view(120,30);

% Define the number of panels to be used to model the sphere 
npanels = 180;  

% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x, y, z] = ellipsoid(0, 0, 0, Rt, Rt, Rt, npanels);

% Create the globe with the surf function
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'none');

% Load Earth image for texture map
cdata = imread(Earth_image);

% Set the transparency of the globe: 1 = opaque, 0 = invisible
alpha = 1; 

% Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
% specify the image data using the 'CData' property with the data loaded 
% from the image. Finally, set the transparency and remove the edges.
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');


end