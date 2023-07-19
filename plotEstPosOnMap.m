function plotEstPosOnMap(pos_ecef, pos_error)

% Convert ECEF to latitude, longitude, and altitude
lla = ecef2lla(pos_ecef');

% Create a figure
figure;

% Create a colormap for the errors
colors = [0 1 0; 0.5 0 0.5; 1 0 0]; % green, purple, red
colormap(colors);

% Create a color scale based on the error
c = discretize(pos_error, [0, 1, 3, inf]);

% Plot the positions with color based on error
geoscatter(lla(:,1), lla(:,2), 10, c, 'filled');

% Add a colorbar
cb = colorbar;
cb.Ticks = [1/6, 0.5, 5/6]; % Adjust these values as needed
cb.TickLabels = {'< 1.0m', '1.0m-3m', '> 3m'};

% Set the map axes
geobasemap('satellite');