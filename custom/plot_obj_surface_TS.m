function plot_obj_surface_TS(Num, Mobj, vartoplot)
% Number of sigma layers.
fz = size(Mobj.siglayz, 2);
oNodes = [Mobj.read_obc_nodes{:}];
% Define the name of variables
% vartoplot = {'temperature','salinity'};
% Num = 24;

obj.salinity = Mobj.salt;
obj.temperature = Mobj.temperature;


nn = 20;   % open boundary index
tt = 1;    % time index
fvz = 1;   % fvcom depth index (currently 1-20)
hyz = 1;   % coarse depth index (1-33)
% Open boundary index.
oNodes = [Mobj.read_obc_nodes{:}];
% Number of sigma layers.
fz = size(Mobj.siglayz, 2);
fields = fieldnames(hycom);
% Find the first 4D array and use it to get the number of vertical levels
% and time steps.
for ff = 1:length(fields)
    if isfield(hycom.(fields{ff}), 'data') && ndims(hycom.(fields{ff}).data) > 3
        [nx, ny, nz, nt] = size(hycom.(fields{ff}).data);
        break
    end
end
hdepth = permute(repmat(hycom.Depth.data, [1, nx, ny]), [2, 3, 1]);
% Find the coarse seabed indices
% [~, hyz] = nanmax(hdepth, [], 3);
% Use the existing rectangular arrays for the nearest point lookup.
[lon, lat] = deal(hycom.lon, hycom.lat);
fvlon = Mobj.lon(oNodes);
fvlat = Mobj.lat(oNodes);
% Get the corresponding indices for the coarse data
[~, idx] = min(sqrt((lon(:) - fvlon(nn)).^2 + (lat(:) - fvlat(nn)).^2));
[xidx, yidx] = ind2sub(size(lon), idx);
zidx = isfinite(hdepth(xidx, yidx, :));
hz = 1:nz;
dx = mean(diff(hycom.lon(:)));
dy = mean(diff(hycom.lat(:)));
vartoplot = {'temperature','salinity'};
Mobj.salinity = Mobj.salt;
% Plot the surface temperature and salinity
for i =1:2
    var = hycom.(vartoplot{i}).data(:, :, :, tt);
    plot_TS_surface(20+i, ...
        hycom.lon - (dx / 2), ...
        hycom.lat - (dy / 2), ...
        squeeze(var(:, :, hyz)), ...
        Mobj.lon(oNodes), Mobj.lat(oNodes), 40, Mobj.(vartoplot{i})(:, fvz, tt), ...
        lon(xidx, yidx), lat(xidx, yidx), ...
        Mobj.lon(oNodes), Mobj.lat(oNodes),...
        vartoplot{i});
    name = ['plot_02_TS_surface_',vartoplot{i},'.png'];
    saveas(gcf, name)
end