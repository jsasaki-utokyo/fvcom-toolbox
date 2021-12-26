%%
load('Mobj_00.mat');

x = Mobj.lon;
y = Mobj.lat;
nodeList = double(cell2mat(Mobj.read_obc_nodes(1)));

plot_MESH(01, [Mobj.lon, Mobj.lat], Mobj.tri, Mobj.h, x(nodeList), y(nodeList));
saveas(gcf,'plot_01_mesh.png')
clear Mobj nodeList x y

%%
load('Mobj_02.mat');
load('hycom.mat');

nn = 40;   % open boundary index
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
temp = hycom.temperature.data(:, :, :, tt);

plot_TS_surface(21, ...
    hycom.lon - (dx / 2), ...
    hycom.lat - (dy / 2), ...
    squeeze(temp(:, :, hyz)), ...
    Mobj.lon(oNodes), Mobj.lat(oNodes), 40, Mobj.temperature(:, fvz, tt), ...
    lon(xidx, yidx), lat(xidx, yidx), ...
    Mobj.lon(oNodes(nn)), Mobj.lat(oNodes(nn)));
saveas(gcf,'plot_02_TS_surface.png')

plot_TS_profile(22, ...
    Mobj.temperature(nn, :, tt), Mobj.siglayz(oNodes(nn), :), ...
    squeeze(hycom.temperature.data(xidx, yidx, zidx, tt)), squeeze(-hdepth(xidx, yidx, zidx)), ...
    Mobj.temperature(nn, :, tt), 1:fz, ...
    squeeze(hycom.temperature.data(xidx, yidx, zidx, tt)), hz(zidx));
saveas(gcf,'plot_02_TS_profile.png')

vartoplot = 'temperature';
tt = 3;

lonsig(1:2,:) = repmat([1:size(oNodes,1)],2,fz);
lonsig(3:4,:) = repmat([2:(size(oNodes,1)+1)],2,fz);
depsig(1,:) = reshape(Mobj.siglayz(1:size(oNodes,1),:),[1,size(oNodes,1)*fz]);
depsig(4,:) = reshape(Mobj.siglayz(1:size(oNodes,1),:),[1,size(oNodes,1)*fz]);
depsig(2,(size(oNodes,1)+1):(size(oNodes,1)*fz)) = reshape(Mobj.siglayz(1:size(oNodes,1),1:(fz-1)),[1,size(oNodes,1)*(fz-1)]);
depsig(3,(size(oNodes,1)+1):(size(oNodes,1)*fz)) = reshape(Mobj.siglayz(1:size(oNodes,1),1:(fz-1)),[1,size(oNodes,1)*(fz-1)]);

sudata = reshape(Mobj.(vartoplot)(:,:,tt),[size(oNodes,1)*fz,1]);
DateString = datestr(Mobj.ts_times(tt,1)+678942, 31);

plot_TS_cross(23, depsig, lonsig, sudata, vartoplot,...
    min(min(min(Mobj.temperature(:,:,:)))),...
    max(max(max(Mobj.temperature(:,:,:)))),...
    DateString);

for tt=1:33
    sudata = reshape(Mobj.(vartoplot)(:,:,tt),[size(oNodes,1)*fz,1]);
    DateString = datestr(Mobj.ts_times(tt,1)+678942, 31);
    plot_TS_cross(23, depsig, lonsig, sudata, vartoplot,...
        min(min(min(Mobj.temperature(:,:,:)))),...
        max(max(max(Mobj.temperature(:,:,:)))),...
        DateString);
    pause(0.5)
end

x = 0:0.01:1;
n = 3;
y = x.^n;
plot(x,y,'LineWidth',3)
title(['y = x^n,  n = ' num2str(n) ])

fig = figure(25);
clf
for tt = 1:33
    sudata = reshape(Mobj.(vartoplot)(:,:,tt),[size(oNodes,1)*fz,1]);
    DateString = datestr(Mobj.ts_times(tt,1)+678942, 31);
    plot_TS_cross(25, depsig, lonsig, sudata, vartoplot,...
        min(min(min(Mobj.temperature(:,:,:)))),...
        max(max(max(Mobj.temperature(:,:,:)))),...
        DateString);
    drawnow
    frame = getframe(fig);
    im{tt} = frame2im(frame);
end
close;

filename = 'plot_02_TS_cross.gif'; % Specify the output file name
for tt = 1:33
    [A,map] = rgb2ind(im{tt},256);
    if tt == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end

% clear Mobj h* f* i* d* n* l* x* y* z* t* s* o* vartoplot* sudata filename

%%
load('forcing.mat');

tt=68;
% vartoplotx = 'uwnd';
% vartoploty = 'vwnd';
% [X, Y] = meshgrid(forcing.(vartoplotx).lon, forcing.(vartoploty).lat);

% for i=1:size(forcing.(vartoplot).data, 3)
%     CData.UV = sqrt(forcing.(vartoplotx).data(:, :, i).^2 + forcing.(vartoploty).data(:, :, i).^2);
%     plot_NCEP_rawdata(31, X, Y, CData.UV, 'winds (m/s)')
% end

% CData.UV = sqrt(forcing.(vartoplotx).data(:, :, tt).^2 + forcing.(vartoploty).data(:, :, tt).^2);
% plot_NCEP_rawdata(31, X, Y, CData.UV, 'winds (m/s)');

vartoplot = 'uwnd';
[X, Y] = meshgrid(forcing.(vartoplot).lon, forcing.(vartoplot).lat);
% for i=1:size(forcing.(vartoplot).data, 3)
%     CData.Var = forcing.(vartoplot).data(:, :, i);
%     plot_NCEP_rawdata(32, X, Y, CData.Var, vartoplot);
% end

CData.Var = forcing.(vartoplot).data(:, :, tt);
plot_NCEP_rawdata(32, X, Y, CData.Var, vartoplot);
saveas(gcf,'plot_03_NCEP_rawdata.png')

clear forcing X Y i CData tt vartoplot*

%%
load('Mobj_03.mat');
load('forcing.mat');
load('forcing_interp.mat');

vartoplot = 'uwnd';
tt = 068; % time index

plot_NCEP_interp(41, forcing.(vartoplot).data(:, :, tt),...
    forcing.lat, forcing.lon, forcing.(vartoplot).data(:, :, tt),...
    [Mobj.lon, Mobj.lat], Mobj.tri, forcing_interp.(vartoplot).data(:, tt),...
    vartoplot);
saveas(gcf,'plot_04_NCEP_interp.png')

% clear forcing* Mobj tt vartoplot*

