% function plotCrnt(Num, strds, lon, lat, ua, va, DateString)

% Import working dic
if ismac
    addpath('/Users/yulong/GitHub/m_map/m_map/');
elseif isunix
    addpath('/home/usr0/n70110d/github/fvcom-toolbox/m_map/');
elseif ispc
    addpath('C:/Users/Yulong WANG/Documents/GitHub/fvcom-toolbox/m_map/');
end

% filename
% ncfile = 'long_box_0001.nc';

% read file
x = ncread(ncfile,'lon');
y = ncread(ncfile,'lat');
nv = ncread(ncfile,'nv');
xc = ncread(ncfile,'lonc');
yc = ncread(ncfile,'latc');
ua = ncread(ncfile,'ua');
va = ncread(ncfile,'va');

% make matrix
xy = [x, y];
nv = double(nv);
xy = double(xy);

% make trianglua mesh
tr = triangulation(nv,xy);

% plot
figure;
clf;
for i = 1:length(xc)
    quiver(xc(i),yc(i),ua(i,25)*10,va(i,25)*10,'r-');
    hold on;
end
triplot(tr);
axis equal;
xlim([min(x)-(max(x)-min(x))*0.2 max(x)+(max(x)-min(x))*0.2]);
ylim([min(y)-(max(y)-min(y))*0.1 max(y)+(max(y)-min(y))*0.1]);
xlabel('Lon')
ylabel('Lat')

