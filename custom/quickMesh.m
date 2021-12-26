function quickMesh(ncfile)
% Quik show the mesh of FVCOM output file.
% quickMesh(ncfile)

% filename
% ncfile = 'long_box_0001.nc';

% read file
x = ncread(ncfile,'lon');
y = ncread(ncfile,'lat');
nv = ncread(ncfile,'nv');

% make matrix
xy = [x, y];
nv = double(nv);
xy = double(xy);

% simpplot(xy,nv);

% make trianglua mesh
tr = triangulation(nv,xy);
% boundaryedges = freeBoundary(tr)'; % Bug because of islands

% plot
figure;
clf;
triplot(tr);
hold on 
% plot(xy(boundaryedges,1),xy(boundaryedges,2),'-r','LineWidth',2) % Bug because of islands
hold off
axis equal;
xlim([min(x)-(max(x)-min(x))*0.2 max(x)+(max(x)-min(x))*0.2]);
ylim([min(y)-(max(y)-min(y))*0.1 max(y)+(max(y)-min(y))*0.1]);
xlabel('Lon')
ylabel('Lat')
