function [val1,val2,val3] = valInPol(x,y,z,shp_path)

% find point in a shapefile polygen.

s = shaperead(shp_path);
% mapshow(s);
sx = [s.X];
sy = [s.Y];
[in_p,~] = inpolygon(x,y,sx,sy);
% plot(x(in_p),y(in_p),'g+');

val1 = x(in_p);
val2 = y(in_p);
val3 = z(in_p);