function [lonsig,depsig] = siglayz2xy(Mobj,crs_node)

fz = size(Mobj.siglayz, 2);

lonsig = [];
depsig = [];
lonsig(1:2,:) = repmat([1:size(crs_node,1)],2,fz);
lonsig(3:4,:) = repmat([2:(size(crs_node,1)+1)],2,fz);
depsig(1,:) = reshape(Mobj.siglayz(crs_node,:),[1,size(crs_node,1)*fz]);
depsig(4,:) = reshape(Mobj.siglayz(crs_node,:),[1,size(crs_node,1)*fz]);
depsig(2,(size(crs_node,1)+1):(size(crs_node,1)*fz)) = reshape(Mobj.siglayz(crs_node,1:(fz-1)),[1,size(crs_node,1)*(fz-1)]);
depsig(3,(size(crs_node,1)+1):(size(crs_node,1)*fz)) = reshape(Mobj.siglayz(crs_node,1:(fz-1)),[1,size(crs_node,1)*(fz-1)]);