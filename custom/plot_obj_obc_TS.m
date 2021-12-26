function plot_obj_obc_TS(Num, Mobj, vartoplot)
% Number of sigma layers.
fz = size(Mobj.siglayz, 2);
oNodes = [Mobj.read_obc_nodes{:}];
% Define the name of variables
% vartoplot = {'temperature','salinity'};
% Num = 24;
[lonsig,depsig] = siglayz2xy(Mobj,oNodes);
% lonsig(1:2,:) = repmat([1:size(oNodes,1)],2,fz);
% lonsig(3:4,:) = repmat([2:(size(oNodes,1)+1)],2,fz);
% depsig(1,:) = reshape(Mobj.siglayz(1:size(oNodes,1),:),[1,size(oNodes,1)*fz]);
% depsig(4,:) = reshape(Mobj.siglayz(1:size(oNodes,1),:),[1,size(oNodes,1)*fz]);
% depsig(2,(size(oNodes,1)+1):(size(oNodes,1)*fz)) = reshape(Mobj.siglayz(1:size(oNodes,1),1:(fz-1)),[1,size(oNodes,1)*(fz-1)]);
% depsig(3,(size(oNodes,1)+1):(size(oNodes,1)*fz)) = reshape(Mobj.siglayz(1:size(oNodes,1),1:(fz-1)),[1,size(oNodes,1)*(fz-1)]);

obj.salinity = Mobj.salt;
obj.temperature = Mobj.temperature;
for i =1:length(vartoplot)
    fig = figure(Num+i);
    %clf
    fprintf([vartoplot{i},': boundary forcing.','\n']);
    minval = min(min(min(obj.(vartoplot{i})(oNodes,:,:))));
    maxval = max(max(max(obj.(vartoplot{i})(oNodes,:,:))));
    for tt = 1:length(Mobj.ts_times)
        sudata = reshape(obj.(vartoplot{i})(:,:,tt),[size(oNodes,1)*fz,1]);
        DateString = datestr(Mobj.ts_times(tt,1)+678942, 31);
        plot_TS_cross(Num+i, depsig, lonsig, sudata, vartoplot{i},'open boundary',...
            minval,...
            maxval,...
            DateString);
        drawnow
        frame = getframe(fig);
        im{tt} = frame2im(frame);
        fprintf([vartoplot{i},': ',num2str(tt),' of time ',num2str(length(Mobj.ts_times)),'\n']);
    end
    % close;
    filename = ['plot_HYCOM_obc_',vartoplot{i},'.gif']; % Specify the output file name
    for tt = 1:length(Mobj.ts_times)
        [A,map] = rgb2ind(im{tt},256);
        if tt == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1/24*6);
        end
    end
end