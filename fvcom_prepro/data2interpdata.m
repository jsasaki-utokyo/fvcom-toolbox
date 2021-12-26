function data_new = data2interpdata(data, timestep)
    % prepare data with blank
    [num_i,num_j,num_k] = size(data);
    x = 1:1:num_j;
    y = 1:1:num_i;
    z = 1:1:num_k;
    [X,Y,Z] = meshgrid(x,y,z);
    V = data;
    % figure
    % slice(X,Y,Z,V,[6 9],9,200);
    % shading flat

    [Xq,Yq,Zq] = meshgrid(1:1:num_j,1:1:num_i,1:timestep:num_k);
    Vq = interp3(X,Y,Z,V,Xq,Yq,Zq);
    % figure
    % slice(Xq,Yq,Zq,Vq,[6 9],9,200+1/24);
    % shading flat
    data_new = Vq;
end