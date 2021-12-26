function data = get_GWO_forcing(time_period)
% Get the forcing data.

gwo.uwnd   = "/Users/yulong/GitHub/pipeline/01_pre_pro/00_data/2014-2017_uwnd.nc";
gwo.vwnd   = "/Users/yulong/GitHub/pipeline/01_pre_pro/00_data/2014-2017_vwnd.nc";
gwo.air    = "/Users/yulong/GitHub/pipeline/01_pre_pro/00_data/2014-2017_air.nc";
gwo.hum    = "/Users/yulong/GitHub/pipeline/01_pre_pro/00_data/2014-2017_hum.nc";
gwo.rin    = "/Users/yulong/GitHub/pipeline/01_pre_pro/00_data/2014-2017_rin.nc";
gwo.prs    = "/Users/yulong/GitHub/pipeline/01_pre_pro/00_data/2014-2017_prs.nc";
gwo.cld    = "/Users/yulong/GitHub/pipeline/01_pre_pro/00_data/2014-2017_cld.nc";
gwo.dsw    = "/Users/yulong/GitHub/pipeline/01_pre_pro/00_data/2014-2017_dsw.nc";

fields = fieldnames(gwo);

for aa = 1:length(fields)
    data.(fields{aa}).data = [];
    data.(fields{aa}).time = [];
    data.(fields{aa}).lat = [];
    data.(fields{aa}).lon = [];
    
    ncid = netcdf.open(gwo.(fields{aa}));
    
    varid = netcdf.inqVarID(ncid, 'time');
    data_time.time = netcdf.getVar(ncid, varid, 'double');   
    varid = netcdf.inqVarID(ncid,'lon');
    data_lon.lon = netcdf.getVar(ncid,varid,'double');
    varid = netcdf.inqVarID(ncid,'lat');
    data_lat.lat = netcdf.getVar(ncid,varid,'double');
    varid = netcdf.inqVarID(ncid, (fields{aa}));
    data_data.data = netcdf.getVar(ncid,varid,'double');
    
    data.(fields{aa}).time = data_time.time;
    data.(fields{aa}).lon = data_lon.lon;
    data.(fields{aa}).lat = data_lat.lat;
    data.(fields{aa}).data = data_data.data(:,:,1:length(time_period));
end

data.time = time_period;
data.lon = data.(fields{1}).lon;
data.lat = data.(fields{1}).lat;
[data.lon,data.lat] = meshgrid(data.lon,data.lat);