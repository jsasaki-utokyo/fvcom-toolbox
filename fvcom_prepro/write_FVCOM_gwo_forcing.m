function write_FVCOM_gwo_forcing(Mobj, data, fileprefix, infos, fver, varargin)
% Write data out to FVCOM netCDF forcing file.
%
% write_FVCOM_forcing(Mobj, fvcom_forcing_file, data, infos, fver)
%
% DESCRIPTION:
%   Takes the given interpolated data (e.g. from grid2fvcom) and writes out
%   to a netCDF file.
%
% INPUT:
%   Mobj - MATLAB mesh object containing fields:
%       tri - triangulation table for the unstructured grid
%       nVerts - number of grid vertices (nodes)
%       nElems - number of grid elements
%       nativeCoords - model coordinate type ('cartesian' or 'spherical')
%       x, y or lon, lat - node positions (depending on nativeCoords value)
%       gwo - gwo data add by Yulong
%   fileprefix - Output netCDF file prefix (plus path) will be
%       fileprefix_{wnd,hfx,evap}.nc if fver is '3.1.0', otherwise output
%       files will be fileprefix_wnd.nc.
%   % data - Struct of the data to be written out.
%   infos - Additional remarks to be written to the "infos" netCDF variable
%   fver - Output for version 3.1.0 or 3.1.6. The latter means all the
%       forcing can go in a single file, the former needs separate files
%       for specific forcing data (wind, heating and precipitation).
%   Optional keyword-argument pairs. These control the time variables. This
%   script defaults to writing 'Times' only.
%   FVCOM needs only one of:
%       1. Times: character string of times
%       2. Itime and Itime2: integer days and milliseconds since midnight
%       3. time: float days.
%   FVCOM checks for these in the order above and this script defaults to
%   writing Times only. Adjust the keyword-argument pairs to your liking:
%
%   'strtime' = set to true to output the 'Times' variable
%   'inttime' = set to true to output the 'Itime' and 'Itime2' variables
%   'floattime' = set to true to output the 'time' variable
%
% The fields in data may be called any of:
%     - 'u10', 'v10' or 'uwnd', 'vwnd' - wind components
%     - 'Et' or 'evap'      - evaporation
%     - 'prate' or 'P_E'    - precipitation
%     - 'nlwrs'             - net longwave radiation*,**
%     - 'nswrs'             - net shortwave radiation*,**,***
%     - 'shtfl'             - sensible heat net flux*,**
%     - 'lhtfl'             - latent heat net flux*,**
%     - 'slp' or 'pres'     - mean sea level pressure***
%     - 'dswrf'             - downward shortwave flux
%     - 'dlwrf'             - downward longwave flux***
%     - 'rhum'              - relative humidity***
%     - 'air'               - air temperature***
%     - 'lon'               - longitude (vector)
%     - 'lat'               - latitude (vector)
%     - 'x'                 - eastings (vector)
%     - 'y'                 - northings (vector)
%     - 'nshf'              - pre-computed net surface heat flux**
%     - 'lcc'               - low cloud cover
%
% Fields marked with an * are combined to form the "surface net heat flux"
% (nshf) as follows:
%
%   nshf = nlwrs         + nswrs         - lhtfl - shtfl;
%   nshf = dlwrf - ulwrf + dswrf - uswrf - lhtfl - shtfl;
%
% ** Alternatively, a new field 'nshf' (net surface heat flux) can be
% supplied, in which case shtfl and lhtfl are not necessary as their only
% function is in the calculation of net surface heat flux. This approach
% eliminates the need to interpolate so many variables onto the FVCOM grid,
% thus decreasing the time needed to generate the FVCOM input files.
%
% *** These fields are required for HEATING_CALCULATED model input files.
%
% OUTPUT:
%   FVCOM forcing netCDF file(s)
%
% EXAMPLE USAGE:
%   windBase = '/path/to/output/casename_wnd.nc';
%   write_FVCOM_forcing(Mobj, windBase, data, ...
%       'FVCOM atmospheric forcing data', '3.1.6');
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%   Karen Amoudry (National Oceanography Centre, Liverpool)
%   Rory O'Hara Murray (Marine Scotland Science)
%
% PWC Revision history:
%   2012-11-01 - First version based on the parts of grid2fvcom_U10V10.m
%   which dealt with writing the netCDF file. This version now dynamically
%   deals with varying numbers of forcing data.
%   2012-11-09 - Add the correct calculation for the net surface heat flux.
%   2013-02-14 - Fix the surface net heat flux sign convention, include the
%   longwave radiation variable and add support for a new field in the
%   input struct ('nshf') which contains the pre-computed net surface heat
%   flux.
%   2013-05-13 - Fix the evaporation to use the correct variable from NCEP
%   (pevpr rather than P_E which is actually the precipitation minus the
%   evaporation in Et). The data in Et are calcaulated from lhtfl whereas
%   pevpr comes directly from NCEP and to me it seems more sensible to use
%   that to maintain consistency.
%   2013-05-14 - Add example usage to the help and specify which fields are
%   required in Mobj.
%   2013-06-21 - Remove support for pevpr (pevpr is in W/m^{2} from the
%   NCEP Reanalysis 2 data (FVCOM wants evaporation in m/s). Update the
%   help accordingly.
%   2013-10-24 - Add support for writing all the variables needed for a
%   HEATING_CALCULATED model run. This essentially makes
%   write_FVCOM_heating redundant, but I'll leave it in as it's a bit
%   simpler to understand what's going on there. Also update the way the
%   net surface heat flux is calculated (sum long and short, subtract
%   latent and sensible).
%   2014-01-08 - Fix the way the wind variables are handled so that both
%   the U10/V10 and uwind_speed/vwind_speed netCDF variables are written if
%   only one of data.u10/data.v10 or data.uwnd/data.vwnd is given.
%   2014-08-11 - (PWC) Add new flags to control which time variables to
%   use. FVCOM reads the 'Times' variable first if present, then falls back
%   to 'Itime' and 'Itime2' and finally 'time'.
%
% KJA Revision history:
%   2013-01-16 - Added support for output of sea level pressure.
%   2013-08-16 - Updated output of Itime2 to avoid rounding errors
%   when converting from double to single format.
%   2013-09-03 - Removed PWC's fix for timestrings. Issue was due to
%   rounding errors caused by mjulian2greg.m, which have now been fixed.
%   2013-09-26 - Added support for output of low cloud cover and specific
%   humidity.
%
% ROM Revision History:
%   2013-12-02 Change the output of tri' to tri, as tri was being written
%   the wrong way around.
%
%==========================================================================

multi_out = false; % default to 3.1.6, single output file
if (nargin < 4 || nargin > 5) && isempty(varargin)
    error('Incorrect number of arguments')
elseif nargin == 5
    if strcmpi(fver, '3.1.0')
        multi_out = true;
    end
end

subname = 'write_FVCOM_forcing';

global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

% Default to string times as FVCOM looks for these first.
% Add time dependent option at 2019-11-21, which means spatial constant.
strtime = false;
inttime = false;
floattime = false;
julian = false;
time_dependent = false;
for vv = 1:2:length(varargin)
    switch varargin{vv}
        case 'strtime'
            strtime = true;
        case 'inttime'
            inttime = true;
        case 'floattime'
            floattime = true;
        case 'julian'
            julian = true;
        case 'time dependent'
            time_dependent = true;
    end
end


% data = Mobj.gwo;
tri = Mobj.tri;
nNodes = Mobj.nVerts;
nElems = Mobj.nElems;
ntimes = numel(data.time);

% if strcmpi(Mobj.nativeCoords, 'cartesian')
%     x = Mobj.x;
%     y = Mobj.y;
% else
%     x = Mobj.lon;
%     y = Mobj.lat;
% end
% Create a string for each variable's coordinate attribute
coordString = sprintf('FVCOM %s coordinates', Mobj.nativeCoords);

% Create element vertices positions
% xc = nodes2elems(x, Mobj);
% yc = nodes2elems(y, Mobj);
xc = Mobj.xc;
yc = Mobj.yc;
x = Mobj.x;
y = Mobj.y;
lonc = Mobj.lonc;
latc = Mobj.latc;
lon = Mobj.lon;
lat = Mobj.lat;
%--------------------------------------------------------------------------
% Create the netCDF header for the FVCOM forcing file
%--------------------------------------------------------------------------

if multi_out
    suffixes = {'_wnd', '_hfx', '_evap', '_air_press'};
else
    suffixes = {'_wnd'};
end

% We use this variable to indicate whether we were given precalculated net
% surface heat flux.
nshf = 0;

for i=1:length(suffixes)
    nc = netcdf.create([fileprefix, suffixes{i}, '.nc'], 'clobber');
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type','FVCOM Forcing File')
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title','FVCOM Forcing File')
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'institution','Plymouth Marine Laboratory')
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'institution','Sasaki Lab, The University of Tokyo')
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'source','FVCOM grid (unstructured) surface forcing')
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history', sprintf('File created with %s from the MATLAB fvcom-toolbox', subname))
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'references','http://fvcom.smast.umassd.edu, http://codfish.smast.umassd.edu')
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.0')
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'infos',infos)
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'CoordinateSystem',Mobj.nativeCoords)
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'CoordinateProjection','init=WGS84') % WGS84?epsg:4326?

    % Dimensions
    nele_dimid=netcdf.defDim(nc,'nele',nElems);
    node_dimid=netcdf.defDim(nc,'node',nNodes);
    three_dimid=netcdf.defDim(nc,'three',3);
    time_dimid=netcdf.defDim(nc,'time',netcdf.getConstant('NC_UNLIMITED'));
    datestrlen_dimid=netcdf.defDim(nc,'DateStrLen',26);

    % Space variables
    x_varid=netcdf.defVar(nc,'x','NC_DOUBLE',node_dimid);
    netcdf.putAtt(nc,x_varid,'long_name','nodal x-coordinate');
    netcdf.putAtt(nc,x_varid,'units','meters');

    y_varid=netcdf.defVar(nc,'y','NC_DOUBLE',node_dimid);
    netcdf.putAtt(nc,y_varid,'long_name','nodal y-coordinate');
    netcdf.putAtt(nc,y_varid,'units','meters');
    
    lon_varid=netcdf.defVar(nc,'lon','NC_DOUBLE',node_dimid);
    netcdf.putAtt(nc,lon_varid,'long_name','nodal longitude');
    netcdf.putAtt(nc,lon_varid,'units','degrees');

    lat_varid=netcdf.defVar(nc,'lat','NC_DOUBLE',node_dimid);
    netcdf.putAtt(nc,lat_varid,'long_name','nodal latitude');
    netcdf.putAtt(nc,lat_varid,'units','degrees');

    xc_varid=netcdf.defVar(nc,'xc','NC_DOUBLE',nele_dimid);
    netcdf.putAtt(nc,xc_varid,'long_name','zonal x-coordinate');
    netcdf.putAtt(nc,xc_varid,'units','meters');

    yc_varid=netcdf.defVar(nc,'yc','NC_DOUBLE',nele_dimid);
    netcdf.putAtt(nc,yc_varid,'long_name','zonal y-coordinate');
    netcdf.putAtt(nc,yc_varid,'units','meters');

    lonc_varid=netcdf.defVar(nc,'lonc','NC_DOUBLE',nele_dimid);
    netcdf.putAtt(nc,lonc_varid,'long_name','zonal longitude');
    netcdf.putAtt(nc,lonc_varid,'units','degrees');

    latc_varid=netcdf.defVar(nc,'latc','NC_DOUBLE',nele_dimid);
    netcdf.putAtt(nc,latc_varid,'long_name','zonal latitude');
    netcdf.putAtt(nc,latc_varid,'units','degrees');
    
    nv_varid=netcdf.defVar(nc,'nv','NC_INT',[nele_dimid, three_dimid]);
    netcdf.putAtt(nc,nv_varid,'long_name','nodes surrounding element');

%    Time variables
    if floattime
       time_varid=netcdf.defVar(nc,'time','NC_FLOAT',time_dimid);
       netcdf.putAtt(nc,time_varid,'long_name','time');
       if julian
           netcdf.putAtt(nc,time_varid,'units','days since 1858-11-17 00:00:00');
           netcdf.putAtt(nc,time_varid,'format','modified julian day (MJD)');
           netcdf.putAtt(nc,time_varid,'time_zone','UTC');
       else
           netcdf.putAtt(nc,time_varid,'units','days since 0.0');
           netcdf.putAtt(nc,time_varid,'time_zone','none');
       end
    end

    if inttime
       itime_varid=netcdf.defVar(nc,'Itime','NC_INT',time_dimid);
       if julian
           netcdf.putAtt(nc,itime_varid,'units','days since 1858-11-17 00:00:00');
           netcdf.putAtt(nc,itime_varid,'format','modified julian day (MJD)');
           netcdf.putAtt(nc,itime_varid,'time_zone','UTC');
       else
           netcdf.putAtt(nc,itime_varid,'units','days since 0.0');
           netcdf.putAtt(nc,itime_varid,'time_zone','none');
       end
       itime2_varid=netcdf.defVar(nc,'Itime2','NC_INT',time_dimid);
       netcdf.putAtt(nc,itime2_varid,'units','msec since 00:00:00');
       netcdf.putAtt(nc,itime2_varid,'time_zone','UTC');
    end

    if strtime
       times_varid=netcdf.defVar(nc,'Times','NC_CHAR',[datestrlen_dimid,time_dimid]);
       netcdf.putAtt(nc,times_varid,'long_name','Calendar Date');
       netcdf.putAtt(nc,times_varid,'format','String: Calendar Time');
       netcdf.putAtt(nc,times_varid,'time_zone','UTC');
    end

    % Since we have a dynamic number of variables in the struct, try to be
    % a bit clever about how to create the output variables.
    fnames = fieldnames(data);
    used_varids = cell(0);
    used_fnames = cell(0);
    used_dims = cell(0); % exclude time (assume all variables vary in time)

    for vv=1:length(fnames)
        % Need to check both whether we have the data but also whether
        % we're outputting to several netCDF files. If so, we drop some
        % variables if we're in the wrong file loop.
        switch fnames{vv}
            case {'uwnd', 'u10'}
                if strcmpi(suffixes{i}, '_wnd') || ~multi_out
                    % wind components (assume we have v if we have u)
                    % On the elements
%                     u10_varid=netcdf.defVar(nc,'U10','NC_FLOAT',[nele_dimid, time_dimid]);
%                     netcdf.putAtt(nc,u10_varid,'long_name','Eastward Wind Speed');
%                     netcdf.putAtt(nc,u10_varid,'units','m/s');
%                     netcdf.putAtt(nc,u10_varid,'grid','fvcom_grid');
%                     netcdf.putAtt(nc,u10_varid,'coordinates',coordString);
%                     netcdf.putAtt(nc,u10_varid,'type','data');
%                     v10_varid=netcdf.defVar(nc,'V10','NC_FLOAT',[nele_dimid, time_dimid]);
%                     netcdf.putAtt(nc,v10_varid,'long_name','Northward Wind Speed');
%                     netcdf.putAtt(nc,v10_varid,'units','m/s');
%                     netcdf.putAtt(nc,v10_varid,'grid','fvcom_grid');
%                     netcdf.putAtt(nc,v10_varid,'coordinates',coordString);
%                     netcdf.putAtt(nc,v10_varid,'type','data');

                    uwind_varid=netcdf.defVar(nc,'uwind_speed','NC_FLOAT',[nele_dimid, time_dimid]);
                    netcdf.putAtt(nc,uwind_varid,'long_name','Eastward Wind Speed');
                    netcdf.putAtt(nc,uwind_varid,'standard_name','Wind Speed');
                    netcdf.putAtt(nc,uwind_varid,'units','m/s');
                    netcdf.putAtt(nc,uwind_varid,'grid','fvcom_grid');
                    netcdf.putAtt(nc,uwind_varid,'type','data');
                    vwind_varid=netcdf.defVar(nc,'vwind_speed','NC_FLOAT',[nele_dimid, time_dimid]);
                    netcdf.putAtt(nc,vwind_varid,'long_name','Northward Wind Speed');
                    netcdf.putAtt(nc,vwind_varid,'standard_name','Wind Speed');
                    netcdf.putAtt(nc,vwind_varid,'units','m/s');
                    netcdf.putAtt(nc,vwind_varid,'grid','fvcom_grid');
                    netcdf.putAtt(nc,vwind_varid,'type','data');

                    % Only on the elements (both U10/V10 and uwind_speed
                    % and vwind_speed).
                    % used_varids = [used_varids, {'u10_varid', 'v10_varid', 'uwind_varid', 'vwind_varid'}];
                    used_varids = [used_varids, {'uwind_varid', 'vwind_varid'}];
                    % We should only have one of {u,v}wnd or {u,v}10 as the
                    % field name and used_fnames needs to reflect that.
                    if isfield(data, 'uwnd') && ~isfield(data, 'u10')
                        % We have uwnd and vwnd variants (no u10/v10).
                        % used_fnames = [used_fnames, {'uwnd', 'vwnd', 'uwnd', 'vwnd'}];
                        used_fnames = [used_fnames, {'uwnd', 'vwnd'}];
                    elseif isfield(data, 'u10') && ~isfield(data, 'uwnd')
                        % We have only u10 and v10 variants (no uwnd/vwnd)
                        % used_fnames = [used_fnames, {'u10', 'v10', 'u10', 'v10'}];
                        used_fnames = [used_fnames, {'u10', 'v10'}];
                    elseif isfield(data, 'u10') && isfield(data, 'uwnd')
                        error('Supply only one set of wind fields: ''uwnd'' and ''vwnd'' or ''u10'' and ''v10''.')
                    else
                        error('Unrecognised wind field names.')
                    end
                    % used_dims = [used_dims, {'nElems', 'nElems', 'nElems', 'nElems'}];
                    used_dims = [used_dims, {'nElems', 'nElems'}];
                end

            case {'vwnd', 'v10'}
                % We dealt with these in the u component section above, so
                % just pass silently.
                true;

            case {'slp', 'pres', 'pressfc','prs'}
                if strcmpi(suffixes{i}, '_wnd') || ~multi_out
                    % Sea level pressure
                    slp_varid=netcdf.defVar(nc,'air_pressure','NC_FLOAT',[node_dimid, time_dimid]);
                    netcdf.putAtt(nc,slp_varid,'long_name','Surface air pressure');
                    netcdf.putAtt(nc,slp_varid,'units','hPa');
                    netcdf.putAtt(nc,slp_varid,'grid','fvcom_grid');
                    netcdf.putAtt(nc,slp_varid,'coordinates',coordString);
                    netcdf.putAtt(nc,slp_varid,'type','data');
                    used_varids = [used_varids, 'slp_varid'];
                    used_fnames = [used_fnames, fnames{vv}];
                    used_dims = [used_dims, 'nNodes'];
                end

            case {'lcc','cld'}'
                if strcmpi(suffixes{i}, '_wnd') || ~multi_out
                    % Cloud cover
                    cld_varid=netcdf.defVar(nc,'cloud_cover','NC_FLOAT',[node_dimid, time_dimid]);
                    netcdf.putAtt(nc,cld_varid,'long_name','Cloud cover');
                    netcdf.putAtt(nc,cld_varid,'units','-');
                    netcdf.putAtt(nc,cld_varid,'grid','fvcom_grid');
                    netcdf.putAtt(nc,cld_varid,'coordinates',coordString);
                    netcdf.putAtt(nc,cld_varid,'type','data');
                    used_varids = [used_varids, 'cld_varid'];
                    used_fnames = [used_fnames, fnames{vv}];
                    used_dims = [used_dims, 'nNodes'];
                end

            case {'prate','rin'}%'prate'
                if strcmpi(suffixes{i}, '_wnd') || ~multi_out
                    prate_varid=netcdf.defVar(nc,'Precipitation','NC_FLOAT',[node_dimid, time_dimid]);
                    netcdf.putAtt(nc,prate_varid,'long_name','Precipitation');
                    netcdf.putAtt(nc,prate_varid,'description','Precipitation, ocean lose water is negative');
                    netcdf.putAtt(nc,prate_varid,'units','m s-1');
                    netcdf.putAtt(nc,prate_varid,'grid','fvcom_grid');
                    netcdf.putAtt(nc,prate_varid,'coordinates',coordString);
                    netcdf.putAtt(nc,prate_varid,'type','data');
                    used_varids = [used_varids, 'prate_varid'];
                    used_fnames = [used_fnames, fnames{vv}];
                    used_dims = [used_dims, 'nNodes'];
                end

            case {'air', 'tmp2m'}
                if strcmpi(suffixes{i}, '_wnd') || ~multi_out
                    % Air temperature.
                    airt_varid = netcdf.defVar(nc, 'air_temperature', 'NC_FLOAT', [node_dimid, time_dimid]);
                    netcdf.putAtt(nc, airt_varid, 'long_name', 'Surface air temperature');
                    netcdf.putAtt(nc, airt_varid, 'units', 'Celsius Degree');
                    netcdf.putAtt(nc, airt_varid, 'grid', 'fvcom_grid');
                    netcdf.putAtt(nc, airt_varid, 'coordinates', coordString);
                    netcdf.putAtt(nc, airt_varid, 'type', 'data');
                    used_varids = [used_varids, 'airt_varid'];
                    used_fnames = [used_fnames, fnames{vv}];
                    used_dims = [used_dims, 'nNodes'];
                end

            case {'rhum','hum'}
                if strcmpi(suffixes{i}, '_wnd') || ~multi_out
                    % Relative humidity
                    rhum_varid = netcdf.defVar(nc, 'relative_humidity', 'NC_FLOAT', [node_dimid, time_dimid]);
                    netcdf.putAtt(nc, rhum_varid, 'long_name', 'surface air relative humidity');
                    netcdf.putAtt(nc, rhum_varid, 'units', 'percentage');
                    netcdf.putAtt(nc, rhum_varid, 'grid', 'fvcom_grid');
                    netcdf.putAtt(nc, rhum_varid, 'coordinates', coordString);
                    netcdf.putAtt(nc, rhum_varid, 'type', 'data');
                    used_varids = [used_varids, 'rhum_varid'];
                    used_fnames = [used_fnames, fnames{vv}];
                    used_dims = [used_dims, 'nNodes'];
                end

%             % Downward solar wave radation
            case {'dswrf', 'dswsfc','dsw'}
                if strcmpi(suffixes{i}, '_wnd') || ~multi_out
                    % Downward shortwave radiation
                    dswrf_varid = netcdf.defVar(nc, 'short_wave', 'NC_FLOAT', [node_dimid, time_dimid]);
                    netcdf.putAtt(nc, dswrf_varid, 'long_name', 'Downward solar shortwave radiation flux');
                    netcdf.putAtt(nc, dswrf_varid, 'units', 'Watts meter-2');
                    netcdf.putAtt(nc, dswrf_varid, 'grid', 'fvcom_grid');
                    netcdf.putAtt(nc, dswrf_varid, 'coordinates', coordString);
                    netcdf.putAtt(nc, dswrf_varid, 'type', 'data');
                    used_varids = [used_varids, 'dswrf_varid'];
                    used_fnames = [used_fnames, fnames{vv}];
                    used_dims = [used_dims, 'nNodes'];
                end
% 
            % case {'time', 'lon', 'lat', 'x', 'y'}
            case {'time'}
                continue

            otherwise
                if ftbverbose
                    warning('Unknown or possibly unused input data type: %s', fnames{vv})
                end
        end
    end

    % End definitions
    netcdf.endDef(nc);

    % Put the easy ones in first.
    netcdf.putVar(nc,x_varid,x);
    netcdf.putVar(nc,y_varid,y);
    netcdf.putVar(nc,xc_varid,xc);
    netcdf.putVar(nc,yc_varid,yc);
    netcdf.putVar(nc,lon_varid,lon);
    netcdf.putVar(nc,lat_varid,lat);
    netcdf.putVar(nc,lonc_varid,lonc);
    netcdf.putVar(nc,latc_varid,latc);
    netcdf.putVar(nc, nv_varid, tri);
    % Do the times.
    if floattime
        if julian
            netcdf.putVar(nc,time_varid,0,ntimes,data.time);
        else
            netcdf.putVar(nc,time_varid,0,ntimes,data.time-data.time(1));
        end
    end

    if inttime
        if julian
            netcdf.putVar(nc,itime_varid,0,ntimes,floor(data.time));
        else
            netcdf.putVar(nc,itime_varid,0,ntimes,floor(data.time-data.time(1)));
        end
        % etcdf.putVar(nc,itime2_varid,0,ntimes,mod(data.time,1)*24*3600*1000); % PWC original
        % KJA edit: avoids rounding errors when converting from double to single
        % Rounds to nearest multiple of the number of msecs in an hour
        netcdf.putVar(nc,itime2_varid,0,ntimes,round((mod(data.time,1)*24*3600*1000)/(3600*1000))*(3600*1000));
    end

    if strtime
        nStringOut = char();
        [nYr, nMon, nDay, nHour, nMin, nSec] = mjulian2greg(data.time);
        for i=1:ntimes
            nDate = [nYr(i), nMon(i), nDay(i), nHour(i), nMin(i), nSec(i)];
            nStringOut = [nStringOut, sprintf('%04i/%02i/%02i %02i:%02i:%09.6f', nDate)];
        end
        netcdf.putVar(nc,times_varid,[0, 0], [26, ntimes], nStringOut);
    end

    % Now do the dynamic ones. Set the heat flux to not done (hf_done = 0)
    % until we hit one of the holy quad (shtfl, lhtfl, nlwrs and nswrs).
    % Also make sure we have wind data, either as u10/v10 or uwnd/vwnd.
    hf_done = 0;
    wnd_done = 0;
    for ff = 1:length(used_fnames)
        % if ftbverbose
            fprintf('write : %s... \n', used_fnames{ff})
        % end
%         if strcmpi(used_fnames{ff}, 'shtfl') || strcmpi(used_fnames{ff}, 'lhtfl') || strcmpi(used_fnames{ff}, 'nlwrs') || strcmpi(used_fnames{ff}, 'nswrs')
% 
%             hf_done = hf_done + 1;
% 
%             if hf_done == 4 && nshf == 0
%                 % if ftbverbose
%                     fprintf('combining heat flux ... \n')
%                 % end
%                 % We've got all four heat parameters, so dump them into the
%                 % file.
%                 %hf = -(data.shtfl.node + data.lhtfl.node + ...
%                 %    data.nlwrs.node + data.nswrs.node);
%                 hf = data.nlwrs.node + data.nswrs.node - ...
%                     data.shtfl.node - data.lhtfl.node;
%                 %%% netcdf.putVar(ncid,varid,start,count,data)
%                 %%% https://jp.mathworks.com/help/matlab/ref/netcdf.putvar.html?lang=zh
%                 if time_dependent
%                     start_row = 1:nNodes;
%                     for ii = 1:length(start_row)
%                         netcdf.putVar(nc,nhf_varid,[start_row(ii)-1,0],[1,ntimes],full(hf));
%                     end
%                 else
%                     netcdf.putVar(nc,nhf_varid,[0,0],[nNodes,ntimes],full(hf));
%                 end
%                 clear hf;
%             elseif strcmpi(used_fnames{ff}, 'nswrs') || strcmpi(used_fnames{ff}, 'nlwrs')
%                 % We've already done the net surface heat flux but we're on
%                 % either of the other fluxes (short/long wave) which we
%                 % need to dump. Do that here.
%                 if strcmpi(used_dims{ff}, 'nNodes')
%                     if time_dependent
%                         % eval(['netcdf.putVar(nc,',used_varids{ff},',[0,0],[1,ntimes],full(data.',used_fnames{ff},'.node));'])
%                         start_row = 1:nNodes;
%                         for ii = 1:length(start_row)
%                             eval(['netcdf.putVar(nc,',used_varids{ff},',[',num2str(start_row(ii)-1),',0],[1,ntimes],full(data.',used_fnames{ff},'.node));'])
%                         end
%                     else
%                         eval(['netcdf.putVar(nc,',used_varids{ff},',[0,0],[',used_dims{ff},',ntimes],full(data.',used_fnames{ff},'.node));'])
%                     end
%                 else
%                     if time_dependent
%                         % eval(['netcdf.putVar(nc,',used_varids{ff},',[0,0],[1,ntimes],full(data.',used_fnames{ff},'.data));'])
%                         start_row = 1:nElems;
%                         for ii = 1:length(start_row)
%                             eval(['netcdf.putVar(nc,',used_varids{ff},',[',num2str(start_row(ii)-1),',0],[1,ntimes],full(data.',used_fnames{ff},'.data));'])
%                         end
%                     else
%                         eval(['netcdf.putVar(nc,',used_varids{ff},',[0,0],[',used_dims{ff},',ntimes],full(data.',used_fnames{ff},'.data));'])
%                     end
%                 end
%             else
%                 % We haven't got the precomputed net surface heat flux but
%                 % we haven't yet got enough of the parameters to export the
%                 % heat flux from the short + long + latent + sensible.
%                 % Essentially this loop just does hf_done = hf_done + 1.
%             end
%         elseif strcmpi(used_fnames{ff}, 'nswsfc') && nshf == 1
%             if ftbverbose
%                 fprintf('existing combined heat flux ... ')
%             end
%             % We have pre-computed net surface heat flux, in which case set
%             % hf_done to 4 and put the data into the netCDF. Also set the
%             % nshf variable 1 to stop the net surface heat flux variable
%             % being overwritten above.
%             hf_done = 4;
%             if time_dependent
%                 % netcdf.putVar(nc, nhf_varid, [0, 0], [1, ntimes],full(data.nshf.node))
%                 start_row = 1:nNodes;
%                 for ii = 1:length(start_row)
%                     netcdf.putVar(nc, nhf_varid, [start_row(ii)-1, 0], [1, ntimes], full(data.nshf.node))
%                 end
%             else
%                 netcdf.putVar(nc, nhf_varid, [0, 0], [nNodes, ntimes], full(data.nshf.node))
%             end
%         else
            % One of the other data sets for which we can simply dump the
            % existing array without waiting for other data.
            if strcmpi(used_dims{ff}, 'nNodes')
                if time_dependent
                    % eval(['netcdf.putVar(nc,',used_varids{ff},',[0,0],[1,ntimes],full(data.',used_fnames{ff},'.node));'])
                    start_row = 1:nNodes;
                    for ii = 1:length(start_row)
                        eval(['netcdf.putVar(nc,',used_varids{ff},',[',num2str(start_row(ii)-1),',0],[1,ntimes],full(data.',used_fnames{ff},'.node));'])
                    end
                else
                    eval(['netcdf.putVar(nc,',used_varids{ff},',[0,0],[',used_dims{ff},',ntimes],full(data.',used_fnames{ff},'.node));'])
                end
            else
                try
                    if time_dependent
                        % eval(['netcdf.putVar(nc,',used_varids{ff},',[0,0],[1,ntimes],full(data.',used_fnames{ff},'.data));'])
                        start_row = 1:nElems;
                        for ii = 1:length(start_row)
                            eval(['netcdf.putVar(nc,',used_varids{ff},',[',num2str(start_row(ii)-1),',0],[1,ntimes],full(data.',used_fnames{ff},'.data));'])
                        end
                    else
                        eval(['netcdf.putVar(nc,',used_varids{ff},',[0,0],[',used_dims{ff},',ntimes],full(data.',used_fnames{ff},'.data));'])
                    end
%                     wnd_done = wnd_done + 1;
                catch err
                    fprintf('%s', err.message)
                end
            end
%         end
        if ftbverbose
            fprintf('done.\n')
        end
    end
%     if hf_done < 4 && nshf == 1
%         % hf_done might be higher than four, but unless it is at least
%         % four, we haven't got everything we need. Only trigger this
%         % warning if we've been given any of the net heat flux components.
%         warning('Did not have all the required heat flux parameters for HEATING_ON. Need ''shtfl'', ''lhtfl'', ''nlwrs'' and ''nwsrs''.')
%     end

%     if wnd_done < 2
%         warning('No wind data was provided (or one component was missing). Expected fields u10 and v10 or uwnd and vwnd.')
%     end

    % Close the netCDF file(s)
    netcdf.close(nc);
end

if ftbverbose
   fprintf('end   : %s\n', subname)
end