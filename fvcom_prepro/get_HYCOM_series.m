function Mobj = get_HYCOM_series(Mobj, timestep, varargin)

for v = 1:2:length(varargin)
    switch varargin{v}
        case 'temperature'
            obs_temp = Mobj.temperature;
            obs_time = Mobj.ts_times;
        case 'salinity'
            obs_salt = Mobj.salt;
            obs_time = Mobj.ts_times;
        case 'u'
            obs_u = Mobj.u;
            obs_time = Mobj.mf_times;
        case 'v'
            obs_v = Mobj.v;
            obs_time = Mobj.mf_times;
        otherwise
            error('Unrecognised ''source'' type. Valid values are ''temperature'', ''salinity'', ''u'' and ''v''.')
    end
end

% fix coninuity problem of the time
for i = 1:length(obs_time)-1
    obs_time(i+1) = obs_time(1) + (obs_time(2)-obs_time(1)) * i;
end

% timestep = 1/24;
ratio = (obs_time(2)-obs_time(1))/timestep;
dims_old = length(obs_time);
dims_new = (dims_old-1)*ratio;

% Time
for i = 1:ratio-1
    obs_time(:,i+1) = obs_time(:,1) + 1/24 * i;
end
obs_time = reshape(obs_time',[],1);
obs_time = obs_time(1:dims_new,1);

for v = 1:2:length(varargin)
    switch varargin{v}
        case 'temperature'
            Mobj.temperature = data2interpdata(obs_temp, timestep);
            Mobj.temperature = Mobj.temperature(:,:,1:end-1);
            Mobj.ts_times = obs_time;
        case 'salinity'
            Mobj.salt = data2interpdata(obs_salt, timestep);
            Mobj.salt = Mobj.salt(:,:,1:end-1);
            Mobj.ts_times = obs_time;
        case 'u'
            Mobj.u = data2interpdata(obs_u, timestep);
            Mobj.u = Mobj.u(:,:,1:end-1);
            Mobj.mf_times = obs_time;
        case 'v'
            Mobj.v = data2interpdata(obs_v, timestep);
            Mobj.v = Mobj.v(:,:,1:end-1);
            Mobj.mf_times = obs_time;
        otherwise
            error('Unrecognised ''source'' type. Valid values are ''temperature'', ''salinity'', ''u'' and ''v''.')
    end
end

end

