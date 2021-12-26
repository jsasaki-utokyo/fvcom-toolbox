function url = get_url(time, threehourly)
% Child function to find the relevant URL to use for a given time step.
%
% url = get_url(time, threehourly);
%
% INPUT:
%   time - Modified Julian Day
%   threehourly - additional string to append for optional three hourly
%   outputs for the 1992/10/02 to 2008/09/18 period. Leave empty for daily
%   data.
%
% OUTPUT:
%   url - string of the approprate URL for the date supplied in time.
%

[t1, t2, t3, t4, t5, t6] = datevec(datestr(now));

if time < greg2mjulian(1992, 10, 2, 0, 0, 0)
    error('No HYCOM data available prior to 1992-10-02. Select a start date from 1992-10-02 onwards.')
elseif time >= greg2mjulian(1992, 10, 2, 0, 0, 0) && time < greg2mjulian(1995, 7, 31, 0, 0, 0)
    warning('Using the HYCOM Global Reanalysis data for dates up to 2008/09/16, thereafter the Global Analysis.')
    url = sprintf('http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.0%s?', threehourly);
elseif time >= greg2mjulian(1995, 7, 31, 0, 0, 0) && time < greg2mjulian(2008, 09, 19, 0, 0, 0)
    warning('Using the HYCOM Global Reanalysis data for dates up to 2008/09/16, thereafter the Global Analysis.')
    url = sprintf('http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1%s?', threehourly);
elseif time >= greg2mjulian(2008, 9, 19, 0, 0, 0) && time < greg2mjulian(2009, 5, 7, 0, 0, 0)
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.6?';
elseif time >= greg2mjulian(2009, 5, 7, 0, 0, 0) && time < greg2mjulian(2011, 1, 3, 0, 0, 0)
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.8?';
elseif time >= greg2mjulian(2011, 1, 3, 0, 0, 0) && time < greg2mjulian(2013, 8, 21, 0, 0, 0)
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9?';
elseif time >= greg2mjulian(2013, 8, 21, 0, 0, 0) && time < greg2mjulian(2014, 4, 21, 0, 0, 0)
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.0?';
elseif time >= greg2mjulian(2014, 4, 21, 0, 0, 0) && time < greg2mjulian(2016, 4, 18, 0, 0, 0)
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1?';
elseif time >= greg2mjulian(2016, 4, 18, 0, 0, 0) && time <= greg2mjulian(t1, t2, t3, t4, t5, t6)
    url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.2?';
elseif time > greg2mjulian(t1, t2, t3, t4, t5, t6)
    error('Given date is in the future.')
else
    error('Date is outside of the known spacetime continuum. See help TARDIS.')
end

