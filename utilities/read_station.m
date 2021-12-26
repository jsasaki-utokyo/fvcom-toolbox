function data=read_station(header)
% This function uses the information read from the GrADS .ctl file
% to read STATION data into the matlab workspace. 
% No reference to the .map file is currently made.
% cf. http://grads.iges.org/grads/gradoc/aboutstationdata.html
% Kristof Sturm, 30.X.02

global l_quiet

fid_bin = header.FID{2};
fseek(fid_bin,0,'eof');
eof=ftell(fid_bin);
frewind(fid_bin);
file_end=0;
istat=1;

data_format={'id','name','lat','lon','var'};
data_tot={};
% Reading the header :

while 1
  
  data_loc={0,'',0,0,[]};
  
  boh=ftell(fid_bin);
  
  stid=fscanf(fid_bin,'%8c',1);
  lat=fread(fid_bin,1,'float32');
  lon=fread(fid_bin,1,'float32');
  tim=fread(fid_bin,1,'float32');
  nlev=fread(fid_bin,1,'int32');
  nflag=fread(fid_bin,1,'int32');
  
  eoh=ftell(fid_bin);
  length_h=eoh-boh;
  
  var=NaN(header.NVAR,header.TDEF.num);
  
  % Reading all variables, skipping the headers :
  for t=1:header.TDEF.num
    
    var(:,t)=fread(fid_bin,header.NVAR,'float32');
    
    if ftell(fid_bin) < (eof-2*length_h)
      file_err=fseek(fid_bin,2*length_h,'cof');
      if file_err==-1, ferror(fid_bin), end
    else
      fprintf('EOF reached for %s.\n\n',header.DATANAME);
      file_end=1;
      break
    end
    var(find(var==header.UNDEF))=NaN;
  end
  
  if file_end, break, end
  
  % Rewinding for reading the first of the next headers :
  file_err=fseek(fid_bin,-length_h,'cof');
  if file_err==-1, ferror(fid_bin), end
  
  % saving all relevant informations into the local cell :
  [data_loc]={istat,deblank(lower(stid)),lat,lon,var};
  data_tot=cat(1,data_tot,data_loc);
  istat=istat+1;
end

data=cell2struct(data_tot,data_format,2);

return