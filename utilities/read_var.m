function data=read_var(header,ivar,t_limits,z_limits,y_limits,x_limits)
% data=read_varseq(header,ivar) reads in 'ivar' from
% 'header.DATANAME' and return it as a 4-D matlab matrix. 
% 'ivar' can be the variable''s short name or its number id.

global l_quiet

vars=header.VARS;
nvar=header.NVAR;
var_size=header.VARSIZE{:};

% different behaviour of the f* functions according to the chosen data precision :
switch lower(header.BINPRECISION)
 case {'uchar','schar','int8','uint8'}
  byte_length=1;
 case {'int16','uint16'}
  byte_length=2;  
 case {'float32','single','int32','uint32'}
  byte_length=4;
 case {'float64','double','int64','uint64'}
  byte_length=8;
 otherwise
  byte_length=1;
end

if isempty(ivar)
  if ~l_quiet
    disp(['Header for ',header.FILENAME,' read.'])
  end
  data=[];
  return
end

if ischar(ivar)
  var_name=ivar;
  ivar=strmatch(ivar,{vars.name},'exact');
  if isempty(ivar)
    if ~l_quiet
      disp(['The variable ',var_name,' cannot be found in ',header.DATANAME])
    end
    data=[];
    return
  end
  clear var_name
end

if ivar==1
  byte_ini=0;
else
  byte_ini = sum(prod(var_size(1:(ivar-1),1:3),2),1) ;
end

if isfield(header,'SEQ') && header.SEQ
  byte_ini=byte_ini+(sum(var_size(1:ivar-1,3))*8+4)/byte_length;
end

if header.FILEHEADER
  byte_ini = byte_ini + header.FILEHEADER/byte_length;
% $$$   if isfield(header,'SEQ') && header.SEQ
% $$$     byte_ini=byte_ini+8/byte_length;
% $$$   end
end

if header.THEADER
  byte_ini = byte_ini + header.THEADER/byte_length;
% $$$   if isfield(header,'SEQ') && header.SEQ
% $$$     byte_ini=byte_ini+8/byte_length;
% $$$   end
end

if header.XYHEADER
  byte_ini=byte_ini+(sum(var_size(1:ivar-1,3))+1)*header.XYHEADER/byte_length;
% $$$   if isfield(header,'SEQ') && header.SEQ
% $$$     byte_ini=byte_ini+(sum(var_size(1:ivar-1,3))+1)*8/byte_length;
% $$$   end
end

block_length = sum(prod(var_size(:,1:3),2),1);

if isfield(header,'SEQ') && header.SEQ
 block_length=block_length+sum(var_size(:,3))*8/byte_length;
end

if header.XYHEADER
  block_length=block_length+sum(var_size(:,3))*header.XYHEADER/byte_length;
end

if header.THEADER
  block_length=block_length+header.THEADER/byte_length;
end

slice_length = header.XDEF.num*header.YDEF.num+header.XYHEADER/byte_length;

if isfield(header,'SEQ') && header.SEQ
  slice_length = slice_length + 8/byte_length ;
end

var_length = prod([header.XDEF.num header.YDEF.num diff(z_limits)+1]);

if header.XYHEADER
  var_length = var_length + (diff(z_limits)+1)*header.XYHEADER/byte_length;
end

if isfield(header,'SEQ') && header.SEQ
  var_length=var_length+(diff(z_limits)+1)*8/byte_length;
end

% disp(byte_ini)
% Selecting the starting time slice :
byte_ini=byte_ini + (t_limits(1)-1)*block_length + (z_limits(1)-1)*slice_length ;

err_stat = fseek(header.FID{2},byte_ini*byte_length,'bof');
if err_stat == -1, ferror(header.FID{2}), end

% data=NaN(diff(x_limits)+1,diff(y_limits)+1,diff(z_limits)+1,diff(t_limits)+1);
data=NaN*ones(diff(x_limits)+1,diff(y_limits)+1,diff(z_limits)+1,diff(t_limits)+1);

%   size(data)
count_byte = 0;
count_byte2 = 0;
for l=1:(t_limits(2)-t_limits(1)+1)
  file_ind = ftell(header.FID{2});
  for k=1:(z_limits(2)-z_limits(1)+1)
    if all([x_limits y_limits]==[1 header.XDEF.num 1 header.YDEF.num])
    try
      [data(:,:,k,l),n]=fread(header.FID{2},[header.XDEF.num header.YDEF.num],header.BINPRECISION);
      count_byte=count_byte+n;
    catch
      warning(['Data read-in interrupted at l=',int2str(l),' and k=',int2str(k)])
      break
    end
    else
    % Advance to the j-th row:
    fseek(header.FID{2},(y_limits(1)-1)*header.XDEF.num*byte_length,'cof');
    for j=1:(diff(y_limits)+1)
      fseek(header.FID{2},(x_limits(1)-1)*byte_length,'cof');
      [data(:,j,k,l),n]=fread(header.FID{2},diff(x_limits)+1,header.BINPRECISION);
      fseek(header.FID{2},(header.XDEF.num-x_limits(2))*byte_length,'cof');
    end
    fseek(header.FID{2},(header.YDEF.num-y_limits(2))*header.XDEF.num*byte_length,'cof');
    end
    if isfield(header,'SEQ') && header.SEQ
      fseek(header.FID{2},8,'cof');
    end
    if header.XYHEADER
      fseek(header.FID{2},header.XYHEADER,'cof');
    end
  end
  count_byte2=count_byte2+ftell(header.FID{2})-file_ind;
  if l < header.TDEF.num
    err_stat=fseek(header.FID{2},(block_length-var_length)*byte_length,'cof');
    if err_stat == -1, ferror(header.FID{2}), end
  end
end

% disp([byte_ini block_length slice_length var_length file_ind/byte_length ftell(header.FID{2})/byte_length])

data(find(data==header.UNDEF))=NaN;   % replacing the missing/undefined values by NaN (not-a-number)

if header.XDEF.rev, data=flipdim(data,1); end
if header.YDEF.rev, data=flipdim(data,2); end
if header.ZDEF.rev, data=flipdim(data,3); end

%if count_byte ~= count_byte2/byte_length
%  warning(['Number of elements read : ',num2str(count_byte),...
%	   ' .NEQ. number of bytes advanced in ',header.DATANAME,': ',num2str(count_byte2/byte_length)])
%end

%fprintf('%i bytes were read in %s .\n size(%s) = [%s], i.e. %i elements.\n',...
%        count_byte2,header.DATANAME,vars(ivar).name,int2str(size(data)),count_byte);
%	count_byte2,header.DATANAME,vars(ivar).name,int2str(size(data)),prod(size(data)));

return