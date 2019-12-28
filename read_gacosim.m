% Function to read GACOS ztd files
% Define the full_path to GACOS zrd path :
% 'home/mgovorcin/GACOS/20190215/20190215.ztd'


function gacos = read_gacosim(infilename)

fid=fopen(infilename,'rb');
fid_header=fopen([infilename '.rsc'],'r');

header_src=textscan(fid_header,'%s%13.8f');
fclose(fid_header);

for i=1:length(header_src{1})
    name=char(header_src{1}(i));
   if(strcmp(name,'WIDTH')==1)
       width=header_src{2}(i);
       continue;
   end
   if(strcmp(name,'FILE_LENGTH')==1)
       len=header_src{2}(i);
       continue;
   end
   if(strcmp(name,'X_FIRST')==1)
       xfirst=header_src{2}(i);
   end
   if(strcmp(name,'Y_FIRST')==1)
       yfirst=header_src{2}(i);
       continue;
   end
   if(strcmp(name,'X_STEP')==1)
       xstep=header_src{2}(i);
       continue;
   end
   if(strcmp(name,'Y_STEP')==1)
       ystep=header_src{2}(i);
       continue;
   end
end

[data, count]=fread(fid,[width, len],'float');
data=data';
fclose(fid);

lat=zeros(len+width,1);
lon=zeros(len+width,1);
data_new = zeros(len+width,1);
m=1;

for i=1:len
    for j=1:width
      lat(m,1)=  yfirst +ystep*(i-1);
      lon(m,1)=  xfirst + xstep*(j-1);
      data_new(m,1) = data(i,j);
      m=m+1;
    end
end

gacos=[lon lat, data_new];