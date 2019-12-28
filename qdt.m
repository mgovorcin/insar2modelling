function down_insar = qdt(insar, mindim, maxdim, thresh)
% Quadtree decomposition of insar_data
% insar : struct array that contain Lon Lat Phase Inc Heading wavelength
% mindim : minimum size of the quadtree cell (32)
% maxdim : maximum size of the quadtree cell (256)
% thresh : variance threshold (0.04)
% modified after Diego Melgar, MudPy quadtree 
h = waitbar(0,'Quadtree downsample');

debug_flag=1;

plot_lims=[-0.1,0.1];

%Load InSAR
lon=double(insar.Lon);
lat=double(insar.Lat);
los=double(insar.Phase);
inc=insar.Inc;
head=insar.Heading;
wavelength = insar.wavelength;

convertedPhase = (los / (4*pi) )  * wavelength;    % Convert phase from radians to m
los = double(-convertedPhase);  % Convert to Line-of-sigth displacement in m

% Calculate three dimensions INSAR displacements
UE = -cosd(head).* sind(inc);
UN = sind(head).* sind(inc);
UV = cosd(inc);

lookE=los.*UE;
lookN=los.*UN;
lookU=los.*UV;

%Prepare grid
dlon=8.34e-4;
dlat=8.34e-4;
waitbar(0.1,h,'Quadtree downsample');
lon_i=linspace(min(lon),max(lon),2048);
lat_i=linspace(min(lat),max(lat),2048);

% Get the InSAR boundary polygon
tic
fprintf('Getting the border of INSAR data\n')
border=boundary(lon,lat);
toc
waitbar(0.3,h,'Quadtree downsample');
tic
fprintf('Interpolating to grid\n')
[X,Y]=meshgrid(lon_i,lat_i);

los_interp_sharp = griddata(lon,lat,los,X,Y);
los_interp=los_interp_sharp;
toc
waitbar(0.5,h,'Quadtree downsample');

display('QT decomp')
s=qtdecomp(los_interp,thresh,[mindim,maxdim]);
%Calculate maximum dimension
maxdim=max(max(full(s)));
%get power of 2 of possible block values
idim=(log(mindim)/log(2)):1:(log(maxdim)/log(2));
%Now loop through
los_out=[];
lon_out=[];
lat_out=[];

tic
for k=1:length(idim)
    [vals, r, c] = qtgetblk(los_interp_sharp, s,2^idim(k));
    icurrent=find(s==2^idim(k));
    %Now get mean and cell center of each grid
    for kgrid=1:length(icurrent)
       %Get values
       new_los=mean(mean(vals(:,:,kgrid)));
       if ~isnan(new_los)  %Add to the list
           los_out=[los_out,mean(mean(vals(:,:,kgrid)))];
           %get indices of center as upepr left plus half the size of the
           %block
           r1=r(kgrid)+(2^idim(k))/2;
           r2=r(kgrid)+(2^idim(k))/2+1;
           c1=c(kgrid)+(2^idim(k))/2;
           c2=c(kgrid)+(2^idim(k))/2+1;
           %Now figure out coordiantes of the center
           lon_out=[lon_out,0.5*(X(r1,c1)+X(r2,c2))];
           lat_out=[lat_out,0.5*(Y(r1,c1)+Y(r2,c2))];
       end
    end
end
toc
waitbar(0.7,h,'Quadtree downsample');
fprintf('Find look direction vector\n')
%Find closest look direction vector
lookE_out=[];
lookN_out=[];
lookU_out=[];
tic
for k=1:length(lon_out)
   d=sqrt((lon_out(k)-lon).^2+(lat_out(k)-lat).^2);
   [a,i]=min(d);
   lookE_out=[lookE_out,lookE(i)];
   lookN_out=[lookN_out,lookN(i)];
   lookU_out=[lookU_out,lookU(i)];
end
toc
%find points withn the InSAR boundaries
index = inpolygon(lon_out,lat_out,lon(border),lat(border));
lon_out = lon_out(index);
lat_out = lat_out(index);
los_out = los_out(index);
lookE_out=lookE_out(index);
lookN_out=lookN_out(index);
lookU_out=lookU_out(index);

waitbar(0.8,h,'Quadtree downsample');

if debug_flag==1
%Plot
qd_fig=figure;
subplot(1,2,1)
mesh(X,Y,los_interp_sharp)
colorbar
view([0,90])
hold on
axis tight
axis equal
scatter3(lon_out,lat_out,100000*ones(size(lat_out)),'kx')
plot(lon(border),lat(border),'r-')
xlabel('Longitude','FontSize',16)
ylabel('Latitude','FontSize',16)

title('LOS (m)','FontSize',16)
caxis(plot_lims)
colormap('jet')
set(gca,'FontSize',14)
subplot(1,2,2)
scatter(lon_out,lat_out,40,los_out,'filled')
colorbar
axis tight
axis equal
view([0,90])
caxis(plot_lims)
colormap('jet')
title('QuadTree Resampled LOS (m)','FontSize',16)
saveas(qd_fig,[insar.name,'_qt.png'])
end
fprintf('Quadtree finished\n')
fprintf('Decrease the point sample size from %g to %g\n',length(lon),length(lon_out))
down_insar = [lon_out' lat_out' los_out' lookE_out' lookN_out' lookU_out'];
close(h)
end
