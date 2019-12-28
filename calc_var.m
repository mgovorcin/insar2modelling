function invCov = calc_var(insar,qt_insar,ramp_flag)

% Fit experimental variogram and calculate variance-covariance matrix
% insar : struck array that contains insar-data (before or after
% corrections) Lon Lat Phase Inc Heading wavelength def_bbox
% insar.def_bbox : [UpperLon UpperLat; LowerLon Lower Lat]
% qt_insar : downsampled insar data 
% ramp_flag : 1 do not apply linear ramp correction (in case already
% appplied in gacos_tropo corr) or 0 apply linear ramp correction


% Find a local reference point
refPoint = [min(insar.Lon), min(insar.Lat)]; % Determine local reference point

% Convert phase to LOS displacement
convertedPhase = (insar.Phase / (4*pi)) * insar.wavelength;   % Convert phase from radians to m
los = single(-convertedPhase);                                % Convert to Line-of-sigth displacement in m

%% Check if the def_box exist and is not empty

    if ~isempty(insar.def_bbox)
        fprintf('Deformation are is defined, calculating the semi-variogram for the other area\n')
        def_polygon = [insar.def_bbox(1,:); [insar.def_bbox(2,1) insar.def_bbox(1,2)]; insar.def_bbox(2,:) ; [insar.def_bbox(1,1) insar.def_bbox(2,2)]; insar.def_bbox(1,:)];
        display(def_polygon)
        index=inpolygon(insar.Lon, insar.Lat, def_polygon(:,1), def_polygon(:,2));
        
        lon = insar.Lon(~index,:);
        lat = insar.Lat(~index,:);
        los = los(~index,:);
    else
        def_flag =1;
    end

if exist('def_flag')
   % Determine subsampling factor for faster plotting
if length(los) > 400000 && length(los) < 1000000
    sampling = 2;
elseif length(los) > 1000000
    sampling = 5;
else
    sampling = 1;
end

% Plot wrapped dataset
figure
scatter(insar.Lon(1:sampling:end), insar.Lat(1:sampling:end), [], mod(los(1:sampling:end), insar.wavelength/2),'.');
colormap('jet')
caxis([0 insar.wavelength/2])
axis xy
axis equal
title('Wrapped Interferogram')
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
colorbar

%% Give choice if semi-variogram should be calculated in rectangular area of interest or over the entire image with a mask

choice = questdlg('Would you like to select a rectangular area or mask out a region?', 'Option:', 'Rectangle', 'Mask','Rectangle');

switch choice
    case 'Rectangle'
        disp('Select rectangular area using mouse.')
        bounds = getrect;
        maxLon = bounds(1);
        minLon = bounds(1)+bounds(3);
        minLat = bounds(2);
        maxLat = bounds(2)+bounds(4);
        ixSubset = find(insar.Lat>minLat & insar.Lat<maxLat & insar.Lon<minLon & insar.Lon>maxLon);
    case 'Mask'
        disp('Draw closed polygon using mouse.')
        polyMask=impoly;
        pos=getPosition(polyMask);
        polyarea(pos(:,1),pos(:,2));
        in = inpolygon(insar.Lon,insar.Lat,pos(:,1),pos(:,2));
        ixSubset = find(in == 0);
end

% Extract subset from rectangular area or after masking
los = los(ixSubset);
lon = insar.Lon(ixSubset);
lat = insar.Lat(ixSubset);
end

%% Remove linear ramp if is not yet removed

sll = [lon'; lat']';
xy = llh2local(sll',refPoint);
xy = xy*1000;


if ramp_flag == 0
A = [xy' ones([length(xy) 1])];

coeff = lscov(A,los);
deramped = los - A*coeff;
end

%% Calculate and display variogram 
figure
if ramp_flag == 1
    variogDtrnd = variogram(xy',double(los),'plotit',false,'subsample',3000,'nrbins',30);
    [a,c,n] = variogramfit(variogDtrnd.distance,variogDtrnd.val,20000,1e-04,variogDtrnd.num, 'model', 'exponential', 'nugget', 1);      
    title('Semi-variogram, CORRECTED PHASE')
elseif ramp_flag ==0
    subplot (1,2,1)
    variog = variogram(xy',double(los),'plotit',true,'subsample',3000);
    title('Semi-variogram, NON-DETRENDED')

% Calculate and display variogram after detrending
    variogDtrnd = variogram(xy',double(deramped),'plotit',false,'subsample',3000,'nrbins',30);
    subplot(1,2,2)
    [a,c,n] = variogramfit(variogDtrnd.distance,variogDtrnd.val,20000,1e-04,variogDtrnd.num, 'model', 'exponential', 'nugget', 1);
    title('Semi-variogram and fit, DETRENDED')
end

fprintf('Fitted semi-variogram exponential with nugget\n')
% Print variogram exponential fit parameters to screen
disp(['Sill:   ',num2str(c)])
disp(['Range:  ',num2str(a)])
disp(['Nugget: ',num2str(n)])


%% Create inverse of covariance matrix

    qt_xy = llh2local(qt_insar(:,1:2)',refPoint)'.*1000;
    
    [X1,X2] = meshgrid(qt_xy(:,1)); % Create square matrices of Xs
    [Y1,Y2] = meshgrid(qt_xy(:,2)); % Create square matrices of Ys
    H = sqrt((X1-X2).^2 + (Y1 - Y2).^2); % Calculate distance between points
    
    covarianceMatrix = c * exp(-H/a) + n*eye(length(qt_insar)); % Calculate covariance matrix for exponential model with nugget
    invCov = inv(covarianceMatrix); % Calculate inverse of covariance matrix
    clear X1 X2 Y1 Y2 obs H covarianceMatrix
%%    




 