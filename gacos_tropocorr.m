function corrected_phase = gacos_tropocorr(master,slave,insar,ref,planar_flag)

%% Read GACOS ztd tropospheric zenith total delay maps and calculate the tropospheric phase correction
% master : Define path to InSAR master date for the GACOS tropospheric zenith delay map
% slave : Define path to InSAR slave date for the GACOS tropospheric zenith delay map
% insar : structy array that contains the path "insar.dataPath" to the mat
% file of insar data, and sensor wavelength (insar.wavelength)
% ref : Reference phase location [Lon Lat] common for tropospheric delay
% map and insar
% planar_flag : [0/1] apply planar correction (1) of not (0)

debug_flag = 1;
wrap_flag = 1;
sampling=10;

insar_data = [insar.Lon insar.Lat insar.Phase];
wavelength =insar.wavelength;

% Read GACOS ztd image

    ztd1=read_gacosim(master);
    ztd2=read_gacosim(slave);

% Difference of ztd maps (master date- slave date) in meters
diff = ztd1(:,3)-ztd2(:,3);

ztd_diff = [ztd1(:,1:2) diff];

% Find and remove 0 or nan

index0=find(ztd_diff(:,3)==0);
ztd_diff(index0,:)=[];

indexNan=isnan(ztd_diff(:,3));
ztd_diff(indexNan,:)=[];


% if debug_flag==1
%     subplot(1,3,1)
%     scatter(ztd1(:,1),ztd1(:,2),[],ztd1(:,3),'.')
%     colorbar
%     title('Master ZTD')
%     subplot(1,3,2)
%     scatter(ztd1(:,1),ztd1(:,2),[],ztd1(:,3),'.')
%     colorbar
%     title('Slave ZTD')
%     subplot(1,3,3)
%     scatter(ztd_diff(:,1),ztd_diff(:,2),[],ztd_diff(:,3),'.')
%     cb1=colorbar;
%     title(cb1,'[m]');
%     title('ZDT difference')
% end

%% Get ZTD_diff extent
ztd_bbox = [min(ztd_diff(:,1)) max(ztd_diff(:,1)) min(ztd_diff(:,2)) max(ztd_diff(:,2))];
fprintf('Tropospheric zenith delay difference map extent\n')
fprintf('Min Lon: %.3f , Max Lon: %.3f\n',ztd_bbox(1,1),ztd_bbox(1,2))
fprintf('Min Lat: %.3f , Max Lat: %.3f\n',ztd_bbox(1,3),ztd_bbox(1,4))

%GET INSAR extent
insar_bbox = [min(insar_data(:,1)) max(insar_data(:,1)) min(insar_data(:,2)) max(insar_data(:,2))];
fprintf('INSAR data extent\n')
fprintf('Min Lon: %.3f , Max Lon: %.3f\n',insar_bbox(1,1),insar_bbox(1,2))
fprintf('Min Lat: %.3f , Max Lat: %.3f\n',insar_bbox(1,3),insar_bbox(1,4))

%Check if GACOS diff map covers InSAR data

if ztd_bbox(:,1) > insar_bbox(:,1)
    fprintf('Download new GACOS ztd maps extended to West: first Lon %.3f\n',insar_bbox(:,1)-0.001)
    bbox(1,1) = ztd_bbox(:,1); % use GACOS min Lon for futher processing, cropping InSAR data
else
    bbox(1,1) = insar_bbox(:,1); % use INSAR min Lon for futher processing
end

if ztd_bbox(:,2) < insar_bbox(:,2)
    fprintf('Download new GACOS ztd maps extended to East: last Lon %.3f\n',insar_bbox(:,2)+0.001)
    bbox(1,2) = ztd_bbox(:,2); % use GACOS max Lon for futher processing, cropping InSAR data
else
    bbox(1,2) = insar_bbox(:,2); % use INSAR max Lon for futher processing
end

if ztd_bbox(:,3) > insar_bbox(:,3)
    fprintf('Download new GACOS ztd maps extended to South: first Lat %.3f\n',insar_bbox(:,3)-0.001)
    bbox(1,3) = ztd_bbox(:,3); % use GACOS min Lat for futher processing, cropping InSAR data
else
    bbox(1,3) = insar_bbox(:,3); % use INSAR min Lat for futher processing
end

if ztd_bbox(:,4) < insar_bbox(:,4)
    fprintf('Download new GACOS ztd maps extended to North: last Lat %.3f\n',insar_bbox(:,4)+0.001)
    bbox(1,4) = ztd_bbox(:,4); % use GACOS max Lat for futher processing, cropping InSAR data
else
    bbox(1,4) = insar_bbox(:,4); % use INSAR max Lat for futher processing
end

%% Crop data to the defined extent

index1 = find(insar_data(:,1)>bbox(:,1) & insar_data(:,1)<bbox(:,2) & insar_data(:,2)>bbox(:,3) & insar_data(:,2)<bbox(:,4));
index2 = find(ztd_diff(:,1)>bbox(:,1) & ztd_diff(:,1)<bbox(:,2) & ztd_diff(:,2)>bbox(:,3) & ztd_diff(:,2)<bbox(:,4));

insar_data = insar_data(index1,:);
inc = insar.Inc(index1,:);
head = insar.Heading(index1,:);
ztd_diff = ztd_diff(index2,:);

% Interpolate ZTD data to InSAR location
fprintf('Interpolating zentih tropospheric phase delay difference data to InSAR points\n')
Fdata=scatteredInterpolant(ztd_diff(:,1),ztd_diff(:,2), ztd_diff(:,3));
Fdata.Method='nearest'; %nearest neighborhood
Fdata.ExtrapolationMethod='nearest'; %nearest neighborhood
Interpolate=Fdata(insar_data(:,1),insar_data(:,2));

% if debug_flag==1
%     figure
%     scatter(insar_data(:,1),insar_data(:,2),[],Interpolate,'.')
%     cb2=colorbar;
%     title('Slant-range delay diff')
%     title(cb2,'[m]');
% end


%% Calculate the tropospheric phase delay in satellite slant range

scale = (-4*pi)/(wavelength);

ztd_diff_slant = (Interpolate.* scale)./cosd(inc);


%% Set common reference point
ref_phase=find(insar_data(:,1)>ref(1,1)-0.0005 & insar_data(:,1)<ref(1,1)+0.0005 & insar_data(:,2)>ref(1,2)-0.0005 & insar_data(:,2)<ref(1,2)+0.0005);

%Check if there ref point is within the data extent

if isempty(ref_phase)
    fprintf('########################################################\n')
    fprintf('ERROR: Reference point is outside of the data extent!!!\n')
    fprintf('Define new location of Reference point\n')
    figure
    scatter(insar_data(:,1),insar_data(:,2),[],insar_data(:,3),'.')
    colormap('jet')
    hold on
    rr = plot(ref(1,1),ref(1,2),'r*');
    legend(rr,'Current reference point')
    return
end

seed_ztd=ztd_diff_slant(ref_phase,1);
seed_phase=insar_data(ref_phase,3);
fprintf('Mean InSAR seed point value: %.4f\n',mean(seed_phase))
fprintf('Mean ZTD seed point value: %.4f\n',mean(seed_ztd))
fprintf('Calculating slant-range tropospheric phase delay difference\n')
ztd_diff_slant(:,1)=ztd_diff_slant(:,1)-mean(seed_ztd);
insar_data(:,3)=insar_data(:,3)-mean(seed_phase);
corr_phase = insar_data(:,3) - ztd_diff_slant;


%% Remove planar ramp
%Remove planar trend of an interferogram
%The planar is defined as:
%          Trend = a0 + a1 * X + a2 * Y

if planar_flag ==1
%conver to local coordinate 
fprintf('Deramping .....\n')
xy = llh2local(insar_data(:,1:2)',[min(insar_data(:,1)) min(insar_data(:,2))])'.*1000;

% p = mean([xy(:,1:2), corr_phase],1);
% 
% [~,~,V]=svd(bsxfun(@minus,[xy(:,1:2), corr_phase],p),0);
% 
% P = V(:,end);
% 
% planeEquation=[P;-dot(P,p(:))];
% 
% planar = (planeEquation(1)*xy(:,1) + planeEquation(2)*xy(:,2) + planeEquation(4) )./(-planeEquation(3));
% 
%Define A matrix

A = [xy ones(length(corr_phase),1)];

% Compute linear relation
coeff =lscov(A,corr_phase);

planar = A*coeff;

corr_planar_phase = corr_phase - planar;

n_corr = 2;
 else
n_corr = 1;
    

end

if wrap_flag==1
   pl_insar_data =  mod((insar_data(:,3)./(-4*pi))* wavelength,wavelength/2);
   pl_ztd_diff_slant = mod((ztd_diff_slant(:,1)./(-4*pi))*wavelength,wavelength/2);
   pl_corr_phase = mod((corr_phase(:,1)./(-4*pi))*wavelength,wavelength/2);
   if planar_flag==1
       pl_planar = mod((planar(:,1)./(-4*pi))*wavelength,wavelength/2);
       pl_corr_planar_phase = mod((corr_planar_phase./(-4*pi))*wavelength,wavelength/2);
   end
else
   pl_insar_data =  (insar_data(:,3)./(-4*pi))*wavelength;
   pl_ztd_diff_slant = (ztd_diff_slant(:,1)./(-4*pi))*wavelength;
   pl_corr_phase = (corr_phase(:,1)./(-4*pi))* wavelength;
   if planar_flag==1
       pl_planar = (planar(:,1)./(-4*pi))* wavelength;
       pl_corr_planar_phase = (corr_planar_phase./(-4*pi))* wavelength;
   end
    
end



%%
if debug_flag==1
    fprintf('Creating figure .....\n')
    tropo_fig=figure('Position',[0 0 2000 600],'Visible','on');
    subplot(2,n_corr+1,1)
    scatter(insar_data(1:sampling:end,1),insar_data(1:sampling:end,2),[],pl_insar_data(1:sampling:end),'.')
    title('InSAR phase')
    axis equal
    axis tight
    colorbar
    subplot(2,n_corr+1,2)
    scatter(insar_data(1:sampling:end,1),insar_data(1:sampling:end,2),[],pl_ztd_diff_slant(1:sampling:end),'.')
    title('SR tropo phase delay')
    colorbar
    axis equal
    axis tight
    
    if planar_flag==1 
            subplot(2,n_corr+1,3)
            scatter(insar_data(1:sampling:end,1),insar_data(1:sampling:end,2),[],pl_planar(1:sampling:end),'.')
            title('Linear phase ramp')
            colorbar
            axis equal
            axis tight
            
            subplot(2,n_corr+1,4)
            scatter(insar_data(1:sampling:end,1),insar_data(1:sampling:end,2),[],pl_insar_data(1:sampling:end),'.')
            title('InSAR phase')
            axis equal
            axis tight
            colorbar
            
            subplot(2,n_corr+1,5)
            scatter(insar_data(1:sampling:end,1),insar_data(1:sampling:end,2),[],pl_corr_phase(1:sampling:end),'.')
            title('InSAR phase - SR tropo')
            axis equal
            axis tight
            colorbar
            
            subplot(2,n_corr+1,6)
            scatter(insar_data(1:sampling:end,1),insar_data(1:sampling:end,2),[],pl_corr_planar_phase(1:sampling:end),'.')
            title('InSAR phase - SR-tropo - Linear ramp')
            axis equal
            axis tight
            colorbar
            
    else
    
            subplot(2,n_corr+1,3)
            scatter(insar_data(1:sampling:end,1),insar_data(1:sampling:end,2),[],pl_insar_data(1:sampling:end),'.')
            title('InSAR phase')
            axis equal
            axis tight
            colorbar
            
            subplot(2,n_corr+1,4)
            scatter(insar_data(1:sampling:end,1),insar_data(1:sampling:end,2),[],pl_corr_phase(1:sampling:end),'.')
            title('IFG phase - SR tropo')
            axis equal
            axis tight
            colorbar
            
    end
    colormap('jet')
    
    %picture resolution
    if ~exist('./tropo','dir')
        mkdir('tropo')
    end
      
    rl = 300; 
    tropo_fig.PaperUnits='centimeters';
    tropo_fig.PaperPosition = [ 0 0 20 8];
    saveas(tropo_fig,['./tropo/',insar.name,'_tropo.png']);
    fprintf('Saving figure as %s\n',['./tropo/',insar.name,'_tropo.png'])
    %print(tropo_fig,[insar.name,'_tropo.png'],'-depsc2',sprintf('-r%d',rl))
end

if planar_flag == 1
    Tdelay = corr_planar_phase;
else
    Tdelay = [insar_data(:,1) insar_data(:,2) corr_phase inc head];
end

corrected_phase.Lon = insar_data(:,1);
corrected_phase.Lat = insar_data(:,2);
corrected_phase.Phase = Tdelay;
corrected_phase.Inc = inc;
corrected_phase.Heading = head;
corrected_phase.wavelength = wavelength;
corrected_phase.name = insar.name;
end






