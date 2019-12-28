% INPUT FILE

%% Reference point and Area of interest
refPoint = [19.619 41.88];    % Longitude and Latitude in degrees for arbitrary reference point of local coordinates system [Lon; Lat;]
boundingBox = [19.1; 42; 20.5; 41;];

tropocorrection_step = 1;    %Flag to do step 2 and apply tropospheric phase correction [1] yes [0] no
covariance_step = 0;         %Flag to do step 4 and calculate the covariance matrix [1] yes [0] no
quadtree_step = 0;           %Flag to do step 3 and perform quadtree downsample [1] yes [0] no


%% InSAR data
% Make sure insarID is unique!
% 
% insarID = 1;                            % InSAR dataset unique identifier
% 
% insar{insarID}.dataPath = '/media/marin/marin/3_Working_DIR/Tirana/MODEL/DATA/A_T73_20191125_20191201.mat'; % Path to data file
% insar{insarID}.wavelength = 0.056;      % Wavelength in m (e.g., Envisat/ERS/Sentinel: 0.056; CSK/TSX/TDX: 0.031)
% insar{insarID}.quadtreeThresh = 0.04; % Quadtree threshold (e.g., 0.04 )
% insar{insarID}.quadtreeMindim = 32; % Quadtree threshold (e.g., 0.04 )
% insar{insarID}.quadtreeMaxdim = 256; % Quadtree threshold (e.g., 0.04 )
% insar{insarID}.gacosmaster = '/media/marin/marin/3_Working_DIR/Tirana/GACOS/A73/20191125.ztd'; % Path to GACOS ztd file on master date
% insar{insarID}.gacosslave = '/media/marin/marin/3_Working_DIR/Tirana/GACOS/A73/20191201.ztd'; % Path to GACOS ztd file on slave date
% insar{insarID}.ramp_flag = 0;
% insar{insarID}.def_bbox = [19.35 41.65;19.65 41.25];
% insar{insarID}.name = 'Tirana_t_73_asc';

% insarID = 1;                            % InSAR dataset unique identifier
% 
% insar{insarID}.dataPath = '/media/marin/marin/3_Working_DIR/Tirana/MODEL/DATA/A_T175_20191114_20191126.mat'; % Path to data file
% insar{insarID}.wavelength = 0.056;      % Wavelength in m (e.g., Envisat/ERS/Sentinel: 0.056; CSK/TSX/TDX: 0.031)
% insar{insarID}.quadtreeThresh = 0.04; % Quadtree threshold (e.g., 0.04 )
% insar{insarID}.quadtreeMindim = 32; % Quadtree threshold (e.g., 0.04 )
% insar{insarID}.quadtreeMaxdim = 256; % Quadtree threshold (e.g., 0.04 )
% insar{insarID}.gacosmaster = '/media/marin/marin/3_Working_DIR/Tirana/GACOS/A175/20191114.ztd'; % Path to GACOS ztd file on master date
% insar{insarID}.gacosslave = '/media/marin/marin/3_Working_DIR/Tirana/GACOS/A175/20191126.ztd'; % Path to GACOS ztd file on slave date
% insar{insarID}.ramp_flag = 0;
% insar{insarID}.def_bbox = [19.35 41.65;19.65 41.25];
% insar{insarID}.name = 'Tirana_t_175_14_26_asc';


% insarID = 1;                            % InSAR dataset unique identifier
% 
% insar{insarID}.dataPath = '/media/marin/marin/3_Working_DIR/Tirana/MODEL/DATA/A_T175_20191120_20191126.mat'; % Path to data file
% insar{insarID}.wavelength = 0.056;      % Wavelength in m (e.g., Envisat/ERS/Sentinel: 0.056; CSK/TSX/TDX: 0.031)
% insar{insarID}.quadtreeThresh = 0.04; % Quadtree threshold (e.g., 0.04 )
% insar{insarID}.quadtreeMindim = 32; % Quadtree threshold (e.g., 0.04 )
% insar{insarID}.quadtreeMaxdim = 256; % Quadtree threshold (e.g., 0.04 )
% insar{insarID}.gacosmaster = '/media/marin/marin/3_Working_DIR/Tirana/GACOS/A175/20191120.ztd'; % Path to GACOS ztd file on master date
% insar{insarID}.gacosslave = '/media/marin/marin/3_Working_DIR/Tirana/GACOS/A175/20191126.ztd'; % Path to GACOS ztd file on slave date
% insar{insarID}.ramp_flag = 0;
% insar{insarID}.def_bbox = [19.35 41.65;19.65 41.25];
% insar{insarID}.name = 'Tirana_t_175_20_26_asc';

% insarID = 1;                            % InSAR dataset unique identifier
% 
% insar{insarID}.dataPath = '/media/marin/marin/3_Working_DIR/Tirana/MODEL/DATA/D_T153_20191125_20191201.mat'; % Path to data file
% insar{insarID}.wavelength = 0.056;      % Wavelength in m (e.g., Envisat/ERS/Sentinel: 0.056; CSK/TSX/TDX: 0.031)
% insar{insarID}.quadtreeThresh = 0.04; % Quadtree threshold (e.g., 0.04 )
% insar{insarID}.quadtreeMindim = 32; % Quadtree threshold (e.g., 0.04 )
% insar{insarID}.quadtreeMaxdim = 256; % Quadtree threshold (e.g., 0.04 )
% insar{insarID}.gacosmaster = '/media/marin/marin/3_Working_DIR/Tirana/GACOS/D153/20191125.ztd'; % Path to GACOS ztd file on master date
% insar{insarID}.gacosslave = '/media/marin/marin/3_Working_DIR/Tirana/GACOS/D153/20191201.ztd'; % Path to GACOS ztd file on slave date
% insar{insarID}.ramp_flag = 0;
% insar{insarID}.def_bbox = [19.35 41.65;19.65 41.25];
% insar{insarID}.name = 'Tirana_t_153_dsc';

insarID = 1;                            % InSAR dataset unique identifier

insar{insarID}.dataPath = '/media/marin/marin/3_Working_DIR/Tirana/MODEL/DATA/A_T175_20191114_20191202.mat'; % Path to data file
insar{insarID}.wavelength = 0.056;      % Wavelength in m (e.g., Envisat/ERS/Sentinel: 0.056; CSK/TSX/TDX: 0.031)
insar{insarID}.quadtreeThresh = 0.04; % Quadtree threshold (e.g., 0.04 )
insar{insarID}.quadtreeMindim = 32; % Quadtree threshold (e.g., 0.04 )
insar{insarID}.quadtreeMaxdim = 256; % Quadtree threshold (e.g., 0.04 )
insar{insarID}.gacosmaster = '/media/marin/marin/3_Working_DIR/Tirana/GACOS/A175/20191114.ztd'; % Path to GACOS ztd file on master date
insar{insarID}.gacosslave = '/media/marin/marin/3_Working_DIR/Tirana/GACOS/A175/20191202.ztd'; % Path to GACOS ztd file on slave date
insar{insarID}.ramp_flag = 0;
insar{insarID}.def_bbox = [19.35 41.65;19.65 41.25];
insar{insarID}.name = 'Tirana_t_175_14_02_asc';


%% Apply GACOS tropospheric and linear ramp correction, then do quadtree
% resample and calculate variance-covariance matrix and export

    if tropocorrection_step == 1
        tropo_choice = 'yes';
    else
        tropo_choice ='no';
    end

    if covariance_step == 1
        covariance_choice = 'yes';
    else
        covariance_choice ='no';
    end
    
    
    if quadtree_step == 1
        quadtree_choice = 'yes';
    else
        quadtree_choice ='no';
    end
     

for i = 1 : length(insar)
    
    % Step 1, Load InSAR data
    
    fprintf('\n################################################\n')
    fprintf('############# LOAD INSAR DATA ##################\n')
    fprintf('################################################\n')
    
    insarData = load(insar{i}.dataPath);
    insarData.wavelength = insar{i}.wavelength;
    insarData.def_bbox = insar{i}.def_bbox;
    insarData.name = insar{i}.name;
    
    
    fprintf('InSAR data : %s\n',insar{i}.dataPath)
    fprintf('Tropo correction : %s\n',tropo_choice)
    fprintf('Quadtree downsample : %s\n',quadtree_choice)
    fprintf('Covariance calculation : %s\n',covariance_choice)
    
    if tropocorrection_step == 1
    
    % Step 2, Apply tropospheric and linear ramp correction
    
    fprintf('\n################################################\n')
    fprintf('########### APPLY TROPO CORRECTION #############\n')
    fprintf('################################################\n')
    
    fprintf('GACOS master date : %s\n',insar{i}.gacosmaster)
    fprintf('GACOS slave date  : %s\n\n',insar{i}.gacosslave)
    tic
    corr_phase{i} = gacos_tropocorr(insar{i}.gacosmaster,insar{i}.gacosslave, insarData,refPoint,insar{i}.ramp_flag);
    toc
    
    td = corr_phase{i};
    save(string(['./tropo/',insar{i}.name,'_topocor.mat']),'-struct','td');
    
    else
        corr_phase{i} = insarData;
    end
    
    % Step 3, Downsample data using Quadtree resampling algorithm
    if quadtree_step == 1
    fprintf('\n################################################\n')
    fprintf('############ QUADTREE DOWNSAMPLING #############\n')
    fprintf('################################################\n')
    
    fprintf('Quadtree threshold : %g\n',insar{i}.quadtreeThresh)
    fprintf('Quadtree min dimensions block : %g\n',insar{i}.quadtreeMindim)
    fprintf('Quadtree max dimensions block : %g\n',insar{i}.quadtreeMaxdim)
    
    fprintf('\nCrop InSAR data to bounding box\n')
    fprintf('Lon: %.2f° - %.2f°\n',boundingBox(1,1),boundingBox(3,1))
    fprintf('Lat: %.2f° - %.2f°\n',boundingBox(4,1),boundingBox(2,1))
    
    %Check if ifg corrected for tropospheric phase delay exists
    if exist(['./tropo/',insar{i}.name,'_topocor.mat'])
        corr_phase{i} = load(['./tropo/',insar{i}.name,'_topocor.mat']);
         fprintf('Loading : %s\n',['./tropo/',insar{i}.name,'_topocor.mat'])
    else
        corr_phase{i} = insarData;
    end
    
    crop_index = find(corr_phase{i}.Lon > boundingBox(1,1) & corr_phase{i}.Lon < boundingBox(3,1) & corr_phase{i}.Lat > boundingBox(4,1) & corr_phase{i}.Lat < boundingBox(2,1));
    corr_phase{i}.Lon = corr_phase{i}.Lon(crop_index,:);
    corr_phase{i}.Lat = corr_phase{i}.Lat(crop_index,:);
    corr_phase{i}.Phase = corr_phase{i}.Phase(crop_index,:);
    corr_phase{i}.Inc = corr_phase{i}.Inc(crop_index,:);
    corr_phase{i}.Heading = corr_phase{i}.Heading(crop_index,:);
     
    tic
    ds_insar{i} = qdt(corr_phase{i},insar{i}.quadtreeMindim,insar{i}.quadtreeMaxdim,insar{i}.quadtreeThresh);
    
    fprintf('\nWriting results as mat and dat file\n')
    if ~exist('./quadtree','dir')
        mkdir('quadtree')
    end
    %InSAR data
    qd = ds_insar{i};
    save(string(['./quadtree/',insar{i}.name,'.mat']),'qd');
    dlmwrite(string(['./quadtree/',insar{i}.name,'.dat']),ds_insar{i},'delimiter','\t','precision',6)
    
    fprintf('\nFinish with writing results\n')
    fprintf('Data stored in : \n%s\n%s\n', string(['./quadtree/',insar{i}.name,'.dat']),string(['./quadtree/',insar{i}.name,'.mat']))
        
    toc
    end
    % Step 4, Calculate covariance matrix, (TO DO)
    
    if covariance_step == 1
    
    fprintf('\n################################################\n')
    fprintf('############ CALCULATE COVARIANCE  #############\n')
    fprintf('################################################\n')
    
    fprintf('\nCalculating covariance matrix\n')
    
    corr_phase{i}.def_bbox = insar{i}.def_bbox;
    
    invCov{i} = calc_var(corr_phase{i},ds_insar{i}, insar{i}.ramp_flag); 
    
    % Step 5, Write the corrected-downsample Insar and covariance matrix
    % into dat file
    
    %Covariance data
    cov= invCov{i};
    if ~exist('./covariance','dir')
        mkdir('covariance')
    end
    save(string(['./covariance/',insar{i}.name,'_cov.mat']),'cov');
    csvwrite(string(['./covariance/',insar{i}.name,'_cov.csv']),invCov{i});
    
    fprintf('\nFinish with calculation of covariance matrix\n')
    
    end
       
end
    
    
    
    
    

