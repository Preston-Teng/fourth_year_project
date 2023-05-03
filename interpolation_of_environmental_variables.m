%% interpolation with environmental variables (temp, DO, MLD, NPP)

%% Temperature
% Read monthly data
% Monthly temperature data go down to depth of 1500 m (we need annual data to get values in deeper depths)
lonMonthly = double(ncread('woa18_decav_t01_01.nc','lon')); % the 12 files have the same lon x lat x depth arrangement, only need to read one month
latMonthly = double(ncread('woa18_decav_t01_01.nc','lat')); 
depthMonthly = double(ncread('woa18_decav_t01_01.nc','depth'));

nx = length(lonMonthly);
ny = length(latMonthly);
nz = length(depthMonthly);

monthlyTemp = zeros([nx ny nz 12]);
for iMonth = 1:12
   monthlyTemp(:,:,:,iMonth) = ncread(['woa18_decav_t' sprintf('%02d',iMonth) '_01.nc'], 't_an');
end

% Read annual data ('depthAnnual' has values > 1500 m, 'depthMonthly' only goes to 1500)
lonAnnual = double(ncread('woa18_decav_t00_01.nc','lon'));
latAnnual = double(ncread('woa18_decav_t00_01.nc','lat'));
depthAnnual = double(ncread('woa18_decav_t00_01.nc','depth'));

annualTemp(:,:,:) = double(ncread('woa18_decav_t00_01.nc','t_an'));

% Create an array with monthly and annual values
iLastDepthMonthlyTemp = find(depthAnnual == depthMonthly(end));
deepAnnualTemp = annualTemp(:,:,(iLastDepthMonthlyTemp+1:end));

nDepths = length(depthAnnual); 
WOA18temp = zeros([nx ny nDepths 12]);
for iMonth = 1:12
    WOA18temp(:,:,:,iMonth) = cat(3,monthlyTemp(:,:,:,iMonth),deepAnnualTemp(:,:,:));
end

% Save the temperatrue array. This will be the temperature array to use from now onwards, it has data until the depth of 5500 m
tempLat = latMonthly;
tempLon = lonMonthly;
tempDepth = cat(1,depthMonthly,depthAnnual(iLastDepthMonthlyTemp+1:end));

save(‘/Users/Preston/Desktop/Fourth_Year_Project/Analysis_Folder/world_ocean_atlas/WOA18temp.mat’,...
    'WOA18temp','tempDepth','tempLon','tempLat')

% Plot to check
figure()
pcolor(flipud(rot90(WOA18temp(:,:,1,1)))); 
caxis([-2 25]); 
cb = colorbar('FontSize', 12); 
cb.Label.String = 'Temperature (ºC)';
shading interp
colormap(jet)
box on

%% 4D integration by lat, lon and depth - getting the weighted temperature average of every unique location

load('/Users/Preston/Desktop/Fourth_Year_Project/Analysis_Folder/world_ocean_atlas/WOA18temp.mat',...
    'WOA18temp','tempDepth','tempLon','tempLat')

% number of unique locations
nd=length(data_unique.lon);

% create a new variable to store the interpolated temperature values
data_unique.Tint = zeros(nd,1); % weighted average temperature that we need
data_unique.Tmax = zeros(nd,1); % just in case
data_unique.Tmin = zeros(nd,1); % just in case

% Define xw, yw, and zw based on the WOA data
xw = tempLon;
yw = tempLat;
zw = tempDepth;

% calculate the thickness of each depth layer
dzw = diff(zw);

% loop over each data point (unique location)
for id=1:nd

    % extract location and depth range for each data point
    xd=data_unique.lon(id);
    yd=data_unique.lat(id);
    zmin=data_unique.depth_min_m(id);
    zmax=data_unique.depth_max_m(id);
    
    % find the available depths that fall within the desired range
    izm=find(zw>=zmin & zw<=zmax);
    
    % to deal with cases where zmin = zmax or if (zmin-zmax) < dzw
    
    if (isempty(izm)) 
        [minDistance1, idxOfMin1] = min(abs(zw-zmin));
        [minDistance2, idxOfMin2] = min(abs(zw-zmax));

        if (minDistance1 < minDistance2)
            izm = idxOfMin1;
        else
            izm = idxOfMin2;
        end
    end
        
    % Interpolate the temperature at the specific location
    nd1=length(izm);
    xd1=repmat(xd,[nd1 1]);
    yd1=repmat(yd,[nd1 1]);
    
    % interpolate temperature values at the available depths (for the 2
    % cases)
        
    if (length(izm) > 1) % for those cases where izm is a vector
        [Xtemp,Ytemp,Ztemp]= ndgrid(xw,yw,zw(izm));
        Tcol = squeeze(interpn(Xtemp,Ytemp,Ztemp,WOA18temp(:,:,izm),xd,yd,zw(izm)) );

    elseif (length(izm) == 1) % for those cases where zw(izm) is just one value
        [Xtemp,Ytemp]= ndgrid(xw,yw);
        Tcol = interpn(Xtemp,Ytemp,WOA18temp(:,:,izm),xd,yd);

    end
    
    % calculate the average and median temperature values
    % removes the last element of izm from the selection because it corresponds to the upper depth limit zmax and should not be included in the average calculation
    Tmed=median(Tcol,'omitnan');
    
    if (length(izm) == 1) % prevent such cases from returning NaN
        Tavg = Tmed;
        
    elseif all(isnan(Tcol)) % prevent such cases from returning zero
    Tavg = NaN;
    Tmed = NaN;
    Tmax = NaN;
    Tmin = NaN;
        
    else
        Tavg=nansum(Tcol.*[dzw(izm(1:end-1)); zmax-tempDepth(izm(end))])/sum([dzw(izm(1:end-1)); zmax-tempDepth(izm(end))]);     
       
    end

    Tmax=max(Tcol);
    Tmin=min(Tcol);
    
    % depth_layers = [dzw(izm(1:end-1)); zmax-tempDepth(izm(end))]; % used to check
    % depth_values = tempDepth(izm); % used to check

    % store the interpolated temperature value
    data_unique.Tint(id)=Tavg;
    data_unique.Tmax(id)=Tmax;
    data_unique.Tmin(id)=Tmin;
end


%% Dissolved Oxygen (DO)

% Read monthly data
% Monthly dissolved oxygen data go down to depth of 1500 m (we need annual data to get values in deeper depths)
lonMonthly = double(ncread('woa18_all_o01_01.nc','lon')); % the 12 files have the same lon x lat x depth arrangement
latMonthly = double(ncread('woa18_all_o01_01.nc','lat')); 
depthMonthly = double(ncread('woa18_all_o01_01.nc','depth'));

nx = length(lonMonthly);
ny = length(latMonthly);
nz = length(depthMonthly);

monthlyDissolvedOxygen = zeros([nx ny nz 12]);
for iMonth = 1:12
   monthlyDissolvedOxygen(:,:,:,iMonth) = ncread(['woa18_all_o' sprintf('%02d',iMonth) '_01.nc'], 'o_an');
end

% Read annual data ('depthAnnual' has values > 1500 m, 'depthMonthly' only goes to 1500)
lonAnnual = double(ncread('woa18_all_o00_01.nc','lon'));
latAnnual = double(ncread('woa18_all_o00_01.nc','lat'));
depthAnnual = double(ncread('woa18_all_o00_01.nc','depth'));

annualDissolvedOxygen(:,:,:) = double(ncread('woa18_all_o00_01.nc','o_an'));

% Create an array with monthly and annual values
iLastDepthMonthlyDissolvedOxygen = find(depthAnnual == depthMonthly(end));
deepAnnualDissolvedOxygen = annualDissolvedOxygen(:,:,(iLastDepthMonthlyDissolvedOxygen+1:end));

nDepths = length(depthAnnual);
WOA18DissolvedOxygen = zeros([nx ny nDepths 12]);
for iMonth = 1:12
    WOA18DissolvedOxygen(:,:,:,iMonth) = cat(3,monthlyDissolvedOxygen(:,:,:,iMonth),deepAnnualDissolvedOxygen(:,:,:));
end

% Save the dissolved oxygen array. This will be the DO array to use from now onwards, it has data until the depth of 5500 m
DissolvedOxygenLat = latMonthly;
DissolvedOxygenLon = lonMonthly;
DissolvedOxygenDepth = cat(1,depthMonthly,depthAnnual(iLastDepthMonthlyDissolvedOxygen+1:end));

save('/Users/Preston/Desktop/Fourth_Year_Project/Analysis_Folder/world_ocean_atlas/WOA18dissolved_oxygen.mat',...
    'WOA18DissolvedOxygen','DissolvedOxygenDepth','DissolvedOxygenLon','DissolvedOxygenLat')

% (5) Plot to check
figure()
pcolor(flipud(rot90(WOA18DissolvedOxygen(:,:,1,1)))); 
caxis([0 450]); 
cb = colorbar('FontSize', 12); 
cb.Label.String = 'Dissolved Oxygen (µmol/kg)';
shading interp
colormap(jet)
box on

%% 4D integration by lat, lon and depth - getting the weighted dissolved oxygen average of every unique location

load('/Users/Preston/Desktop/Fourth_Year_Project/Analysis_Folder/world_ocean_atlas/WOA18dissolved_oxygen.mat',...
    'WOA18DissolvedOxygen','DissolvedOxygenDepth','DissolvedOxygenLon','DissolvedOxygenLat')

nd=length(data_unique.lon);

% create a new variable to store the interpolated dissolved oxygen values
data_unique.DOint = zeros(nd,1);
data_unique.DOmax = zeros(nd,1);
data_unique.DOmin = zeros(nd,1);

% Define xw, yw, and zw based on the WOA data
xw = DissolvedOxygenLon;
yw = DissolvedOxygenLat;
zw = DissolvedOxygenDepth;

% calculate the thickness of each depth layer
dzw = diff(zw);

% loop over each data point
for id=1:nd
    % extract location and depth range for this data point
    xd=data_unique.lon(id);
    yd=data_unique.lat(id);
    zmin=data_unique.depth_min_m(id);
    zmax=data_unique.depth_max_m(id);
    
    % find the available depths that fall within the desired range
    izm=find(zw>=zmin & zw<=zmax);
    
    % to deal with cases where zmin = zmax or if (zmin-zmax) < dzw
    
    if (isempty(izm)) 
        [minDistance1, idxOfMin1] = min(abs(zw-zmin));
        [minDistance2, idxOfMin2] = min(abs(zw-zmax));

        if (minDistance1 < minDistance2)
            izm = idxOfMin1;
        else
            izm = idxOfMin2;
        end

    end
        
    % Interpolate the dissolved oxygen at a specific location
    nd1=length(izm);
    xd1=repmat(xd,[nd1 1]);
    yd1=repmat(yd,[nd1 1]);
    
    % interpolate dissolved oxygen values at the available depths (for the 2 cases)
        
    if (length(izm) > 1) % for those cases where izm is a vector
        [XDO,YDO,ZDO]= ndgrid(xw,yw,zw(izm));
        DOCol = squeeze(interpn(XDO,YDO,ZDO,WOA18DissolvedOxygen(:,:,izm),xd,yd,zw(izm)) );

    elseif (length(izm) == 1) % for those cases where zw(izm) is just one value
        [XDO,YDO]= ndgrid(xw,yw);
        DOCol = interpn(XDO,YDO,WOA18DissolvedOxygen(:,:,izm),xd,yd);

    end
    
    % calculate the average and median dissolved oxygen values
    % removes the last element of izm from the selection because it corresponds to the upper depth limit zmax and should not be included in the average calculation
    DOMed=median(DOCol,'omitnan');
    
    if (length(izm) == 1) % prevent such cases from returning NaN
        DOavg = DOMed;
        
    elseif all(isnan(DOCol)) % prevent such cases from returning zero
    DOavg = NaN;
    DOMed = NaN;
    DOMax = NaN;
    DOMin = NaN;
        
    else
        DOavg=nansum(DOCol.*[dzw(izm(1:end-1)); zmax-DissolvedOxygenDepth(izm(end))])/sum([dzw(izm(1:end-1)); zmax-DissolvedOxygenDepth(izm(end))]);     
       
    end

    DOMax=max(DOCol);
    DOMin=min(DOCol);
    
    % depth_layers = [dzw(izm(1:end-1)); zmax-DissolvedOxygenDepth(izm(end))]; % used to check
    % depth_values = DissolvedOxygenDepth(izm); % used to check

    % store the interpolated dissolved oxygen value
    data_unique.DOint(id)=DOavg;
    data_unique.DOmax(id)=DOMax;
    data_unique.DOmin(id)=DOMin;
end

%% Mixed Layer Depth (MLD)
% Read monthly data
lonMonthly = double(ncread('woa18_A5B7_M0201_01.nc','lon')); % the 12 files have the same lon x lat x depth arrangement
latMonthly = double(ncread('woa18_A5B7_M0201_01.nc','lat')); 

nx = length(lonMonthly);
ny = length(latMonthly);

WOA18mld = zeros([nx ny 12]);
for iMonth = 1:12
   WOA18mld(:,:,iMonth) = ncread(['woa18_A5B7_M02' sprintf('%02d',iMonth) '_01.nc'], 'M_an');
end

% Save the mld array
mldLat = latMonthly;
mldLon = lonMonthly;

save('/Users/Preston/Desktop/Fourth_Year_Project/Analysis_Folder/world_ocean_atlas/WOA18mld.mat',...
    'WOA18mld','mldLat','mldLon')

% Plot to check
figure()
pcolor(flipud(rot90(WOA18mld(:,:,1)))); 
caxis([-2 25]); 
cb = colorbar('FontSize', 12); 
cb.Label.String = 'Mixed Layer Depth (m)';
shading interp
colormap(jet)
box on

%% 3D integration by lat and lon - getting the mixed layer depth of every unique location

load('/Users/Preston/Desktop/Fourth_Year_Project/Analysis_Folder/world_ocean_atlas/WOA18mld.mat',...
    'WOA18mld','mldLat','mldLon')

% mld array for later interpolation 
[Xt, Yt, Tt] = ndgrid(mldLon, mldLat, (1:12));
Fmld = griddedInterpolant(Xt, Yt, Tt, WOA18mld, 'linear'); % linear interpolation

% Create interpolation array
nInterpolationPoints = length(data_unique.lat);
qMLD = NaN(nInterpolationPoints,1); % create a column of NaN values to store the interpolated temperature

% Loop over the cell array
for i = 1:nInterpolationPoints

    Lat = data_unique.lat(i);
    Lon = data_unique.lon(i);
    Month = data_unique.month(i);
    
    % Generate the query points for interpolation
    [qX, qY, qN] = ndgrid(Lon, Lat, Month);
 
    % Interpolate the mixed layer depth at the query points
    qMLD(i) = Fmld(qX, qY, qN);

end

% store the mld value
data_unique.mld = qMLD;

%% Net Primary Production

%% 3D integration by lat and lon - getting the mixed layer depth of every unique location
%% Load VGPM dataset for NPP
load('/Users/Preston/Desktop/Fourth_Year_Project/Analysis_Folder/npp_vgpm.mat')

% NPP array for later interpolation 
[Yn, Xn, Tn] = ndgrid(npp_lat, npp_lon, (1:12));
Fnpp = griddedInterpolant(Yn, Xn, Tn, npp_avg, 'linear');

% Create interpolation array
nInterpolationPoints = length(data_unique.lat);
qNPP = NaN(nInterpolationPoints,1); % create a column of NaN values to store the interpolated NPP

% Loop over the cell array
for i = 1:nInterpolationPoints

    Lat = data_unique.lat(i);
    Lon = data_unique.lon(i);
    Month = data_unique.month(i);
    
    % Generate the query points for interpolation
    [qX, qY, qN] = ndgrid(Lat, Lon, Month);
 
    % Interpolate the NPP at the query points
    qNPP(i) = Fnpp(qX, qY, qN);
end

% store the NPP value
data_unique.npp = qNPP;