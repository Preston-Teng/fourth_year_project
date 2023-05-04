%% biomass unit conversion

% zooplankton dry mass
biomass_dry = readtable('zooplankton_biomass_dry_combined.xlsx');

% zooplankton wet mass
biomass_wet = readtable('zooplankton_biomass_wet_combined.xlsx');

% zooplankton displacement volume
biovolume = readtable('zooplankton_biovolume_combined.xlsx');

% zooplankton C mass
biomass_carbon = readtable('zooplankton_biomass_c_combined.xlsx');

%% Wiebe (1988) conversion units

% conversion of dry mass to carbon mass
biomass_dry_array = table2array(biomass_dry(:,12));
biomass_c_from_dry = 10.^((log10(biomass_dry_array) - 0.499)/0.991); % using the conversion factor from Wiebe (1988)
biomass_c_from_dry = array2table(biomass_c_from_dry); % converting back to table
biomass_dry(:,12) = biomass_c_from_dry; % appending data back to dataset

% conversion of wet mass to carbon mass
biomass_wet_array = table2array(biomass_wet(:,12));
biomass_wet_array = biomass_wet_array/1000; % wet biomass values to be converted from mg to g for Wiebe's formula
biomass_c_from_wet = 10.^((log10(biomass_wet_array) + 1.537)/0.852); % using the conversion factor from Wiebe (1988)
biomass_c_from_wet = array2table(biomass_c_from_wet); % converting back to table
biomass_wet(:,12) = biomass_c_from_wet; % appending data back to dataset

%% total mass data

total_biomass = vertcat(biomass_dry, biomass_wet, biomass_carbon); % appending both wet, dry and carbon biomass into one table (after all converted to carbon biomass)
total_biomass_array = table2array(total_biomass(:,12));

% conversion of displacement volume to carbon mass
biovolume_array = table2array(biovolume(:,12));
biomass_c_from_bv = 10.^((log10(biovolume_array) + 1.434)/0.820); % using the conversion factor from Wiebe (1988)
biomass_c_from_bv = array2table(biomass_c_from_bv); % converting back to table
biovolume(:,12) = biomass_c_from_bv; % appending data back to dataset

%% total zooplankton data (after all converted to carbon mass)

biovolume = renamevars(biovolume, {'biovolume_mlm2'}, {'biomass_integrated_mgm2'}); %% renaming headers
biovolume = renamevars(biovolume, {'biovolume_concentration_mlm3'}, {'biomass_concentration_mgm3'});
zooplankton_data = vertcat(total_biomass, biovolume);

%% unique count of reported_units

% Get unique unit values and count
[unique_units, ~, unit_idx] = unique(zooplankton_data.reported_units);
unit_count = accumarray(unit_idx, 1);

% Calculate percentage
unit_percent = (unit_count / sum(unit_count)) * 100;

% Display results
fprintf('Unit\tCount\tPercent\n');
for ii = 1:numel(unique_units)
    fprintf('%s\t%d\t%.2f%%\n', unique_units{ii}, unit_count(ii), unit_percent(ii));
end

%% unique count of units_before_conversion

% Get unique unit values and count
[unique_units, ~, unit_idx] = unique(zooplankton_data.unit_before_conversion);
unit_count = accumarray(unit_idx, 1);

% Calculate percentage
unit_percent = (unit_count / sum(unit_count)) * 100;

% Display results
fprintf('Unit\tCount\tPercent\n');

for ii = 1:numel(unique_units)
    fprintf('%s\t%d\t%.2f%%\n', unique_units{ii}, unit_count(ii), unit_percent(ii));
end

%% mesh size standardisation

% remove the two studies which used a CTD instead of zooplankton mesh
zooplankton_data = zooplankton_data(~isnan(zooplankton_data.mesh_size), :);

% extract mesh sizes from zooplankton_data into a variable
mesh_sizes = zooplankton_data.mesh_size;

% calculate the bin edges for the histogram
min_mesh = min(mesh_sizes);
max_mesh = max(mesh_sizes);
bin_width = 1;
bin_edges = min_mesh:bin_width:max_mesh;

% create the histogram
figure;
histogram(mesh_sizes, bin_edges);

% create the bar chart
[counts, edges] = histcounts(mesh_sizes, bin_edges);
bar(bin_edges(1:end-1), counts, bin_width, 'FaceColor', 'b', 'EdgeColor', 'k');

% Set the x-axis limit
xlim([100 650]); % mesozooplankton size range we set

% labels
xlabel('Mesh size');
ylabel('Count');
title('Distribution of mesozooplankton mesh sizes');

% create the table of mesh sizes and their frequency
t = tabulate(mesh_sizes);
t = array2table(t, 'VariableNames', {'Mesh_Size', 'Count', 'Percentage'});
t = sortrows(t, 'Mesh_Size', 'ascend');
disp(t);

%% removing mesh size <100 and >650 for mesozooplankton size fraction

% Create a logical array to identify rows with mesh size less than 100 or greater than 650
idx = (zooplankton_data.mesh_size < 100) | (zooplankton_data.mesh_size > 650);

% Use the logical array to create a new table with the filtered rows
mesozooplankton_data = zooplankton_data(~idx, :);

%% mesh unit standardisation

% Remove rows with NaN values in the mesh_size column
mesozooplankton_data = mesozooplankton_data(~isnan(mesozooplankton_data.mesh_size), :);

% Define the bin centers for the mesh sizes
bin_centers = [125, 150, 200, 280, 330, 500];

% Use the discretize function to group the mesh sizes into bins based on which bin center they are closest to
[~, group_idx] = min(abs(mesozooplankton_data.mesh_size - bin_centers), [], 2);

% Replace the group index values with the corresponding bin centers
group_centers = bin_centers(group_idx);

% Add the new column to the table
mesozooplankton_data.group = group_centers';

% to check if code works correctly
% [za, ~, zb] = unique([mesozooplankton_data.mesh_size mesozooplankton_data.group], 'rows');

%% mesh conversion using the equations provided by Yang (2022)

% Define the factors for each group
group_factors = [1.25623, 1.14296, 1, 0.83435, 0.75899, 0.51169];

% Use the discretize function to group the mesh sizes into bins based on which bin center they are closest to
[~, group_idx] = min(abs(mesozooplankton_data.mesh_size - bin_centers), [], 2);

% Get the corresponding factor for each row in the table
group_factors_adj = group_factors(group_idx);

% Add the adjusted biomass concentration column to the table
mesozooplankton_data.biomass_concentration_adj_mgm3 = mesozooplankton_data.biomass_concentration_mgm3 ./ group_factors_adj';

%% ETOPO
% set up ETOPO grid
bathy = ncread('ETOPO_2022_v1_60s_N90W180_bed.nc','z');
lat = ncread('ETOPO_2022_v1_60s_N90W180_bed.nc','lat');
lon = ncread('ETOPO_2022_v1_60s_N90W180_bed.nc','lon');

save('/Users/Preston/Desktop/Fourth_Year_Project/Analysis_Folder/datasets_matlab/ETOPO.mat','bathy','lat','lon','-v7.3');

%%
% load ETOPO dataset (after grid was set up and saved previously)
load('/Users/Preston/Desktop/Fourth_Year_Project/Analysis_Folder/datasets_matlab/ETOPO.mat')

% set up bathymetry grid for interpolation
[Xbathy, Ybathy] = ndgrid(lon, lat);
Fbathy = griddedInterpolant(Xbathy, Ybathy, bathy);

% interpolation
v = mesozooplankton_data;
nLocs = size(v,1);
nLocs_new_filtered = v;

iLocToDelete = zeros(nLocs,1);
i = 1;

for iLoc = 1:nLocs

    qLat = v{iLoc,6};   
    qLon = v{iLoc,7};  

    [qX, qY] = ndgrid(qLon, qLat); % query points for interpolation

    qZooDepth = Fbathy(qX, qY);

    if (abs(qZooDepth) <= 200)
        iLocToDelete(i) = iLoc;
        i = i + 1;
    end
end

nLocsToDelete = sum(iLocToDelete > 0);

% Delete rows
nLocs_new_filtered(iLocToDelete(1:nLocsToDelete),:) = []; 

% How many locations in the zooarray did we delete?
fracLocsToDelete = 100*(nLocsToDelete/nLocs); 

% saving the filtered ETOPO data
mesozooplankton_data_filtered = nLocs_new_filtered;

%% Assigning Longhurst provinces

fullpathConfigDir    = '/Users/Preston/Desktop/Fourth_Year_Project/config/';
addpath(genpath(fullpathConfigDir));
addpath(strcat(fullpathConfigDir,'longhurst_v4_2010'));

% Assign Longhurst provinces to the mesozooplankton sampling locations
% Note: the longitudes used to define the polygons are signed. Do NOT change them as it will mess up the polygon          % definitions! Instead, we fix the data location
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
fixLong = @(l) iif( l(1)>180, @() [l(1)-360 l(2)], ...
                   true, @() l, ...
                   true, @() l );

obsloc=mat2cell([lon lat],ones(1,length(lon)),[2]);
obsloc=cellfun(fixLong,obsloc,'UniformOutput',0);

% This assigns a province index to each data point
provinces=cellfun(@(l) find(cellfun(@(p)inpolygon(l(1),l(2),p(:,1),p(:,2)),provs.ncst)),obsloc,'UniformOutput',0);

% Some points are not assigned an index; they will be assigned an [] value below
province_name=cellfun(@(p)provs.dbf.ProvCode(p),provinces,'UniformOutput',0);
province_descr=cellfun(@(p)provs.dbf.ProvDescr(p),provinces,'UniformOutput',0);

% Adding new columns to the zooplankton array, which will contain the Longhurst province information
province_name = cell(size(provinces));
province_descr = cell(size(provinces));

% this loop is to remove points exactly on the boundaries of the Longhurst provinces
for i = 1:length(provinces)
    if length(provinces{i}) > 1
        province_name{i} = provs.dbf.ProvCode(provinces{i}(1));
        province_descr{i} = provs.dbf.ProvDescr(provinces{i}(1));
    else
        province_name{i} = provs.dbf.ProvCode(provinces{i});
        province_descr{i} = provs.dbf.ProvDescr(provinces{i});
    end
end

% append new columns
mesozooplankton_data_filtered.province_Longhurst2010 = province_name;
mesozooplankton_data_filtered.provincedescr_Longhurst2010 = province_descr;

% replacing the empty cells (from the Longhurst data) with NaN
locEmpties_prov = cellfun('isempty',mesozooplankton_data_filtered.province_Longhurst2010);
mesozooplankton_data_filtered.province_Longhurst2010(locEmpties_prov) = {'NaN'};

locEmpties_desc = cellfun('isempty',mesozooplankton_data_filtered.provincedescr_Longhurst2010);
mesozooplankton_data_filtered.provincedescr_Longhurst2010(locEmpties_desc) = {'NaN'};

%% removing all the points that were assigned 'NaN' as a province
% Create a logical index for rows with NaN in the province_Longhurst2010 column
nan_idx_longhurst = strcmp(mesozooplankton_data_filtered.province_Longhurst2010, ‘NaN');

% Remove rows with NaN from the table
mesozooplankton_data_filtered = mesozooplankton_data_filtered(~nan_idx_longhurst, :);

% remove MEDI and REDS
remove_provinces = cellstr({'REDS', 'MEDI'});
mesozooplankton_data_filtered.province_Longhurst2010 = cellstr(string(mesozooplankton_data_filtered.province_Longhurst2010));

% logical index to remove REDS and MEDI
logical_index_longhurst_province = ~ismember(mesozooplankton_data_filtered.province_Longhurst2010, remove_provinces); 
mesozooplankton_data_filtered = mesozooplankton_data_filtered(logical_index_longhurst_province, :);

%% Quality Control Flagging System

% Define the mesozooplankton biomass data, province data, and month data
biomass = mesozooplankton_data_filtered.biomass_concentration_adj_mgm3;
province = mesozooplankton_data_filtered.province_Longhurst2010;
month = mesozooplankton_data_filtered.month;

% Define seasons (F3 flag)
season = zeros(size(month));
season(month==12 | month<=2) = 1; 
season(month>=3 & month<=5) = 2; 
season(month>=6 & month<=8) = 3; 
season(month>=9 & month<=11) = 4; 

% Add to original table
mesozooplankton_data_filtered.season = season;

% Initialize flag arrays
F1_flag = zeros(size(biomass));
F2_flag = zeros(size(biomass));
F3_flag = zeros(size(biomass));

% Loop over all data points for F1 and F2 flag
for ii = 1:numel(biomass)
    % Calculate F1 flag
    F1_flag(ii) = sum(biomass(ii) > biomass) / numel(biomass);

    % Calculate F2 flag for the same province
    same_province = strcmp(province, province(ii)) & ~isnan(biomass);
    F2_flag(ii) = sum(biomass(ii) > biomass(same_province)) / sum(same_province);
    
    % Calculate F3 flag for the same province and season
    same_season_province = strcmp(province, province(ii)) & season == season(ii) & ~isnan(biomass);
    F3_flag(ii) = sum(biomass(ii) > biomass(same_season_province)) / sum(same_season_province);
end

% copepod 2007 threshold
QC_F1_idx = F1_flag < 0.01; 
QC_F2_idx = F2_flag < 0.001;
QC_F3_idx = F3_flag < 0.0001;

num_flaggedF1 = sum(QC_F1_idx == 1);
num_flaggedF2 = sum(QC_F2_idx == 1);
num_flaggedF3 = sum(QC_F3_idx == 1);

QC_tableF1 = mesozooplankton_data_filtered(QC_F1_idx,:);
QC_tableF2 = mesozooplankton_data_filtered(QC_F2_idx,:);
QC_tableF3 = mesozooplankton_data_filtered(QC_F3_idx,:);

% number of data points flagged multiple times
% Find data points flagged by both F1 and F2 and F3
QC_F1_F2_F3_idx = QC_F1_idx & QC_F2_idx & QC_F3_idx;
num_flaggedF1_F2_F3 = sum(QC_F1_F2_F3_idx);

% Find data points flagged by both F1 and F2
QC_F1_F2_idx = QC_F1_idx & QC_F2_idx;
num_flaggedF1_F2 = sum(QC_F1_F2_idx);
num_flaggedF1_F2_only = num_flaggedF1_F2 - num_flaggedF1_F2_F3;

% Find data points flagged by both F1 and F3
QC_F1_F3_idx = QC_F1_idx & QC_F3_idx;
num_flaggedF1_F3 = sum(QC_F1_F3_idx);
num_flaggedF1_F3_only = num_flaggedF1_F3 - num_flaggedF1_F2_F3;

% Find data points flagged by both F2 and F3
QC_F2_F3_idx = QC_F2_idx & QC_F3_idx;
num_flaggedF2_F3 = sum(QC_F2_F3_idx);
num_flaggedF2_F3_only = num_flaggedF2_F3 - num_flaggedF1_F2_F3;

% To display results
disp("Number of data points flagged by both F1 and F2: " + num_flaggedF1_F2_only)
disp("Number of data points flagged by both F1 and F3: " + num_flaggedF1_F3_only)
disp("Number of data points flagged by both F2 and F3: " + num_flaggedF2_F3_only)
disp("Number of data points flagged by all three flags: " + num_flaggedF1_F2_F3)

% table for all the flagged data points (for future reference)
% Add flag type column to each QC table
flagF1 = cell(size(QC_tableF1, 1), 1);
flagF2 = cell(size(QC_tableF2, 1), 1);
flagF3 = cell(size(QC_tableF3, 1), 1);

flagF1(:) = {'F1'};
flagF2(:) = {'F2'};
flagF3(:) = {'F3'};

QC_tableF1.flag = flagF1;
QC_tableF2.flag = flagF2;
QC_tableF3.flag = flagF3;

% Concatenate QC tables into one
QC_table = vertcat(QC_tableF1, QC_tableF2, QC_tableF3);

%% Removing the flagged points from dataset

% Find unique rows
% Replace empty time cells with "NA"
QC_table(:, 3) = fillmissing(QC_table(:, 3), 'constant', 'NA');
mesozooplankton_data_filtered(:,3) = fillmissing(mesozooplankton_data_filtered(:, 3), 'constant', 'NA');  

% Fill up NaN rows in integrated biomass
mesozooplankton_data_filtered.biomass_integrated_mgm2 = mesozooplankton_data_filtered.biomass_concentration_mgm3 .* mesozooplankton_data_filtered.average_depth_m .* 2;
QC_table.biomass_integrated_mgm2 = QC_table.biomass_concentration_mgm3 .* QC_table.average_depth_m .* 2; 

% Remove flag column so that we can get rid of duplicate data points that are flagged separately
QC_table = removevars(QC_table, 'flag'); 

% Find and remove unique rows
[unique_rows, ~, idx_duplicates] = unique(QC_table, 'rows', 'stable');

QC_table_rows = size(QC_table, 1); % number of rows before removing duplicates % QC_table has 1617 rows
QC_table_unique_rows = size(unique_rows, 1); % number of unique rows % QC_table_unique has 1284 rows

% Find the rows to remove from mesozooplankton_data_filtered
rows_to_remove = ismember(mesozooplankton_data_filtered(:,1:20), QC_table(:,1:20), 'rows');

% Remove the rows from mesozooplankton_data_filtered
mesozooplankton_data_filtered(rows_to_remove, :) = [];

%% average biomass by unique locations to reduce bias
% struct creation for data that we need for analysis (before getting unique combinations)
data.paper = string(mesozooplankton_data_filtered.paper)
data.month = mesozooplankton_data_filtered.month
data.lat = mesozooplankton_data_filtered.lat
data.lon = mesozooplankton_data_filtered.long
data.depth_min = mesozooplankton_data_filtered.top_water_depth_m
data.depth_max = mesozooplankton_data_filtered.bottom_water_depth_m
data.depth = (data.depth_max - data.depth_min)/2
data.province = string(mesozooplankton_data_filtered.province_Longhurst2010)
data.description = string(mesozooplankton_data_filtered.provincedescr_Longhurst2010)
data.season = mesozooplankton_data_filtered.season
data.biomass_carbon_concentration_mgCm3 = mesozooplankton_data_filtered.biomass_concentration_adj_mgm3

% Get unique combinations of 'lat', 'long', 'top_water_depth_m', 'bottom_water_depth_m', and 'month' 
% note: the remaining variables doesn't affect the data extraction
[unique_locs, ~, loc_indices] = unique([data.lat data.lon data.depth_min data.depth_max data.month data.paper data.depth data.province data.description data.season], 'rows');

% Compute average biomass carbon for each unique combination
biomass_carbon_concentration = accumarray(loc_indices, data.biomass_carbon_concentration_mgCm3, [], @mean);

% extract data that we need and put into a struct
% struct for data_unique (each unique location) - this is the struct to be used for future analysis 

% Add locs to unique_data
data_unique.lat = double(unique_locs(:,1));
data_unique.lon = double(unique_locs(:,2));
data_unique.depth_min_m = double(unique_locs(:,3));
data_unique.depth_max_m = double(unique_locs(:,4));
data_unique.month = double(unique_locs(:,5));
data_unique.paper = unique_locs(:,6);
data_unique.depth_m = double(unique_locs(:,7));
data_unique.province = unique_locs(:,8);
data_unique.description = unique_locs(:,9);
data_unique.season = double(unique_locs(:,10));

% adding to struct
data_unique.biomass_concentration_mgCm3 = biomass_carbon_concentration;
