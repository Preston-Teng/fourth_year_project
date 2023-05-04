%% biome classification (grouping the Longhurst Provinces into biomes)

%% doing the index match
% raw data set
province = data_unique.province;

% data set to index match to
biome_grouping = readtable('Longhurst_Classification.xlsx'); % this is an excel of the classification based on Longhurst (2010)

province_names = biome_grouping.ProvCode;
biome_names = biome_grouping.Biome1;

% Finding indices of provinces that match with biome_names
[~, idx] = ismember(province, province_names);

% Create a new colum that contains the corresponding biomes
biome = cell(size(province));
biome(idx > 0) = biome_names(idx(idx>0));

% Append the new column
data_unique.biome = string(biome); 

%% Spearman Rank

% Global
% Define the variables and their corresponding labels
vars = {data_unique.npp, data_unique.DOint, data_unique.mld, data_unique.Tint};
labels = {'npp', 'DO', 'mld', 'temp'};

% Define the correlation types to calculate — do Pearson in addition to Spearman just in case
corr_types = {'Pearson', 'Spearman'};

% Loop over the different variables
for i = 1:length(vars)
    var = vars{i};
    label = labels{i};

    % Extract biomass data
    biomass = data_unique.biomass_concentration_mgCm3;

    % Remove NaN values from both biomass and var
    nan_idx = isnan(var);
    biomass(nan_idx) = [];
    var(nan_idx) = [];
    
    % Display variable name
    disp(['Biomass and ', label]);
    
    % Calculate correlation coefficients
    r = corrcoef(biomass, var);
    correlation_coefficient = r(1,2);

    % Display correlation coefficient
    disp(['Correlation coefficient: ', num2str(correlation_coefficient)]);

    % Loop over the correlation types
    for j = 1:length(corr_types)
        corr_type = corr_types{j};

        % Calculate and display correlation coefficient of the specified type
        r_corr = corr(biomass, var, 'type', corr_type);
        disp([corr_type, ' correlation coefficient: ', num2str(r_corr)]);
    end
    
    % Add separator between different variables
    disp(repmat('-',1,40));
end

% For all the biomes
% Define the variables and their corresponding labels
vars = {data_unique.Tint, data_unique.DOint, data_unique.mld, data_unique.npp};
labels = {'temp', 'DO', 'mld', 'NPP'};

% Define the correlation type
corr_type = 'Spearman';

% Get unique biomes
biomes = unique(data_unique.biome);

% Loop over the different biomes
for i = 1:length(biomes)
    biome = biomes{i};
    disp(['Results for ', biome, ' biome:']);
    
    % Extract data for the current biome
    biome_idx = strcmp(data_unique.biome, biome);
    biomass = data_unique.biomass_concentration_mgCm3(biome_idx);

    % Loop over the different variables
    for j = 1:length(vars)
        var = vars{j}(biome_idx);
        label = labels{j};
        biomass = data_unique.biomass_concentration_mgCm3(biome_idx);

        % Remove NaN values from both biomass and var
        nan_idx = isnan(biomass) | isnan(var);
        biomass(nan_idx) =[];
        var(nan_idx) = [];
        
        % Calculate and display correlation coefficient
        r_corr = corr(biomass, var, 'type', corr_type);
        disp([label, ' correlation coefficient: ', num2str(r_corr)]);
    end
    
    % Add separator between different variables
    disp(repmat('-',1,40));
end

%% Machine Learning - Random Forest

% Preparing data for modelling
mesozooplankton_data = struct2table(data_unique);

% Extract variables
biomass = mesozooplankton_data.biomass_concentration_mgCm3;
DO = mesozooplankton_data.DOint;
temperature = mesozooplankton_data.Tint;
mld = mesozooplankton_data.mld;
npp = mesozooplankton_data.npp;

% Remove NaN values from all variables
nan_idx = isnan(biomass) | isnan(DO) | isnan(temperature) | isnan(mld) | isnan(npp);
mesozooplankton_data(nan_idx, :) = [];

% Extract variables after removing missing data
biomass = mesozooplankton_data.biomass_concentration_mgCm3;
DO = mesozooplankton_data.DOint;
temperature = mesozooplankton_data.Tint;
mld = mesozooplankton_data.mld;
npp = mesozooplankton_data.npp;

%Combine predictor variables into a matrix
X = [DO, temperature, mld, npp];

% Define the response variable
y = biomass;

% Set the random seed for reproducibility
rng(1);

% Split data into training and validation sets
cv = cvpartition(size(X,1),'HoldOut',0.3);
idx = cv.test;

X_train = X(~idx,:);
y_train = y(~idx,:);

X_val = X(idx,:);
y_val = y(idx,:);

%% Random Forest

% Define number of cross-validation folds
num_folds = 10;

% Create a cross-validation partition object
cv_partition = cvpartition(size(X,1), 'KFold', num_folds);

% Pre-allocate variables to store the results
R2_cv = zeros(num_folds, 1);
RMSE_cv = zeros(num_folds, 1);
variable_importance_cv = zeros(num_folds, size(X,2));

% Loop over cross-validation folds
for i = 1:num_folds
    % Split the data into training and test sets for this fold
    X_train = X(cv_partition.training(i),:);
    y_train = y(cv_partition.training(i));
    X_test = X(cv_partition.test(i),:);
    y_test = y(cv_partition.test(i));
    
    % Train Random Forest model on the training set
    mdl_rf = TreeBagger(300,X_train,y_train,'Method','regression','MinLeafSize',10,'OOBPredictorImportance','on');
    
    % Make predictions on the test set
    y_pred_rf = predict(mdl_rf,X_test);
    
    % Evaluate the model performance on the test set
    R2_cv(i) = corr(y_pred_rf,y_test)^2;
    RMSE_cv(i) = sqrt(mean((y_pred_rf-y_test).^2));
    
    % Compute variable importance for this fold
    variable_importance_cv(i,:) = mdl_rf.OOBPermutedPredictorDeltaError';

end

% Compute mean and standard deviation of performance metrics across folds
R2_mean_RF = mean(R2_cv);
R2_std_RF = std(R2_cv);

RMSE_mean_RF = mean(RMSE_cv);
RMSE_std_RF = std(RMSE_cv);

% Compute mean and standard deviation of variable importance across folds
variable_importance_mean = mean(variable_importance_cv);
variable_importance_std = std(variable_importance_cv);

%% variable importance

% Define predictor names
predictor_names = {'Temperature', 'Dissolved Oxygen', 'Mixed Layer Depth', 'NPP'};

% Compute mean and standard deviation of variable importance across folds
variable_importance_mean = mean(variable_importance_cv);
variable_importance_std = std(variable_importance_cv);

% Create a bar plot of variable importance scores
figure;
bar(variable_importance_mean);
hold on;
errorbar(variable_importance_mean, variable_importance_std, 'LineStyle', 'none', 'Color', 'black');
xticklabels(predictor_names);
xtickangle(20);
set(gca, 'FontSize', 12);
ylabel('Variable Importance Score','FontSize',18);
title('Random Forest Variable Importance Scores’,'FontSize',21);

% Train Multiple Linear Regression model
mdl_rf = TreeBagger(300,X_train,y_train,'Method','regression','MinLeafSize',10,'OOBPredictorImportance','on');

% Make predictions on validation set
y_pred_RF = predict(mdl_rf, X_val);

% Calculate R2 and RMSE
R2_RF = corr(y_pred_RF,y_val)^2;
RMSE_RF = sqrt(mean((y_pred_RF-y_val).^2));

%% Multiple Linear Regression (MLR)
% Train Multiple Linear Regression model
mdl_mlr = fitlm(X_train,y_train);

% Make predictions on validation set
y_pred_mlr = predict(mdl_mlr,X_val);

% Calculate R2 and RMSE
R2_mlr = corr(y_pred_mlr,y_val)^2;
RMSE_mlr = sqrt(mean((y_pred_mlr-y_val).^2))

%% PCA

% preparing data for PCA
mesozooplankton_data = struct2table(data_unique);

% Extract variables
biomass = mesozooplankton_data.biomass_concentration_mgCm3;
DO = mesozooplankton_data.DOint;
temperature = mesozooplankton_data.Tint;
mld = mesozooplankton_data.mld;
npp = mesozooplankton_data.npp;
biome = mesozooplankton_data.biome;

% Remove NaN values from all variables
nan_idx = isnan(biomass) | isnan(DO) | isnan(temperature) | isnan(mld) | isnan(npp);
mesozooplankton_data(nan_idx, :) = [];

% variables for PCA
biomass = mesozooplankton_data.biomass_concentration_mgCm3;
DO = mesozooplankton_data.DOint;
temperature = mesozooplankton_data.Tint;
mld = mesozooplankton_data.mld;
npp = mesozooplankton_data.npp;
biome = mesozooplankton_data.biome;

% Combine data into a matrix
X = [temperature, DO, mld, npp];

% standardize the data
X_std = zscore(X);

% perform PCA
[coeff, score, latent, ~, explained] = pca(X_std);

% plot the scree plot
figure;
plot(cumsum(explained), 'bo-')
xlabel('Number of principal components')
ylabel('Cumulative explained variance (%)')

% calculate component loadings
loadings = coeff .* sqrt(latent');

% calculate component scores for each observation —> this step is just to check if “scores” are correct
comp_scores = X_std * coeff;

% project biomass onto datapoints
figure;
hold on;
scatter(comp_scores(:,1), comp_scores(:,2), 10, biomass, 'filled');
colormap('jet');
xlabel('PC1');
ylabel('PC2');
title('Biomass (mg C m-3)');
colorbar;
set(gca,'CLim',[0.01 100]); % Set the axis limits
set(gca,'ColorScale','log'); % Set the log scale
hold off;

%% biplot
figure;
biplot(coeff(:,1:2), 'Scores',comp_scores(:,1:2),'VarLabels', {'temperature', 'DO', 'MLD', 'NPP'});
xlabel('PC1');
ylabel('PC2');
title(‘Biplot');

%% project biomes onto data points
biome_labels = unique(biome);

% define marker symbols
markers = {'o', '+', '*', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};

figure;
hold on;
for i = 1:length(biome_labels);
    group_idx = (biome == biome_labels(i));
    scatter(comp_scores(group_idx,1), comp_scores(group_idx,2), [], markers{i}, 'DisplayName', biome_labels{i}, 'SizeData', 50);
end
hold off;
xlabel('PC1');
ylabel('PC2');
title('PCA Plot with Labels');


