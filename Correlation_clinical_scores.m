%% Correlation with clinical scores 

clear
close all
clc

%% Add-ons

% The fdr_bh() function for implementing the FDR correction for multiple
% comparisons can be easily dowloaded as an add-on
% REFERENCE: David Groppe (2023). fdr_bh 
% (https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh), 
% MATLAB Central File Exchange. Retrieved June 2, 2023.

%% Loading of data

load Homework_Dataset.mat 
% The dataset was previously converted in .mat format to make the loading
% easier. Homework_Dataset.mat was made available with this code.

rng('default')
rng(1)

PD = (Demographics.Group1PD0Controls == 1); % Logical array: PD = 1, HC = 0
Radiomics_mat = table2array(Radiomics(:, 2 : end)); % The first column
                                                    % contains IDs
Radiomics_table = Radiomics(:, 2 : end);

%% Extracting selected features from the dataset
features_labels = Radiomics_table.Properties.VariableNames;
selected_features_labels = {'stat_skew', 'stat_qcod', 'ivh_v75', 'ivh_diff_v25_v75'};
% These are the selected features from the Python code for feature
% selection

% We want to extract the indexes of the selected features
selected_features_indexes = zeros(length(selected_features_labels), 1);
for i = 1 : length(selected_features_indexes)
    selected_features_indexes(i) = find(strcmp(selected_features_labels(i), features_labels));
end

selected_features = Radiomics_mat(:, selected_features_indexes);


%% Correlation between selected features and clinical scales

% N.B. This analysis is performed only on PD subjects

patients = selected_features(PD, :);
Clinical_mat = table2array(Clinical(:, 3 : end)); % Clinical features
Clinical_labels = Clinical.Properties.VariableNames(3 : end);

[rho, pVal] = corr(patients, Clinical_mat, 'type', 'Pearson');

figure()
subplot(1, 2, 1)
heatmap(Clinical_labels, selected_features_labels, rho)
title('Correlation between clinical scales and selected features')
subplot(1, 2, 2)
heatmap(Clinical_labels, selected_features_labels, pVal)
title('p-values')

% Before correction, there are some significant correlations:
% - stat_qcod        vs. UPDRSIII (r = -0.3701, p = 0.03399)
% - ivh_v75          vs. UPDRSI   (r =  0.4747, p = 0.005251)
% - ivh_v75          vs. MMSE     (r = -0.4591, p = 0.007206)
% - ivh_diff_v25_v75 vs. UPDRSIV  (r =  0.4351, p = 0.0114)

%% FDR correction

% We corrected independently for each selected feature considering the
% comparisons with all the clinical features
h_fdr = zeros(size(pVal));
crit_p = zeros(length(selected_features_labels), 1);
adj_p = zeros(size(pVal));
for ii = 1 : length(selected_features_labels)
    [h_fdr(ii, :), crit_p(ii), ~, adj_p(ii, :)] = fdr_bh(pVal(ii, :));
end

% After correction:
% - ivh_v75          vs. UPDRSI   (r =  0.4747, p_adj = 0.0288)
% - ivh_v75          vs. MMSE     (r = -0.4591, p_adj = 0.0288)

%% Plot

figure()
s = scatter(table(patients(:, 3), Clinical_mat(:, 1), 'VariableNames', {'ivh_v75', 'UPDRS I'}), 'ivh_v75', 'UPDRS I', 'filled');
s.MarkerEdgeColor = [0.7 0.0313 0.1055];
s.MarkerFaceColor = [0.7 0.0313 0.1055];
s.MarkerEdgeAlpha = 0.75;
s.MarkerFaceAlpha = 0.2;
s.LineWidth = 1;
ylim([-1 25])
r = lsline;
r.Color = [0.3125 0.3125 0.3125];
r.LineWidth = 1;
title('r = 0.47, p_{adj} = 0.03')

figure()
s = scatter(table(patients(:, 3), Clinical_mat(:, 7), 'VariableNames', {'ivh_v75', 'MMSE'}), 'ivh_v75', 'MMSE', 'filled');
s.MarkerEdgeColor = [0.7 0.0313 0.1055];
s.MarkerFaceColor = [0.7 0.0313 0.1055];
s.MarkerEdgeAlpha = 0.75;
s.MarkerFaceAlpha = 0.2;
s.LineWidth = 1;
r = lsline;
r.Color = [0.3125 0.3125 0.3125];
r.LineWidth = 1;
title('r = -0.46, p_{adj} = 0.03')