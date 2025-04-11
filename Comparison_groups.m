%% Comparison between groups

clear
close all
clc

%% Add-ons

% The swtest() function for implementing the Shapiro-Wilk gaussianity test
% is part of the 'Shapiro-Wilk and Shapiro-Francia normality tests' package
% that can be easility imported in MATLAB as an add-on
% REFERENCE: Ahmed BenSaïda (2023). Shapiro-Wilk and Shapiro-Francia 
% normality tests. 
% (https://www.mathworks.com/matlabcentral/fileexchange/13964-shapiro-wilk-and-shapiro-francia-normality-tests), 
% MATLAB Central File Exchange. Retrieved June 1, 2023.


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

%% Check for gaussianity of selected features 

disp('%%%%%%%%%% GAUSSIANITY CHECK %%%%%%%%%%')
fprintf('\n')

for ii = 1 : size(selected_features, 2)
    [h_swtest, p_swtest] = swtest(selected_features(:, ii), 0.05);
    disp(char(selected_features_labels(ii)))
    if h_swtest == 1
        disp('We CAN reject the null hypothesis: the distribution is significantly not gaussian')
    else
        disp("We CANNOT reject the null hypothesis: there's no significant proof of non-gaussianity")
    end
    disp (['p = ' num2str(p_swtest)])
    fprintf('\n')
end
fprintf('\n\n\n')

% stat_qcod is gaussian and the others are not

%% Separating patients from controls

patients = selected_features(PD, :);
controls = selected_features(not(PD), :);

%% Correlation between demographics and selected features

Demographics_mat = table2array(Demographics(:, 3 : end)); % Demographic features
Demographics_labels = Demographics.Properties.VariableNames(3 : end);

[rho, pVal] = corr(selected_features, Demographics_mat, 'type', 'Spearman');

figure()
subplot(1, 2, 1)
h = heatmap(Demographics_labels, selected_features_labels, rho);
title('Correlation between demographics and selected features')
h.CellLabelFormat = '%.3f';
subplot(1, 2, 2)
h = heatmap(Demographics_labels, selected_features_labels, pVal);
title('p-values')
h.CellLabelFormat = '%.3f';

% p-values are extremely high, there's not even need to perform multiple 
% comparison correction
% No significant correlations -> We won't include the demographics as
%                                covariates

%% t-test for comparison using gaussian variables

disp('%%%%%%%%%% T-TEST FOR COMPARISON USING GAUSSIAN VARIABLES %%%%%%%%%%')
fprintf('\n')

disp(char(selected_features_labels(2))) % The second variable is the only 
                                        % gaussian

[h_ttest, p_ttest, ~, stats_ttest] = ttest2(patients(:, 2), controls(:, 2));
if h_ttest == 1
    fprintf("We CAN reject the null hypothesis: there's a significant difference\n" + ...
        "between the means of the this variable between the two groups\n")
else
    fprintf("We CANNOT reject the null hypothesis: there's not a significant difference\n" + ...
        "between the means of the this variable between the two groups\n")
end
disp (['p = ' num2str(p_ttest)])
fprintf('\n\n\n')

% The following steps are needed to plot the boxplots in the same space
labels = {'PD', 'HC'};
groups = categorical([repmat({'PD'}, size(patients, 1), 1); repmat({'HC'}, size(controls, 1), 1)]);

% Plot
figure()
subplot(2, 2, 1)
b = boxchart([patients(:, 2); controls(:, 2)], 'GroupByColor', groups);
a = gca;
b(1).JitterOutliers = 'on';
b(1).JitterOutliers = 'on';
b(1).BoxFaceColor = 'k';
b(2).BoxFaceColor = [0.7 0.0313 0.1055];
b(1).MarkerColor = 'k';
b(2).MarkerColor = [0.7 0.0313 0.1055];
a.XTick = [];
legend(labels)
a.YLabel.Interpreter = 'none';
ylabel(selected_features_labels(2))
xtext = 1;
ytext = max([controls(:, 2); patients(:, 2)] + 0.05);
text(a, xtext, ytext, '****', 'HorizontalAlignment', 'center')
ylim([min([controls(:, 2); patients(:, 2)]) - 0.05, max([controls(:, 2); patients(:, 2)]) + 0.05])

%% Wilcoxon-Mann-Whitney Univariate Test for comparison using non-gaussian variables

disp('%%%%%%%%%% WILCOXON-MANN-WHITNEY U TEST FOR COMPARISON USING NON-GAUSSIAN VARIABLES %%%%%%%%%%')
fprintf('\n')

p_ranksum = zeros(3, 1);
h_ranksum = zeros(3, 1);
stats_ranksum = struct();
idx_non_gaussian = [1 3 4]; % stat_qcod, which is the second variable in
                            % the selected_features matrix, is gaussian, so
                            % it doesn't have to be accounted for in this
                            % section

for i = 1 : 3
    [p_ranksum(i), h_ranksum(i), stats_ranksum(i).stats] = ranksum(controls(:, idx_non_gaussian(i)), patients(:, idx_non_gaussian(i)));
    %stats_ranksum = [stats_ranksum stats_ranksum_tmp];
    disp(char(selected_features_labels(idx_non_gaussian(i))))
    if h_ranksum(i) == 1
        fprintf("We CAN reject the null hypothesis: there's a significant difference\nbetween the medians of the this variable between the two groups\n")
    else
        fprintf("We CANNOT reject the null hypothesis: there's not a significant difference\nbetween the medians of the this variable between the two groups\n")
    end
    disp (['p = ' num2str(p_ranksum(i))])
    fprintf('\n')

    % Plot
    subplot(2, 2, i + 1)
    b = boxchart([patients(:, idx_non_gaussian(i)); controls(:, idx_non_gaussian(i))], 'GroupByColor', groups);
    a = gca;
    b(1).JitterOutliers = 'on';
    b(1).JitterOutliers = 'on';
    b(1).BoxFaceColor = 'k';
    b(2).BoxFaceColor = [0.7 0.0313 0.1055];
    b(1).MarkerColor = 'k';
    b(2).MarkerColor = [0.7 0.0313 0.1055];
    a.XTick = [];
    a.YLabel.Interpreter = 'none';
    legend(labels)
    ylabel(selected_features_labels(idx_non_gaussian(i)))
    xtext = 1;
    ytext = max([controls(:, idx_non_gaussian(i)); patients(:, idx_non_gaussian(i))]) + 0.05;
    text(a, xtext, ytext, '****', 'HorizontalAlignment', 'center')
    ylim([min([controls(:, idx_non_gaussian(i)); patients(:, idx_non_gaussian(i))]) - 0.05, max([controls(:, idx_non_gaussian(i)); patients(:, idx_non_gaussian(i))]) + 0.05])
end
