%% Group matching

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

%% Loading of data

load Homework_Dataset.mat 
% The dataset was previously converted in .mat format to make the loading
% easier. Homework_Dataset.mat was made available with this code.

rng('default')
rng(1)

PD = (Demographics.Group1PD0Controls == 1); % Logical array: PD = 1, HC = 0
Demographics_mat_PD = table2array(Demographics(PD, [3 5 : end]));
Demographics_mat_HC = table2array(Demographics(not(PD), [3 5 : end]));
labels = Demographics.Properties.VariableNames([3 5 : end]);

plot_labels = {'Average age [years]', 'Average education duration [years]', ...
    'Average height [cm]', 'Average weight [kg]', 'Average BMI [kg/cm^2]'};

%% Gaussianity check 

disp('%%%%%%%%%% GAUSSIANITY CHECK (PD) %%%%%%%%%%')
fprintf('\n')

for ii = 1 : size(Demographics_mat_PD, 2)
    [h_swtest, p_swtest] = swtest(Demographics_mat_PD(:, ii), 0.01);
    disp(char(labels(ii)))
    if h_swtest == 1
        disp('We CAN reject the null hypothesis: the distribution is significantly not gaussian')
    else
        disp("We CANNOT reject the null hypothesis: there's no significant proof of non-gaussianity")
    end
    disp (['p = ' num2str(p_swtest)])
    fprintf('\n')
end

disp('%%%%%%%%%% GAUSSIANITY CHECK (HC) %%%%%%%%%%')
fprintf('\n')

for ii = 1 : size(Demographics_mat_PD, 2)
    [h_swtest, p_swtest] = swtest(Demographics_mat_HC(:, ii), 0.05);
    disp(char(labels(ii)))
    if h_swtest == 1
        disp('We CAN reject the null hypothesis: the distribution is significantly not gaussian')
    else
        disp("We CANNOT reject the null hypothesis: there's no significant proof of non-gaussianity")
    end
    disp (['p = ' num2str(p_swtest)])
    fprintf('\n')
end

% Age, education years and height are gaussian for both the groups. Weight
% and BMI are non-gaussian in at least one of the groups.


%% t-test for comparison using gaussian variables

disp('%%%%%%%%%% T-TEST FOR COMPARISON USING GAUSSIAN VARIABLES %%%%%%%%%%')
fprintf('\n')

p_ttest = zeros(3, 1);
h_ttest = zeros(3, 1);
stats_ttest = struct();
idx_gaussian = [1 2 3]; % Age, education years and height

figure()

for i = 1 : 3
    [h_ttest(i), p_ttest(i), ~, stats_ttest(i).stats] = ttest2(Demographics_mat_PD(:, idx_gaussian(i)), ...
        Demographics_mat_HC(:, idx_gaussian(i)));
    disp(char(labels(idx_gaussian(i))))
    if h_ttest(i) == 1
        fprintf("We CAN reject the null hypothesis: there's a significant difference\n" + ...
            "between the medians of the this variable between the two groups\n")
    else
        fprintf("We CANNOT reject the null hypothesis: there's not a significant difference\n" + ...
            "between the medians of the this variable between the two groups\n")
    end
    disp (['p = ' num2str(p_ttest(i))])
    fprintf('\n')

    % Plot
    subplot(2, 3, i)
    X = categorical({'PD', 'HC'});
    X = reordercats(X, {'PD','HC'});
    Y = [mean(Demographics_mat_PD(:, idx_gaussian(i))) mean(Demographics_mat_HC(:, idx_gaussian(i)))];
    Z = [std(Demographics_mat_PD(:, idx_gaussian(i))) std(Demographics_mat_HC(:, idx_gaussian(i)))];
    b = bar(X, Y);
    b.FaceColor = 'flat';
    b.EdgeColor = 'none';
    b.CData(1, :) = [0 0 0];
    b.CData(2, :) = [0.7 0.0313 0.1055];
    b.FaceAlpha = 0.2;
    hold on
    er = errorbar(X(1), Y(1), Z(1), Z(1));
    er.Color = [0 0 0];                           
    er.LineStyle = 'none';
    er.LineWidth = 1;
    er = errorbar(X(2), Y(2), Z(2), Z(2));
    er.Color = [0.7 0.0313 0.1055];                          
    er.LineStyle = 'none';
    er.LineWidth = 1;
    b = bar(X(1), Y(1));
    set(b,'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceColor', 'none', 'EdgeAlpha', 0.75);
    b = bar(X(2), Y(2));
    set(b,'EdgeColor', [0.7 0.0313 0.1055], 'LineWidth', 1.5, 'FaceColor', 'none', 'EdgeAlpha', 0.75);
    ylabel(char(plot_labels(i)))
end


%% Wilcoxon-Mann-Whitney Univariate Test for comparison using non-gaussian variables

disp('%%%%%%%%%% WILCOXON-MANN-WHITNEY U TEST FOR COMPARISON USING NON-GAUSSIAN VARIABLES %%%%%%%%%%')
fprintf('\n')

p_ranksum = zeros(2, 1);
h_ranksum = zeros(2, 1);
stats_ranksum = struct();
idx_non_gaussian = [4 5]; % Weigt and BMI

for i = 1 : 2
    [p_ranksum(i), h_ranksum(i), stats_ranksum(i).stats] = ranksum(Demographics_mat_PD(:, idx_non_gaussian(i)), ...
        Demographics_mat_HC(:, idx_non_gaussian(i)));
    disp(char(labels(idx_non_gaussian(i))))
    if h_ttest(i) == 1
        fprintf("We CAN reject the null hypothesis: there's a significant difference\n" + ...
            "between the medians of the this variable between the two groups\n")
    else
        fprintf("We CANNOT reject the null hypothesis: there's not a significant " + ...
            "difference\nbetween the medians of the this variable between the two groups\n")
    end
    disp (['p = ' num2str(p_ranksum(i))])
    fprintf('\n')

    % Plot
    subplot(2, 3, i + 3)
    X = categorical({'PD', 'HC'});
    X = reordercats(X, {'PD','HC'});
    Y = [mean(Demographics_mat_PD(:, idx_gaussian(i))) mean(Demographics_mat_HC(:, idx_gaussian(i)))];
    Z = [std(Demographics_mat_PD(:, idx_gaussian(i))) std(Demographics_mat_HC(:, idx_gaussian(i)))];
    b = bar(X, Y);
    b.FaceColor = 'flat';
    b.EdgeColor = 'none';
    b.CData(1, :) = [0 0 0];
    b.CData(2, :) = [0.7 0.0313 0.1055];
    b.FaceAlpha = 0.2;
    hold on
    er = errorbar(X(1), Y(1), Z(1), Z(1));
    er.Color = [0 0 0];                           
    er.LineStyle = 'none';
    er.LineWidth = 1;
    er = errorbar(X(2), Y(2), Z(2), Z(2));
    er.Color = [0.7 0.0313 0.1055];                          
    er.LineStyle = 'none';
    er.LineWidth = 1;
    b = bar(X(1), Y(1));
    set(b,'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceColor', 'none', 'EdgeAlpha', 0.75);
    b = bar(X(2), Y(2));
    set(b,'EdgeColor', [0.7 0.0313 0.1055], 'LineWidth', 1.5, 'FaceColor', 'none', 'EdgeAlpha', 0.75);
    ylabel(char(plot_labels(i + 3)))
end

%% Gender

disp('%%%%%%%%%% CHI-SQUARED TEST FOR COMPARISON OF GENDER RATIOS %%%%%%%%%%')
fprintf('\n')

Gender = table2array(Demographics(:, 4));

Gender_PD = Gender(PD);
Gender_HC = Gender(not(PD));

% Observed data
n1 = sum(Gender(PD) == 1); 
N1 = sum(PD);
n2 = sum(Gender(not(PD)) == 1); 
N2 = sum(not(PD));

% Pooled estimate of proportion
p0 = (n1 + n2) / (N1 + N2);

% Expected counts under H0 (null hypothesis)
n10 = N1 * p0;
n20 = N2 * p0;

% Chi-square test, by hand
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
chi2stat = sum((observed - expected) .^ 2 ./ expected);
p_chi = 1 - chi2cdf(chi2stat, 1);

if p_chi <= 0.05
    fprintf("We CAN reject the null hypothesis: there's a significant difference\n" + ...
        "between the gender frequencies between the two groups\n\n")
else
    fprintf("We CANNOT reject the null hypothesis: there's not a significant " + ...
        "difference\nbetween the gender frequencies between the two groups\n\n")
end

% Plot
subplot(2, 3, 6)
X = categorical({'PD', 'HC'});
X = reordercats(X, {'PD','HC'});
Y = [sum(Gender(PD) == 1) sum(Gender(PD) == 0);  sum(Gender(not(PD)) == 1) sum(Gender(not(PD)) == 0)];
b = bar(X, Y, 'stacked', 'FaceColor', 'flat');
b(1, 1).CData = [0 0 0];
b(1, 2).CData = [0.7 0.0313 0.1055];
b(1, 1).FaceAlpha = 0.2;
b(1, 2).FaceAlpha = 0.2;
set(b(1, 1), 'EdgeColor', 'k', 'LineWidth', 1.5, 'EdgeAlpha', 0.75);
set(b(1, 2), 'EdgeColor', [0.7 0.0313 0.1055], 'LineWidth', 1.5, 'EdgeAlpha', 0.75);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1, [9, 2], labels1, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize', 9, 'FontWeight', 'bold')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2, [30, 17], labels2, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize', 9, 'Color', [0.7 0.0313 0.1055], 'FontWeight', 'bold')
legend('Males', 'Females')
ylabel('Number of subjects')