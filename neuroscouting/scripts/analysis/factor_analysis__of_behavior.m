%% This is a script to run a Factor Analysis on the behavioral data and determine the most appropriate model

clear all

behav_arr = load('/Users/Eichenbaum/Desktop/behav_arr_n74_ALL.txt');

nFactors = 12;

[Loadings,specVar,T,stats] = factoran(behav_arr, nFactors,'rotate','promax');
behavInFASpace = behav_arr * Loadings;

eigVals = sum(Loadings.^2,1)
eigPVar = eigVals / sum(eigVals)
eigPVar_csum = cumsum(eigVals / sum(eigVals))

%% Plotting FA loadings
figure
ax = gca;
ax.FontSize = 16; 
imagesc(Loadings)
colorbar
yticks([1:18])
yticklabels({'sRT-mean','sRT-Stdev', ...
            'rRT-Coef', 'rRT-mean', 'rRT-dPrime', 'rRT-std', 'rRT-stoprate', ...
            'IC-PSE', 'IC-dPrime', 'IC-stoprate', ...
            'pRT-Acc-1', 'pRT-Acc-2', 'pRT-Acc-3', 'pRT-Acc-4', 'pRT-RT-1', 'pRT-RT-2', 'pRT-RT-3', 'pRT-RT-4'})
communalities = sum(Loadings.^2,2)


textStrings = num2str(Loadings(:), '%0.2f');
textStrings = strtrim(cellstr(textStrings));
[x, y] = meshgrid(1:5, 1:18);
hStrings = text(x(:), y(:), textStrings(:), 'Color', 'k', 'FontSize', 14, 'HorizontalAlignment', 'center');
hStrings = text(x(6,2), y(6,2), textStrings(24), 'Color', 'w', 'FontSize', 14, 'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));