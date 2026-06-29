clc, clearvars, close all;

% Read the CSV file
% T = readtable('all_THE_fit.csv');
T = readtable('all_fits_v2.csv');

% For each unique Name, find the row with the minimum nss
[G, uniqueNames] = findgroups(T.Name);

idx_keep = splitapply(@(x) x(find(x == min(x), 1, 'first')), ...
                      (1:height(T))', G);

% Extract the selected rows
T_min_nss = T(idx_keep, :);

% Save to a new CSV file
% writetable(T_min_nss, 'less_nss_all_THE_fit.csv');
writetable(T_min_nss, 'v2_less_nss_all_THE_fit.csv');

% disp('Saved: less_nss_all_THE_fit.csv');
disp('Saved: v2_less_nss_all_THE_fit.csv');

%%% %%% %%%
% SUBSEQUENT STEP
%%% %%% %%%

clc, clearvars, close all;
%% Load data

% fitTbl = readtable('less_nss_all_THE_fit.csv');
fitTbl = readtable('v2_less_nss_all_THE_fit.csv');

S = load('participants_updrs_KANG.mat');
kangTbl = S.T;

%% Create matching IDs from fitTbl.Name

fitID = string(fitTbl.Name);

% Remove .mat
fitID = erase(fitID, ".mat");

% Remove trailing _f<number>
fitID = regexprep(fitID, '_f\d+$', '');

%% Create matching IDs from kangTbl.participant_id

kangID = string(kangTbl.participant_id);

% Remove "sub-"
kangID = erase(kangID, "sub-");

%% Prepare output columns

TD_AR_ratio = nan(height(fitTbl),1);

if iscell(kangTbl.subtype)
    subtype = cell(height(fitTbl),1);
else
    subtype = strings(height(fitTbl),1);
end

%% Match and copy values

for i = 1:height(fitTbl)

    idx = find(strcmp(fitID(i), kangID));

    if ~isempty(idx)

        idx = idx(1);   % use first match if somehow multiple exist

        TD_AR_ratio(i) = kangTbl.TD_AR_ratio(idx);

        if iscell(kangTbl.subtype)
            subtype{i} = kangTbl.subtype{idx};
        else
            subtype(i) = kangTbl.subtype(idx);
        end

    end

end

%% Append columns

fitTbl.TD_AR_ratio = TD_AR_ratio;
fitTbl.subtype     = subtype;

%% Save

writetable(fitTbl, 'v2_less_nss_all_THE_fit_with_KANG.csv');

fprintf('Saved %s\n', 'v2_less_nss_all_THE_fit_with_KANG.csv');