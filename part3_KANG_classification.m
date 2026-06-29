% %% UPDRS-III Parkinson's Disease Subtype Classification — KANG METHOD
% % Kang et al. (2005): TD / Mixed / Akinetic-Rigid (AR)
% %
% % Tremor score = mean of rest + postural/action tremor items
% % AR score     = mean of rigidity + bradykinesia items
% % Ratio        = Tremor score / AR score
% %
% % Classification (Kang cutoffs):
% %     ratio > 1.0          → Tremor Dominant (TD)
% %     0.8 <= ratio <= 1.0  → Mixed
% %     ratio < 0.8          → Akinetic-Rigid (AR)
% %
% % NOTE 1: Kang was originally defined on the OLD UPDRS (items 20-21 for
% %         tremor; 22-27,31 for AR). The item lists below map those domains
% %         onto the equivalent MDS-UPDRS Part III items present in this file.
% % NOTE 2: MATLAB prepends 'x' to column names starting with a digit,
% %         so '3_15_a' becomes 'x3_15_a'.
% 

%% UPDRS-III Parkinson's Disease Subtype Classification — KANG METHOD
% Modified version:
% 1. Replaces NaNs in every numeric column with that column's median
% 2. Saves output to user-defined path "y"

clc;
clear;
close all;

% ─────────────────────────────────────────────
% 1. Load data
% ─────────────────────────────────────────────
data_path = "/home/amur/Documents/Nikhil/alfons_dataset_full/ds004998-download/participants_updrs_off.tsv";

if ~isfile(data_path)
    error("File not found: %s", data_path);
end

T = readtable(data_path, 'FileType', 'text', 'Delimiter', '\t');
fprintf("Loaded %d participants, %d columns.\n\n", height(T), width(T));

% ─────────────────────────────────────────────
% 2. Replace NaNs with column medians
% ─────────────────────────────────────────────
fprintf("Replacing NaNs with column medians...\n");

for c = 1:width(T)

    if isnumeric(T.(c))

        col = T.(c);

        if any(isnan(col))

            med_val = median(col, 'omitnan');

            % If entire column is NaN, leave unchanged
            if ~isnan(med_val)
                col(isnan(col)) = med_val;
                T.(c) = col;
            end

        end
    end
end

fprintf("Done.\n\n");

% ─────────────────────────────────────────────
% 3. Define item sets (Kang domains -> MDS-UPDRS items)
% ─────────────────────────────────────────────

tremor_items = { ...
    'x3_15_a','x3_15_b', ...
    'x3_16_a','x3_16_b', ...
    'x3_17_a','x3_17_b','x3_17_c','x3_17_d','x3_17_e', ...
    'x3_18'};

ar_items = { ...
    'x3_3_a','x3_3_b','x3_3_c','x3_3_d','x3_3_e', ...
    'x3_4_a','x3_4_b', ...
    'x3_5_a','x3_5_b', ...
    'x3_6_a','x3_6_b', ...
    'x3_7_a','x3_7_b', ...
    'x3_8_a','x3_8_b', ...
    'x3_9', ...
    'x3_14'};

% ─────────────────────────────────────────────
% 4. Verify items are present
% ─────────────────────────────────────────────
col_names = T.Properties.VariableNames;

missing_tremor = tremor_items(~ismember(tremor_items, col_names));
missing_ar     = ar_items(~ismember(ar_items, col_names));

if ~isempty(missing_tremor)
    fprintf("WARNING — missing tremor items: %s\n", ...
        strjoin(missing_tremor, ', '));
end

if ~isempty(missing_ar)
    fprintf("WARNING — missing AR items: %s\n", ...
        strjoin(missing_ar, ', '));
end

% ─────────────────────────────────────────────
% 5. Compute mean subscores
% ─────────────────────────────────────────────
tremor_data = T{:, tremor_items};
ar_data     = T{:, ar_items};

Tremor_score = mean(tremor_data, 2);
AR_score     = mean(ar_data, 2);

% ─────────────────────────────────────────────
% 6. Compute ratio = Tremor / AR
% ─────────────────────────────────────────────
n_patients = height(T);
ratio      = NaN(n_patients,1);

for i = 1:n_patients

    t = Tremor_score(i);
    a = AR_score(i);

    if a > 0
        ratio(i) = t/a;

    elseif t > 0 && a == 0
        ratio(i) = Inf;

    end

end

% ─────────────────────────────────────────────
% 7. Classify (Kang thresholds)
% ─────────────────────────────────────────────
subtype = cell(n_patients,1);

for i = 1:n_patients

    t = Tremor_score(i);
    a = AR_score(i);
    r = ratio(i);

    if (t == 0) && (a == 0)

        subtype{i} = 'Unclassifiable';

    elseif isinf(r)

        subtype{i} = 'TD';

    elseif isnan(r)

        subtype{i} = 'Unclassifiable';

    elseif r > 1.0

        subtype{i} = 'TD';

    elseif r < 0.8

        subtype{i} = 'AR';

    else

        subtype{i} = 'Mixed';

    end

end

% ─────────────────────────────────────────────
% 8. Add results to table
% ─────────────────────────────────────────────
T.Tremor_score = Tremor_score;
T.AR_score     = AR_score;
T.TD_AR_ratio  = ratio;
T.subtype      = categorical(subtype);

% ─────────────────────────────────────────────
% 9. Display summary
% ─────────────────────────────────────────────
fprintf("=== Subtype Distribution (Kang TD / Mixed / AR) ===\n");
fprintf("Cutoffs: ratio>1.0 = TD | 0.8-1.0 = Mixed | <0.8 = AR\n\n");

subtypes_list = {'TD','Mixed','AR','Unclassifiable'};

for i = 1:numel(subtypes_list)

    n = sum(strcmp(subtype, subtypes_list{i}));

    if n > 0
        fprintf("  %-16s %d\n", subtypes_list{i}, n);
    end

end

fprintf("\n");

disp(T(:,{'participant_id', ...
          'Tremor_score', ...
          'AR_score', ...
          'TD_AR_ratio', ...
          'subtype'}));

% ─────────────────────────────────────────────
% 10. Save output
% ─────────────────────────────────────────────

y = "/home/amur/Documents/Nikhil/Codes/full_extracted/converted/export to BIG COMP";

writetable(T, y, ...
    'FileType', 'text', ...
    'Delimiter', '\t');

fprintf("\nSaved classified table → %s\n", y);