clc, clearvars, close all;

%% =========================================================================
%  Build analysis table
%
%  One row per STN side (up to 2 per file).
%  Columns: Name | High_beta_electrode | full_time_series | side | best_epoch
%% =========================================================================

% ----------------------------- CONFIG ------------------------------------
folderpath  = '/home/amur/Documents/Nikhil/Codes/full_extracted/converted';
BETA_BAND   = [12 33];   % Hz
EPOCH_SEC   = 30;        % sliding window length (s)
EPOCH_STEP  = 0.5;       % step size (s)
% -------------------------------------------------------------------------

files  = dir(fullfile(folderpath, '*.mat'));
nFiles = numel(files);

maxN                = 2 * nFiles;
Name                = strings(maxN, 1);
High_beta_electrode = strings(maxN, 1);
full_time_series    = cell(maxN, 1);   % struct: .signal  .time
side                = strings(maxN, 1);
best_epoch          = cell(maxN, 1);   % struct: .signal  .time  .beta_power
row                 = 0;

%% ===================== loop over .mat files =============================
for iFile = 1:nFiles

    fname = files(iFile).name;
    fprintf('Processing %d/%d : %s\n', iFile, nFiles, fname);

    loaded = load(fullfile(folderpath, fname));
    D      = loaded.data;
    Time   = D.Time(:);
    Fs     = 1 / mean(diff(Time));

    % --- identify channel fields -----------------------------------------
    chans = fieldnames(D);
    chans = chans(~strcmp(chans, 'Time'));

    elec         = cell2mat(cellfun(@parseElec, chans, 'uni', 0));
    montageIsTwo = any(elec(:) >= 33);
    if montageIsTwo
        rightSet = 1:8;    leftSet = 33:40;
    else
        rightSet = 1:4;    leftSet = 5:8;
    end

    % --- split channels into left / right structs ------------------------
    D_left  = struct();
    D_right = struct();
    for i = 1:numel(chans)
        ab = parseElec(chans{i});
        if     all(ismember(ab, leftSet)),  D_left.(chans{i})  = D.(chans{i});
        elseif all(ismember(ab, rightSet)), D_right.(chans{i}) = D.(chans{i});
        end
    end

    % --- for each side, pick the highest-beta-power channel --------------
    for s = ["left" "right"]
        if s == "left", Dside = D_left; else, Dside = D_right; end

        f = fieldnames(Dside);
        if isempty(f), continue; end

        bp = nan(numel(f), 1);
        for c = 1:numel(f)
            bp(c) = bandPower_(Dside.(f{c}), Fs, BETA_BAND);
        end
        [~, best] = max(bp);
        bestField = f{best};

        row = row + 1;
        Name(row)                = string(fname);
        High_beta_electrode(row) = string(bestField);
        ts = struct('signal', Dside.(bestField)(:), 'time', Time);
        full_time_series{row}    = ts;
        side(row)                = s;

        % --- find best 30-sec epoch for this channel ---------------------
        best_epoch{row} = findBestEpoch(ts, Fs, BETA_BAND, EPOCH_SEC, EPOCH_STEP);
    end
end

% --- trim pre-allocated rows and assemble table --------------------------
keep                = 1:row;
Name                = Name(keep);
High_beta_electrode = High_beta_electrode(keep);
full_time_series    = full_time_series(keep);
side                = side(keep);
best_epoch          = best_epoch(keep);

Results = table(Name, High_beta_electrode, full_time_series, side, best_epoch);

fprintf('\nBuilt table with %d rows.\n', height(Results));
disp(Results(:, {'Name', 'High_beta_electrode', 'side'}));


%% ========================= HELPERS ======================================

function ep = findBestEpoch(ts, Fs, betaBand, epochSec, stepSec)
% Slide a window of epochSec seconds (step stepSec s) over ts.signal,
% compute beta power in each window, return the window with the max power.
%
%   ts        : struct with .signal (column) and .time (column)
%   Fs        : sampling rate (Hz)
%   betaBand  : [lo hi] Hz
%   epochSec  : window length in seconds   (30)
%   stepSec   : step size in seconds       (0.5)
%
%   ep        : struct with fields
%                 .signal      — beta-bandpassed epoch snippet
%                 .time        — corresponding time vector
%                 .beta_power  — integrated beta power of this epoch

    signal    = ts.signal(:);
    time      = ts.time(:);
    winLen    = round(epochSec * Fs);   % samples in one window
    stepLen   = round(stepSec  * Fs);   % samples per step
    nSamples  = numel(signal);

    % --- bandpass full signal once before windowing ----------------------
    %     (filter on the full length for clean edge behaviour)
    sig_bp = matlab_bpass_filt(signal, betaBand(1), betaBand(2), Fs);

    % start indices of every valid window
    starts = 1 : stepLen : (nSamples - winLen + 1);

    if isempty(starts)
        warning('Recording shorter than epoch length — returning full signal as epoch.');
        ep = struct('signal', sig_bp, 'time', time, 'beta_power', bandPower_(sig_bp, Fs, betaBand));
        return
    end

    nWins   = numel(starts);
    bpWin   = nan(nWins, 1);

    for k = 1:nWins
        idx      = starts(k) : starts(k) + winLen - 1;
        bpWin(k) = bandPower_(sig_bp(idx), Fs, betaBand);
    end

    [maxPow, bestK] = max(bpWin);
    idx_best = starts(bestK) : starts(bestK) + winLen - 1;

    ep = struct( ...
        'signal',     sig_bp(idx_best), ...   % beta-bandpassed snippet
        'time',       time(idx_best),    ...
        'beta_power', maxPow             ...
    );
end

function ab = parseElec(name)
% Return [a b] electrode numbers from 'Ea_Eb', else [].
    tok = regexp(name, '^E(\d+)_E(\d+)$', 'tokens');
    if isempty(tok)
        ab = [];
    else
        ab = [str2double(tok{1}{1}), str2double(tok{1}{2})];
    end
end

function P = bandPower_(x, Fs, band)
% Integrated Welch PSD over [band(1) band(2)] Hz.
    win = round(Fs); ov = round(Fs/2);
    [Pxx, fr] = pwelch(x, win, ov, [], Fs);
    idx = (fr >= band(1)) & (fr <= band(2));
    P   = trapz(fr(idx), Pxx(idx));
end

function bp_signal = matlab_bpass_filt(signal, low_cut, high_cut, Fs)
% Zero-phase FIR bandpass filter (pure MATLAB).
%   Filter order = 3 * Fs / low_cut (rounded to nearest even number).
%   Signal is mirror-padded before filtfilt to suppress edge artefacts.

    signal   = signal(:);
    Wn       = [low_cut, high_cut] / (Fs / 2);
    bp_order = round(3 * Fs / low_cut);
    if mod(bp_order, 2) ~= 0; bp_order = bp_order + 1; end

    b = fir1(bp_order, Wn, 'bandpass', hamming(bp_order + 1));

    pad_samples = min(round(2 * Fs), floor((length(signal) - 1) / 2));
    signal_pad  = [flipud(signal(1:pad_samples)); signal; ...
                   flipud(signal(end - pad_samples + 1 : end))];

    bp_pad    = filtfilt(b, 1, signal_pad);
    bp_signal = bp_pad(pad_samples + 1 : end - pad_samples);
end