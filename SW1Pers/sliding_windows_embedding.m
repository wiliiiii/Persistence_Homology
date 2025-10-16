function [Y, meta] = sliding_windows_embedding(data_file_path, dopt, tauopt, varargin)
% ============================================================
% File: sliding_windows_embedding.m
% Author: Yuzhou He
% Email: ribosomehyz@gmail.com
%
% Description:
% Read a TXT produced by create_random_time_series (1st row = time ticks;
% subsequent rows = "<id>\t v1 \t v2 ..."), select a series, and build a
% sliding-window embedding with given optimal (d*, tau*). Return the m×d*
% point cloud Y and metadata.
%
% Usage:
%   [Y, meta] = sliding_windows_embedding('random_ts.txt', dopt, tauopt, ...
%                     'SeriesIndex', 1, 'Stride', 1);
%
% Name-Value options:
%   'SeriesIndex'  : 1-based row index after the header (default 1)
%   'SeriesID'     : string id to select a specific row (overrides SeriesIndex)
%   'Stride'       : window stride (default 1; same notion as in Takens search)
%   'Standardize'  : z-score the series (default false)
%   'Center'       : mean-center only (default false; if true, overrides Standardize)
%
% Returns:
%   Y    : m × dopt embedded point cloud; m = 1 + floor((n - 1 - (dopt-1)*tauopt)/Stride)
%   meta : struct with file path, chosen id, n, m, d*, tau*, stride
% ============================================================

    % ---- Parse inputs ----
    ip = inputParser;
    addParameter(ip,'SeriesIndex',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip,'SeriesID',"",@(s)ischar(s)||isstring(s));
    addParameter(ip,'Stride',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip,'Standardize',false,@islogical);
    addParameter(ip,'Center',false,@islogical);
    parse(ip,varargin{:});
    P = ip.Results;

    % ---- Basic checks ----
    if ~exist(data_file_path,'file')
        error('sliding_windows_embedding: File not found: %s', data_file_path);
    end
    if ~(isscalar(dopt) && dopt>=1 && dopt==round(dopt))
        error('sliding_windows_embedding: dopt must be a positive integer.');
    end
    if ~(isscalar(tauopt) && tauopt>=1 && tauopt==round(tauopt))
        error('sliding_windows_embedding: tauopt must be a positive integer.');
    end

    % ---- Read TXT: skip first line (time); each subsequent line: id \t values ----
    fid = fopen(data_file_path,'r');
    header = fgetl(fid); %#ok<NASGU>
    lines = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
    fclose(fid);
    lines = lines{1};
    if isempty(lines)
        error('sliding_windows_embedding: No series rows found in file.');
    end

    % ---- Select target line by SeriesID (if provided) or SeriesIndex ----
    targetLine = '';
    chosenID = "";
    if strlength(P.SeriesID) > 0
        for i = 1:numel(lines)
            parts = split(string(lines{i}), sprintf('\t'));
            if parts(1) == string(P.SeriesID)
                targetLine = lines{i};
                chosenID = parts(1);
                break;
            end
        end
        if isempty(targetLine)
            error('sliding_windows_embedding: SeriesID="%s" not found.', string(P.SeriesID));
        end
    else
        idx = P.SeriesIndex;
        if idx > numel(lines)
            error('sliding_windows_embedding: SeriesIndex=%d out of range (total %d).', idx, numel(lines));
        end
        targetLine = lines{idx};
        chosenID = split(string(targetLine), sprintf('\t')); 
        chosenID = chosenID(1);
    end

    % ---- Parse the selected row ----
    parts = split(string(targetLine), sprintf('\t'));
    if numel(parts) < 2
        error('sliding_windows_embedding: Malformed row: %s', targetLine);
    end
    vals = str2double(parts(2:end));
    if any(~isfinite(vals))
        error('sliding_windows_embedding: Series contains NaN/Inf.');
    end
    y = vals(:);
    n = numel(y);

    % ---- Optional pre-processing ----
    if P.Center
        y = y - mean(y);
    elseif P.Standardize
        mu = mean(y); sg = std(y);
        if sg == 0, sg = 1; end
        y = (y - mu) / sg;
    end

    % ---- Sliding-window embedding ----
    % Window length = (d-1)*tau + 1
    win = (dopt-1)*tauopt + 1;
    if n < win
        error(['sliding_windows_embedding: Series length n=%d is smaller than window length win=%d ', ...
               '(determined by d*=%d, tau*=%d). Increase n or decrease d*/tau*.'], n, win, dopt, tauopt);
    end
    stride = P.Stride;
    starts = 1:stride:(n - win + 1);
    m = numel(starts);

    Y = zeros(m, dopt);
    for i = 1:m
        s = starts(i);
        idxs = s + (0:(dopt-1)) * tauopt;   % y[s + (0:d-1)*tau]
        Y(i,:) = y(idxs).';
    end

    % ---- Metadata ----
    meta = struct( ...
        'data_file_path', data_file_path, ...
        'series_id', chosenID, ...
        'n', n, 'm', m, ...
        'dopt', dopt, 'tauopt', tauopt, ...
        'stride', stride ...
    );
end
