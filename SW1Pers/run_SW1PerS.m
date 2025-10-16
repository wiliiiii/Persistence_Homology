function [score, details] = run_SW1PerS(y, varargin)
% ============================================================
% File: run_SW1PerS.m
% Author: Yuzhou He
% Email: ribosomehyz@gmail.com
%
% Description:
% Friendly wrapper around SW1PerS core (positional-only).
% Supports three call styles:
%   [score,details] = run_SW1PerS(y, paramsStruct)
%   [score,details] = run_SW1PerS(y, 'num_cycles',2, 'feature_type',3, ...)
%   [score,details] = run_SW1PerS(y, num_cycles, feature_type, num_points, ...
%                                 allow_trending, ...
%                                 use_meanshift, meanshift_epsilon, ...
%                                 use_expaverage, expaverage_alpha, ...
%                                 use_movingaverage, movingaverage_window, ...
%                                 use_smoothspline, smoothspline_sigma)
% ============================================================

    % ===== 1) Locate core and avoid shadowing =====
    % Clear same-named variables that might shadow the function
    if exist('SW1PerS_v1','var'); clear SW1PerS_v1; end
    if evalin('base','exist(''SW1PerS_v1'',''var'')'); evalin('base','clear(''SW1PerS_v1'')'); end

    % Try common folders if core is not on path
    if isempty(which('SW1PerS_v1')) && isempty(which('SW1PerS_v0')) && isempty(which('SW1PerS'))
        here = fileparts(mfilename('fullpath'));
        cands = { here, fullfile(here,'sw1pers'), fullfile(here,'..','sw1pers'), fullfile(here,'../sw1pers') };
        for i = 1:numel(cands)
            if exist(fullfile(cands{i},'SW1PerS_v1.m'),'file') || ...
               exist(fullfile(cands{i},'SW1PerS_v0.m'),'file') || ...
               exist(fullfile(cands{i},'SW1PerS.m'),'file')
                addpath(cands{i});
                break;
            end
        end
    end

    % Pick an available core (prefer v1)
    if ~isempty(which('SW1PerS_v1'))
        coreFun = @SW1PerS_v1;
    elseif ~isempty(which('SW1PerS_v0'))
        coreFun = @SW1PerS_v0;
    elseif ~isempty(which('SW1PerS'))
        coreFun = @SW1PerS;
    else
        error(['Core not found: SW1PerS_v1.m / SW1PerS_v0.m / SW1PerS.m. ', ...
               'Add the sw1pers/ folder to path or place a core file next to run_SW1PerS.m.']);
    end

    % ===== 2) Basic checks =====
    if ~isvector(y) || ~isnumeric(y)
        error('run_SW1PerS:BadSignal', 'Input y must be a numeric vector.');
    end
    y = y(:)';  % enforce row vector

    % ===== 3) Defaults (aligned with reference settings) =====
    P.num_cycles           = 2;
    P.feature_type         = 3;
    P.num_points           = 200;
    P.allow_trending       = true;

    P.use_meanshift        = true;
    P.meanshift_epsilon    = 1 - cos(pi/16);

    P.use_expaverage       = false;
    P.expaverage_alpha     = NaN;

    P.use_movingaverage    = true;
    P.movingaverage_window = 10;

    P.use_smoothspline     = false;
    P.smoothspline_sigma   = NaN;

    % ===== 4) Parse varargin (struct / 12-position / name-value) =====
    if numel(varargin) == 1 && isstruct(varargin{1})
        % (A) struct-style
        S = varargin{1};
        P = mergeStruct(P, S);

    elseif numel(varargin) == 12 && all(~cellfun(@ischar, varargin)) && all(~cellfun(@isstring, varargin))
        % (B) 12 positional args
        keys = {'num_cycles','feature_type','num_points','allow_trending', ...
                'use_meanshift','meanshift_epsilon', ...
                'use_expaverage','expaverage_alpha', ...
                'use_movingaverage','movingaverage_window', ...
                'use_smoothspline','smoothspline_sigma'};
        for i = 1:12
            P.(keys{i}) = varargin{i};
        end

    else
        % (C) name-value pairs
        if mod(numel(varargin),2) ~= 0
            error('run_SW1PerS:BadArgs', 'Name-value pairs must be provided in an even number.');
        end
        for i = 1:2:numel(varargin)
            k = normalizeKey(varargin{i});
            v = varargin{i+1};
            if isfield(P, k)
                P.(k) = v;
            else
                error('run_SW1PerS:UnknownParam', 'Unknown parameter name: %s', string(varargin{i}));
            end
        end
    end

    % ===== 5) Call core (positional-only) =====
    score = coreFun( ...
        y, ...
        P.num_cycles, ...
        P.feature_type, ...
        P.num_points, ...
        P.allow_trending, ...
        P.use_meanshift, P.meanshift_epsilon, ...
        P.use_expaverage, P.expaverage_alpha, ...
        P.use_movingaverage, P.movingaverage_window, ...
        P.use_smoothspline, P.smoothspline_sigma);

    % ===== 6) Optional details =====
    if nargout > 1
        details = P;
        details.signal_length = numel(y);
    end
end

% ===== helpers =====
function P = mergeStruct(P, S)
    fS = fieldnames(S);
    for i = 1:numel(fS)
        k = normalizeKey(fS{i});
        if isfield(P, k)
            P.(k) = S.(fS{i});
        end
    end
end

function k = normalizeKey(k)
    % Normalize keys like 'NumPoints' / 'num_points' / 'numpoints'
    if isstring(k); k = char(k); end
    if ~ischar(k); error('Parameter names must be char or string.'); end
    k = lower(strrep(k, '-', '_'));
    k = strrep(k, ' ', '_');
    aliases = containers.Map( ...
        {'numpoints','movingaveragewindow','meanshiftepsilon', ...
         'expaveragealpha','smoothsplinesigma'}, ...
        {'num_points','movingaverage_window','meanshift_epsilon', ...
         'expaverage_alpha','smoothspline_sigma'});
    k_nosym = regexprep(k, '_', '');
    if isKey(aliases, k_nosym)
        k = aliases(k_nosym);
    end
end
