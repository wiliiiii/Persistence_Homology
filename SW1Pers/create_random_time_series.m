function create_random_time_series(output_file, n_points, n_series, t_range, varargin)
% ============================================================
% File: create_random_time_series.m
% Author: Yuzhou He
% Email: ribosomehyz@gmail.com
%
% Description:
% Generate one or more random time series y = f(t) over a given time range
% and save them to a TXT file.
%
% Usage:
%   create_random_time_series('random_ts.txt', 50, 3, [-5 5]);
%
% Name-Value options:
%   'allow_trending'   (logical, default=false)  - add mild linear trend
%   'use_meanshift'    (logical, default=true)   - shift mean to zero
%   'amplitude_range'  ([a_min,a_max], default=[0.5,1.5])
%   'frequency_range'  ([f_min,f_max], default=[0.5,2])
%   'noise_std'        (scalar >=0, default=0.1)
% ============================================================

    % ---- Parse options ----
    ip = inputParser;
    addParameter(ip,'allow_trending',false,@islogical);
    addParameter(ip,'use_meanshift',true,@islogical);
    addParameter(ip,'amplitude_range',[0.5,1.5],@(v)isnumeric(v)&&numel(v)==2);
    addParameter(ip,'frequency_range',[0.5,2],@(v)isnumeric(v)&&numel(v)==2);
    addParameter(ip,'noise_std',0.1,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    parse(ip,varargin{:});
    P = ip.Results;

    % ---- Time vector ----
    t_min = t_range(1);
    t_max = t_range(2);
    t = linspace(t_min, t_max, n_points)';   % column vector

    % ---- Synthesize series ----
    Y = zeros(n_points, n_series);
    for s = 1:n_series
        % Sum 2–4 sinusoidal components with random amps/freqs/phases
        k = randi([2,4]);
        y = zeros(size(t));
        for j = 1:k
            A   = rand_range(P.amplitude_range);
            f   = rand_range(P.frequency_range);
            phi = 2*pi*rand;
            y   = y + A * cos(2*pi*f*t + phi);
        end
        % Optional trend & noise
        if P.allow_trending
            y = y + (randn*0.1)*t;
        end
        if P.noise_std > 0
            y = y + P.noise_std*randn(size(y));
        end
        if P.use_meanshift
            y = y - mean(y);
        end
        Y(:,s) = y;
    end

    % ---- Save TXT (matrix of size n_points × n_series) ----
    writematrix(Y, output_file, 'Delimiter', ' ');
    fprintf('[OK] Generated random time series %s (t ∈ [%.2f, %.2f], n=%d)\n', ...
            output_file, t_min, t_max, n_points);

    % ---- Quick plot ----
    figure('Color','w','Name','Random Time Series');
    plot(t, Y, 'LineWidth',1.1);
    xlabel('t'); ylabel('y');
    title(sprintf('Random Time Series (n=%d, range=[%.2f, %.2f])', n_points, t_min, t_max));
    grid on;
end

function v = rand_range(r)
    v = r(1) + rand*(r(2)-r(1));
end
