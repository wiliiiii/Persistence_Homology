function score = SW1PerS_v1(...
    signal, ...
    num_cycles, ...
    feature_type, ...
    num_points, ... 
    allow_trending, ...
    use_meanshift, ... 
    meanshift_epsilon, ... 
    use_expaverage, ...
    expaverage_alpha, ...
    use_movingaverage, ... 
    movingaverage_window, ...
    use_smoothspline, ...
    smoothspline_sigma ...
    )
% ============================================================
% File: SW1PerS_v1.m
% Author: Yuzhou He
% Email: ribosomehyz@gmail.com
%
% Description:
% Sliding Windows + 1-Persistence Scoring (SW1PerS).
% Given an evenly sampled 1D signal, build a sliding-window point cloud,
% (optionally) denoise/smooth, then compute a 1D persistence-based score.
%
% Inputs (key):
%   signal              - 1D signal, evenly sampled (row/column vector)
%   num_cycles          - number of cycles to search (e.g., 2 for daily rhythm)
%   feature_type        - integer in {1..5}; recommended: 3
%   num_points          - number of samples (windows) in the point cloud
%   allow_trending      - if true, center each point by its own mean
%   use_meanshift       - if true, mean-shift smoothing on the cloud
%   meanshift_epsilon   - mean-shift bandwidth in cosine distance
%                         (typical: 1 - cos(pi/16))
%   use_expaverage      - if true, exponential moving average denoising
%   expaverage_alpha    - EMA parameter in (0,1)
%   use_movingaverage   - if true, simple moving average denoising
%   movingaverage_window- moving average window size (e.g., ~length/5)
%   use_smoothspline    - if true, smoothing spline denoising on signal
%   smoothspline_sigma  - smoothing spline scale parameter
%
% Output:
%   score               - SW1PerS score (higher = more periodic)
% ============================================================

    nS = length(signal);
    signal = reshape(signal,1,[]);

    %% Part II — Parameters
    N   = 7;                         % highest Fourier harmonic captured
    p   = 11;                        % prime for Z/p; choose p > N
    M   = 2 * N;                     % window dimension = M+1
    tau = (2*pi) / ((M+1) * num_cycles);   % step; window width = M*tau

    %% Part III — Signal pre-processing

    % Exponential moving average (low-pass)
    if use_expaverage
        sigLowPass = signal;
        for i = 2:nS
            sigLowPass(i) = sigLowPass(i-1) + expaverage_alpha * (signal(i) - sigLowPass(i-1));
        end
        signal = sigLowPass;
    end

    % Simple moving average
    if use_movingaverage
        signal = smooth(signal', movingaverage_window, 'moving')';
    end

    % Detrending (disabled by default)
    use_detrending = 0;
    if (use_detrending == 1)
        if (detrending_type == 1)
            % piecewise linear detrend
            signal = detrend(signal, floor(nS/num_cycles));
        elseif (detrending_type == 2)
            % global linear detrend
            signal = detrend(signal);
        end
    end

    %% Part IV — Sliding-window point cloud

    t  = 2*pi*linspace(0,1,nS);
    T  = (2*pi - M*tau) * linspace(0,1,num_points);
    tt = tau*(0:M)'*ones(1,num_points) + ones(M+1,1)*T;

    if use_smoothspline
        ss.weights = ones(1,nS);
        tolerance  = nS * (smoothspline_sigma * peak2peak(signal)).^2;
        [sp, vals] = spaps(t, signal, tolerance, ss.weights);
        cloud_raw  = fnval(sp, tt);
        signal     = vals;
    else
        cloud_raw = spline(t, signal, tt);
    end

    % Mean-centering
    if allow_trending
        cloud_centered = cloud_raw - ones(M+1,1) * mean(cloud_raw);
    else
        cloud_centered = cloud_raw - mean(signal);
    end

    % Point-wise normalization
    SW_cloud = cloud_centered ./ (ones(M+1,1) * sqrt(sum(cloud_centered.^2)));
    SW_cloud = SW_cloud';

    %% Part V — Cloud post-processing

    % Mean-shift smoothing (cosine distance)
    if use_meanshift
        cloud = zeros(size(SW_cloud));
        D     = squareform(pdist(SW_cloud,'cosine'));
        indD  = (D <= meanshift_epsilon);
        for k = 1:num_points
            cloud(k,:) = mean(SW_cloud(indD(k,:),:));
        end
        SW_cloud = cloud;
    end

    %% Part VI — 1D persistence & score

    % Java TDA (ripser-like) interface
    javaclasspath('tda/jars/tda_1.0.0_unionFind3.jar');
    import shell.*

    tda = Tda();
    tda.RCA1({'settingsFile=tda/cts.txt', ['zp_value=',num2str(p)], 'distanceBoundOnEdges=2'}, SW_cloud);

    intervals0 = tda.getResultsRCA1(0).getIntervals;
    intervals1 = tda.getResultsRCA1(1).getIntervals;

    if (isempty(intervals0) || isempty(intervals1))
        score = 1;
    else
        % Feature variants
        if (feature_type == 1)
            % 1 - max length(H1)/sqrt(3)
            score = 1 - max(intervals1(:,2) - intervals1(:,1)) / sqrt(3);
        elseif (feature_type == 2)
            % 1 - max(b^2 - a * max(H0_death)) / 3
            score = 1 - max(intervals1(:,2).^2 - intervals1(:,1) * max(intervals0(:,2))) / 3;
        elseif (feature_type == 3)
            % 1 - max(b^2 - a^2) / 3   (default)
            score = 1 - max(intervals1(:,2).^2 - intervals1(:,1).^2) / 3;
        elseif (feature_type == 4)
            % 1 - max(b^3 - a^2) / (3*sqrt(3))
            score = 1 - max(intervals1(:,2).^3 - intervals1(:,1).^2) / (3*sqrt(3));
        elseif (feature_type == 5)
            % 1 - max(b^4 - a^2) / 9
            score = 1 - max(intervals1(:,2).^4 - intervals1(:,1).^2) / 9;
        else
            % Fallback to type-3 if invalid input
            score = 1 - max(intervals1(:,2).^2 - intervals1(:,1).^2) / 3;
        end
    end
end
