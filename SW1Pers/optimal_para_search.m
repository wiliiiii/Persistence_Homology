function [Y_embedded, tau_opt, d_opt] = optimal_para_search(y, max_dim, max_delay, stride)
% ============================================================
% File: optimal_para_search.m
% Author: Yuzhou He
% Email: ribosomehyz@gmail.com
%
% Description:
% Search optimal time delay (tau) via first local minimum of mutual
% information, and optimal embedding dimension (d) via simplified FNN
% (Kennel ’92). Return Takens-embedded matrix with downsampling by stride.
%
% Usage:
%   [Y_embedded, tau_opt, d_opt] = optimal_para_search(y, max_dim, max_delay, stride)
%
% Notes:
%   Feasibility condition: n - (d-1)*tau >= 1. The search space is clipped
%   internally to avoid infeasible settings.
% ============================================================

    % ---- Preprocess & defaults ----
    y = y(:);
    n = numel(y);
    if nargin < 4 || isempty(stride),    stride    = 1;  end
    if nargin < 3 || isempty(max_delay), max_delay = min(30, max(1, floor((n-1)/2))); end
    if nargin < 2 || isempty(max_dim),   max_dim   = 30; end

    % ---- 1) Choose tau by first local minimum of I(tau) ----
    % Ensure at least d>=3 is feasible; clip tau upper bound
    tau_max_feas = max(1, floor((n-1)/(3-1)));  % max tau s.t. d=3 feasible
    tau_upper    = min(max_delay, tau_max_feas);
    if tau_upper < 1
        error('Sequence too short: cannot embed with d>=3. Increase length or lower max_delay/max_dim.');
    end

    nbins = max(16, round(sqrt(n)));     % histogram bins (heuristic)
    Ivals = zeros(tau_upper,1);
    for tau = 1:tau_upper
        Ivals(tau) = mutual_information_hist(y, tau, nbins);
    end
    tau_opt = first_local_minimum(Ivals);
    if isempty(tau_opt), [~,tau_opt] = min(Ivals); end

    % ---- 2) Choose d by simplified FNN ----
    d_upper_feas = min(max_dim, floor((n-1)/tau_opt) + 1);
    if d_upper_feas < 2
        error('At tau=%d the sequence is too short for d>=2.', tau_opt);
    end
    R_th    = 10;    % FNN ratio threshold (Kennel)
    fnn_tol = 0.02;  % acceptable false-neighbor rate (2%)
    d_opt   = choose_dim_by_fnn(y, tau_opt, d_upper_feas, R_th, fnn_tol);

    % ---- 3) Final embedding with stride ----
    Y_full     = takens_embed(y, d_opt, tau_opt);   % full sampling
    Y_embedded = Y_full(1:stride:end, :);           % downsample

    % ---- Print (Python-style) ----
    fprintf('Shape of embedded time series: [%d, %d]\n', size(Y_embedded,1), size(Y_embedded,2));
    fprintf('Optimal embedding dimension is %d and time delay is %d\n', d_opt, tau_opt);
end


% ================== Helpers ==================

function I = mutual_information_hist(y, tau, nbins)
% Histogram-based MI: I(X;Y) = sum p(x,y) log ( p(x,y) / (p(x)p(y)) )
    y1 = y(1:end-tau);
    y2 = y(1+tau:end);

    ymin = min(y); ymax = max(y);
    if ymin == ymax
        I = 0; return;
    end
    edges = linspace(ymin, ymax, nbins+1);

    px  = histcounts(y1, edges, 'Normalization','probability');
    py  = histcounts(y2, edges, 'Normalization','probability');
    pxy = histcounts2(y1, y2, edges, edges, 'Normalization','probability');

    [ix, iy] = find(pxy > 0);
    vals = zeros(numel(ix),1);
    for k = 1:numel(ix)
        j = ix(k); l = iy(k);
        denom = px(j)*py(l) + eps;
        vals(k) = pxy(j,l) * log( pxy(j,l)/denom + eps );
    end
    I = sum(vals);
end

function idx = first_local_minimum(v)
% First strict local minimum index (1-based); empty if none
    idx = [];
    if numel(v) < 3, return; end
    for t = 2:numel(v)-1
        if v(t-1) > v(t) && v(t) <= v(t+1)
            idx = t; return;
        end
    end
end

function Y = takens_embed(y, d, tau)
% Classic Takens embedding (no stride)
% Output size: (n - (d-1)*tau) x d
    y = y(:);
    n = numel(y);
    m = n - (d-1)*tau;
    if m < 1
        error('Takens embedding infeasible: n - (d-1)*tau < 1 (n=%d, d=%d, tau=%d)', n, d, tau);
    end
    Y = zeros(m, d);
    for j = 1:d
        Y(:,j) = y( (1:m) + (j-1)*tau );
    end
end

function d_opt = choose_dim_by_fnn(y, tau, d_upper, R_th, fnn_tol)
% Simplified FNN (Kennel ’92): pick smallest d with FNN ratio < fnn_tol.
% Uses Euclidean distances; nearest neighbor via knnsearch if available,
% otherwise a pure-MATLAB O(m^2) fallback.

    y = y(:);
    n = numel(y); %#ok<NASGU>

    % Precompute max-d embedding once; reuse slices
    Ymax = takens_embed(y, d_upper, tau);   % [m x d_upper]
    m = size(Ymax,1);

    getY  = @(d)  Ymax(:, 1:d);
    getY1 = @(d)  Ymax(:, 1:(d+1));

    % Nearest neighbor (exclude self)
    function [idx_nn, dist_nn] = nearest_neighbor(P)
        try
            [idx_nn, dist_nn] = knnsearch(P, P, 'K', 2);
            idx_nn  = idx_nn(:,2);          % second is the true NN
            dist_nn = dist_nn(:,2);
        catch
            D2 = squareform_pd(P);
            D2(1:m+1:end) = inf;            % remove diagonal
            [dist_nn, idx_nn] = min(D2, [], 2);
        end
    end

    d_opt = d_upper;   % fallback
    for d = 1:(d_upper-1)
        Pd   = getY(d);
        Pd1  = getY1(d);

        [~,  dist_d]  = nearest_neighbor(Pd);
        [~, dist_d1]  = nearest_neighbor(Pd1);

        % Kennel ratio: R = ||x^{d+1}_i - x^{d+1}_{j(i)}|| / ||x^{d}_i - x^{d}_{j(i)}||
        nume = dist_d1;
        deno = max(dist_d, eps);
        R    = nume ./ deno;

        fnn_ratio = mean(R > R_th);
        if fnn_ratio < fnn_tol
            d_opt = d + 1;
            break;
        end
    end

    d_opt = max(2, d_opt);  % ensure at least 2
end

function D = squareform_pd(P)
% Pairwise Euclidean distance matrix (simple O(m^2))
    m = size(P,1);
    D = zeros(m,m);
    for i = 1:m
        Pi = P(i,:);
        dif = P - Pi;
        D(i,:) = sqrt(sum(dif.*dif, 2)).';
    end
end
