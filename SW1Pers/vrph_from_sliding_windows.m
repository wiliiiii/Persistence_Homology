function [dgms, Y] = vrph_from_sliding_windows(data_file_path, dopt, tauopt, varargin)
% ============================================================
% File: vrph_from_sliding_windows.m
% Author: Yuzhou He
% Email: ribosomehyz@gmail.com
%
% Description:
% Read a TXT time series and (d*, tau*), perform sliding window embedding,
% call ripser to compute Vietoris–Rips persistent homology (H0/H1/H2), and plot.
% GUI usage:
%   [dgms, Y] = vrph_from_sliding_windows('random_ts.txt', dopt, tauopt, ...
%       'SeriesIndex',1, 'Stride',1, 'coeff',2, 'maxdim',2, 'CapInf',true, 'Axes', ax);
% ============================================================

    % ---- Parse parameters ----
    ip = inputParser;
    addParameter(ip,'SeriesIndex',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip,'SeriesID',"",@(s)ischar(s)||isstring(s));
    addParameter(ip,'Stride',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip,'coeff',2,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
    addParameter(ip,'maxdim',2,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(ip,'CapInf',true,@(x)islogical(x)||ismember(x,[0 1]));
    addParameter(ip,'Axes',[],@(h) isempty(h) || isa(h,'matlab.graphics.axis.Axes'));
    parse(ip,varargin{:});
    P = ip.Results;

    % ---- 1) Sliding-window embedding ----
    [Y, meta] = sliding_windows_embedding(data_file_path, dopt, tauopt, ...
        'SeriesIndex',P.SeriesIndex, 'SeriesID',P.SeriesID, 'Stride',P.Stride);

    % ---- 2) Run ripser (exe in the same folder) ----
    dgms = run_ripser_vr(Y, P.maxdim, P.coeff);

    % ---- 3) Plot (robust legend & diagonal excluded from legend) ----
    ttl_suffix = sprintf(' (d*=%d, tau*=%d)', meta.dopt, meta.tauopt);
    if isempty(P.Axes)
        fig = figure('Color','w','Name','Persistence Diagram'); %#ok<NASGU>
        ax = axes();
    else
        ax = P.Axes; cla(ax);
    end

    if exist('plotPD','file') == 2
        % Use user-provided plotPD if available
        plotPD(ax, dgms, P.CapInf, ttl_suffix);
    else
        % ---- Safe plotting fallback ----
        hold(ax,'on'); axis(ax,'equal');

        % Extract diagrams for each homology dimension
        D0 = iGetDgm(dgms,0);
        D1 = iGetDgm(dgms,1);
        D2 = iGetDgm(dgms,2);

        % Optionally cap Inf deaths before computing axes
        if P.CapInf
            [D0c, D1c, D2c] = iCapInfAll(D0, D1, D2);
        else
            D0c = D0; D1c = D1; D2c = D2;
        end

        % Axis limits
        allXY = [];
        if ~isempty(D0c), allXY = [allXY; D0c]; end %#ok<AGROW>
        if ~isempty(D1c), allXY = [allXY; D1c]; end %#ok<AGROW>
        if ~isempty(D2c), allXY = [allXY; D2c]; end %#ok<AGROW>
        if isempty(allXY)
            lo = 0; hi = 1;
        else
            lo = min(allXY(:)); hi = max(allXY(:));
            if ~isfinite(lo) || ~isfinite(hi), lo = 0; hi = 1; end
            if hi <= lo, hi = lo + 1; end
        end

        % Diagonal (not in legend)
        plot(ax, [lo hi],[lo hi], 'k:', 'HandleVisibility','off');

        % Colors
        col0 = [0.00 0.45 0.74];  % H0
        col1 = [0.85 0.33 0.10];  % H1
        col2 = [0.30 0.30 0.30];  % H2

        % Scatter plots (only if non-empty)
        if ~isempty(D0c)
            scatter(ax, D0c(:,1), D0c(:,2), 24, 'filled', ...
                'MarkerFaceColor', col0, 'MarkerEdgeColor','none', 'DisplayName','H0');
        end
        if ~isempty(D1c)
            scatter(ax, D1c(:,1), D1c(:,2), 24, 'filled', ...
                'MarkerFaceColor', col1, 'MarkerEdgeColor','none', 'DisplayName','H1');
        end
        if ~isempty(D2c)
            scatter(ax, D2c(:,1), D2c(:,2), 24, 'filled', ...
                'MarkerFaceColor', col2, 'MarkerEdgeColor','none', 'DisplayName','H2');
        end

        % Legend with counts (ensure empty dims are still listed)
        n0 = size(D0c,1); n1 = size(D1c,1); n2 = size(D2c,1);

        if n0==0
            h0 = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor',col0, ...
                'MarkerEdgeColor','none', 'DisplayName',sprintf('H0 (%d)',n0));
        else
            h0 = findobj(ax,'-depth',1,'Type','Scatter','-and','DisplayName','H0');
            set(h0,'DisplayName',sprintf('H0 (%d)',n0));
        end
        if n1==0
            h1 = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor',col1, ...
                'MarkerEdgeColor','none', 'DisplayName',sprintf('H1 (%d)',n1));
        else
            h1 = findobj(ax,'-depth',1,'Type','Scatter','-and','DisplayName','H1');
            set(h1,'DisplayName',sprintf('H1 (%d)',n1));
        end
        if n2==0
            h2 = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor',col2, ...
                'MarkerEdgeColor','none', 'DisplayName',sprintf('H2 (%d)',n2));
        else
            h2 = findobj(ax,'-depth',1,'Type','Scatter','-and','DisplayName','H2');
            set(h2,'DisplayName',sprintf('H2 (%d)',n2));
        end

        legend(ax, [h0 h1 h2], ...
            {get(h0,'DisplayName'), get(h1,'DisplayName'), get(h2,'DisplayName')}, ...
            'Location','southeast');

        xlim(ax,[lo hi]); ylim(ax,[lo hi]);
        grid(ax,'on'); box(ax,'on');
        title(ax, ['Persistence Diagram (H0/H1/H2)' ttl_suffix]);
        xlabel(ax,'birth'); ylabel(ax,'death');
        hold(ax,'off');
    end
end

% ================= Helpers =================

function dgms = run_ripser_vr(Y, maxdim, coeff)
    % Find ripser executable next to this .m; error if not found
    local_dir = fileparts(mfilename('fullpath'));
    cand = { fullfile(local_dir,'ripser_static.exe'), ...
             fullfile(local_dir,'ripser.exe'), ...
             'ripser_static.exe', 'ripser.exe' };
    ripser_path = '';
    for i = 1:numel(cand)
        if exist(cand{i}, 'file'), ripser_path = cand{i}; break; end
    end
    if isempty(ripser_path)
        error('Ripser executable not found. Place ripser_static.exe or ripser.exe in the same folder as this .m.');
    end

    % Write point cloud (one point per line). Let ripser compute Euclidean distances.
    tmpDir = tempname; mkdir(tmpDir);
    c = onCleanup(@() (exist(tmpDir,'dir')==7 && rmdir(tmpDir,'s'))); %#ok<NASGU>
    pcFile = fullfile(tmpDir,'pointcloud.txt');
    try
        writematrix(Y, pcFile, 'Delimiter',' ');
    catch
        dlmwrite(pcFile, Y, 'delimiter',' ');
    end

    % Call ripser with point-cloud input
    if coeff ~= 2
        warning('Current ripser build does not support arbitrary coefficient fields; using Z/2.');
    end
    cmd = sprintf('"%s" --format point-cloud --dim %d "%s"', ripser_path, maxdim, pcFile);

    [status, out] = system(cmd);
    if status ~= 0
        error('Ripser call failed: %s', out);
    end

    % Parse ripser stdout
    dgms = parse_ripser_stdout(out, maxdim);
end

function dgms = parse_ripser_stdout(s, maxdim)
    % Parse blocks "persistence intervals in dim k"
    dgms = cell(1, maxdim+1);
    for k = 0:maxdim, dgms{k+1} = zeros(0,2); end
    lines = regexp(s,'\r?\n','split');
    dim = -1; data = zeros(0,2);
    for i = 1:numel(lines)
        L = strtrim(lines{i});
        m = regexp(L, '^persistence intervals in dim (\d+):$', 'tokens');
        if ~isempty(m)
            if dim>=0, dgms{dim+1} = data; end
            dim = str2double(m{1}{1}); data = zeros(0,2); continue;
        end
        t = regexp(L,'^\[\s*([^\s,]+)\s*,\s*([^\s\)]+)\s*\)','tokens');
        if ~isempty(t)
            a = str2double(t{1}{1}); b = str2double(t{1}{2});
            if ~isfinite(b), b = inf; end
            data(end+1,:) = [a,b]; %#ok<AGROW>
        end
    end
    if dim>=0, dgms{dim+1} = data; end
end

% -------- Small utilities --------
function D = iGetDgm(dgms, k)
    if numel(dgms) >= k+1
        D = dgms{k+1};
    else
        D = zeros(0,2);
    end
end

function [D0c, D1c, D2c] = iCapInfAll(D0, D1, D2)
    % Cap death=Inf to a common upper bound across dims for stable axes/legend
    allFiniteMax = -inf;
    if ~isempty(D0), allFiniteMax = max(allFiniteMax, max(D0(isfinite(D0(:,2)),2))); end
    if ~isempty(D1), allFiniteMax = max(allFiniteMax, max(D1(isfinite(D1(:,2)),2))); end
    if ~isempty(D2), allFiniteMax = max(allFiniteMax, max(D2(isfinite(D2(:,2)),2))); end
    if ~isfinite(allFiniteMax), allFiniteMax = 0; end
    cap = allFiniteMax*1.05 + 1e-6;

    D0c = iCapOne(D0, cap);
    D1c = iCapOne(D1, cap);
    D2c = iCapOne(D2, cap);
end

function D = iCapOne(D, cap)
    if isempty(D), return; end
    infMask = ~isfinite(D(:,2));
    if any(infMask)
        D(infMask,2) = cap;
    end
end
