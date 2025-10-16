function TDA_TimeSeries_GUI
% ============================================================
% File: TDA_TimeSeries_GUI.m
% Author: Yuzhou He
% Email: ribosomehyz@gmail.com
%
% Description:
% Minimal GUI to (1) generate/load a time series, (2) search Takens params
% and compute VR persistent homology (H0/H1/H2) with ripser, and (3) compute
% SW1PerS score. Allows passing stride/coeff to the pipeline.
%
% Dependencies:
%   create_random_time_series.m, optimal_para_search.m,
%   vrph_from_sliding_windows.m, run_SW1PerS.m, computePH.m, plotPD.m
% ============================================================

    defaultFont = 'Microsoft YaHei UI';

    % ---- Main figure ----
    f = figure('Name','TDA Time-Series (Minimal GUI)', ...
               'Color','w','NumberTitle','off', ...
               'Position',[80 80 1200 680]);

    % ---- Left control panel ----
    pnl = uipanel(f,'Title','Controls','FontWeight','bold', ...
        'BackgroundColor',[0.96 0.97 1.00], ...
        'Position',[0.02 0.05 0.30 0.90],'FontName',defaultFont);

    % === Sub-panel: giotto-tda input ===
    pnlGiotto = uipanel(pnl,'Title','giotto-tda Input', ...
        'FontWeight','bold','BackgroundColor',pnl.BackgroundColor, ...
        'Position',[0.03 0.62 0.94 0.37],'FontName',defaultFont);

    % ==== File path ====
    uicontrol(pnlGiotto,'Style','text','String','TXT path','BackgroundColor',pnl.BackgroundColor, ...
        'HorizontalAlignment','left','Position',[10 175 80 20],'FontName',defaultFont);
    edtPath = uicontrol(pnlGiotto,'Style','edit','String','random_ts.txt', ...
        'Position',[90 175 180 26],'FontName',defaultFont);
    uicontrol(pnlGiotto,'Style','pushbutton','String','Browse...','Position',[275 175 70 26], ...
        'Callback',@onBrowse,'FontName',defaultFont);

    % ==== Random series params ====
    uicontrol(pnlGiotto,'Style','text','String','n (points)','BackgroundColor',pnl.BackgroundColor, ...
        'HorizontalAlignment','left','Position',[10 125 60 20],'FontName',defaultFont);
    edtN = uicontrol(pnlGiotto,'Style','edit','String','800','Position',[70 125 60 24],'FontName',defaultFont);

    uicontrol(pnlGiotto,'Style','text','String','range [min,max]','BackgroundColor',pnl.BackgroundColor, ...
        'HorizontalAlignment','left','Position',[140 125 110 20],'FontName',defaultFont);
    edtMin = uicontrol(pnlGiotto,'Style','edit','String','-1','Position',[250 125 45 24],'FontName',defaultFont);
    edtMax = uicontrol(pnlGiotto,'Style','edit','String','1','Position',[300 125 45 24],'FontName',defaultFont);

    % (kept hidden as per your previous logic)
    uicontrol(pnlGiotto,'Style','text','String','NumSeries','BackgroundColor',pnl.BackgroundColor, ...
        'HorizontalAlignment','left','Position',[10 61 70 20],'FontName',defaultFont,'Visible','off');
    edtNumSeries = uicontrol(pnlGiotto,'Style','edit','String','1','Position',[80 57 50 24], ...
        'FontName',defaultFont,'Visible','off');

    uicontrol(pnlGiotto,'Style','text','String','Step (tick of row-1)','BackgroundColor',pnl.BackgroundColor, ...
        'HorizontalAlignment','left','Position',[140 61 110 20],'FontName',defaultFont,'Visible','off');
    edtStep = uicontrol(pnlGiotto,'Style','edit','String','4','Position',[250 57 45 24], ...
        'FontName',defaultFont,'Visible','off');

    uicontrol(pnlGiotto,'Style','text','String','Series index','BackgroundColor',pnl.BackgroundColor, ...
        'HorizontalAlignment','left','Position',[10 31 70 20],'FontName',defaultFont,'Visible','off');
    edtSeriesIdx = uicontrol(pnlGiotto,'Style','edit','String','1','Position',[80 27 50 24], ...
        'FontName',defaultFont,'Visible','off');

    % ==== Takens search & sliding window params ====
    uicontrol(pnlGiotto,'Style','text','String','max_dim','BackgroundColor',pnl.BackgroundColor, ...
        'HorizontalAlignment','left','Position',[10 90 60 20],'FontName',defaultFont);
    edtMaxDim = uicontrol(pnlGiotto,'Style','edit','String','30','Position',[70 90 50 24],'FontName',defaultFont);

    uicontrol(pnlGiotto,'Style','text','String','max_delay','BackgroundColor',pnl.BackgroundColor, ...
        'HorizontalAlignment','left','Position',[130 90 70 20],'FontName',defaultFont);
    edtMaxDelay = uicontrol(pnlGiotto,'Style','edit','String','30','Position',[200 90 50 24],'FontName',defaultFont);

    uicontrol(pnlGiotto,'Style','text','String','stride','BackgroundColor',pnl.BackgroundColor, ...
        'HorizontalAlignment','left','Position',[10 55 50 20],'FontName',defaultFont);
    edtStride = uicontrol(pnlGiotto,'Style','edit','String','4','Position',[70 55 45 24],'FontName',defaultFont);

    % ==== VRPH option (hidden, value=1) ====
    chkCapInf = uicontrol(pnlGiotto,'Style','checkbox','String','Cap Inf','Value',1, ...
        'BackgroundColor',pnl.BackgroundColor,'Position',[230 -30 110 22], ...
        'FontName',defaultFont,'Visible','off');

    % === Sub-panel: SW1PerS input ===
    pnlSw1 = uipanel(pnl,'Title','SW1PerS Input', ...
        'FontWeight','bold','BackgroundColor',pnl.BackgroundColor, ...
        'Position',[0.03 0.32 0.94 0.30],'FontName',defaultFont);

    uicontrol(pnlSw1,'Style','text','String','num_cycles','BackgroundColor',pnl.BackgroundColor, ...
        'HorizontalAlignment','left','Position',[10 130 80 20],'FontName',defaultFont);
    edtNumCycles = uicontrol(pnlSw1,'Style','edit','String','2','Position',[90 130 50 24],'FontName',defaultFont);

    uicontrol(pnlSw1,'Style','text','String','feature_type','BackgroundColor',pnl.BackgroundColor, ...
        'HorizontalAlignment','left','Position',[150 130 80 20],'FontName',defaultFont);
    edtFeatType = uicontrol(pnlSw1,'Style','edit','String','3','Position',[230 130 50 24],'FontName',defaultFont);

    uicontrol(pnlSw1,'Style','text','String','num_points','BackgroundColor',pnl.BackgroundColor, ...
        'HorizontalAlignment','left','Position',[10 100 80 20],'FontName',defaultFont);
    edtNumPoints = uicontrol(pnlSw1,'Style','edit','String','200','Position',[90 100 50 24],'FontName',defaultFont);

    uicontrol(pnlSw1,'Style','text','String','movingavg win','BackgroundColor',pnl.BackgroundColor, ...
        'HorizontalAlignment','left','Position',[150 100 90 20],'FontName',defaultFont);
    edtMAWin = uicontrol(pnlSw1,'Style','edit','String','10','Position',[240 100 40 24],'FontName',defaultFont);

    chkAllowTrend = uicontrol(pnlSw1,'Style','checkbox','String','allow_trending','Value',1, ...
        'BackgroundColor',pnl.BackgroundColor,'Position',[10 60 120 20],'FontName',defaultFont);
    chkMeanShift  = uicontrol(pnlSw1,'Style','checkbox','String','use_meanshift','Value',1, ...
        'BackgroundColor',pnl.BackgroundColor,'Position',[140 60 110 20],'FontName',defaultFont);
    uicontrol(pnlSw1,'Style','text','String','ms epsilon','BackgroundColor',pnl.BackgroundColor, ...
        'HorizontalAlignment','left','Position',[260 60 100 20],'FontName',defaultFont);
    edtMSEps = uicontrol(pnlSw1,'Style','edit','String',num2str(1 - cos(pi/16),'%.8f'), ...
        'Position',[340 60 60 24],'FontName',defaultFont);

    chkExpAvg = uicontrol(pnlSw1,'Style','checkbox','String','expavg','Value',0, ...
        'BackgroundColor',pnl.BackgroundColor,'Position',[10 16 70 20],'FontName',defaultFont);
    edtExpAlpha = uicontrol(pnlSw1,'Style','edit','String','NaN','Position',[80 12 50 24],'FontName',defaultFont);

    chkMAvg = uicontrol(pnlSw1,'Style','checkbox','String','movingavg','Value',1, ...
        'BackgroundColor',pnl.BackgroundColor,'Position',[140 16 90 20],'FontName',defaultFont);

    chkSSpline = uicontrol(pnlSw1,'Style','checkbox','String','smoothspline','Value',0, ...
        'BackgroundColor',pnl.BackgroundColor,'Position',[240 16 100 20],'FontName',defaultFont);
    edtSSigma   = uicontrol(pnlSw1,'Style','edit','String','NaN','Position',[340 12 30 24],'FontName',defaultFont);

    % === Sub-panel: actions ===
    pnlOps = uipanel(pnl,'Title','Actions', ...
        'FontWeight','bold','BackgroundColor',pnl.BackgroundColor, ...
        'Position',[0.03 0.06 0.94 0.25],'FontName',defaultFont);

    uicontrol(pnlOps,'Style','pushbutton','String','① Generate random time series', ...
        'Position',[10 90 360 30],'FontName',defaultFont, 'Callback',@onGenerateTS);

    uicontrol(pnlOps,'Style','pushbutton','String','② Compute Persistence Diagram (Takens → VR)', ...
        'Position',[10 50 360 30],'FontName',defaultFont, 'Callback',@onComputePD);

    uicontrol(pnlOps,'Style','pushbutton','String','SW1PerS Score', ...
        'Position', [10 10 360 30],'FontName',defaultFont, 'Callback',@onSw1PerSBtn);

    % Status text
    txtInfo = uicontrol(pnl,'Style','text','String','Status: Ready','BackgroundColor',pnl.BackgroundColor, ...
        'HorizontalAlignment','left','Position',[10 5 360 20],'FontName',defaultFont);

    % ---- Right plots ----
    axTS = axes(f,'Position',[0.36 0.58 0.28 0.36]); title(axTS,'Random Time Series');
    axPD = axes(f,'Position',[0.67 0.58 0.28 0.36]); title(axPD,'Persistence Diagram (H0/H1/H2)'); axis(axPD,'equal');
    uit  = uitable(f,'Position',[0.36 0.08 0.59 0.42],'Data',cell(0,3), ...
                   'ColumnName',{'id','score','per'},'ColumnEditable',[false,false,false]);

    % ---- Bottom-right hints ----
    txtInstr = uicontrol(f,'Style','text', ...
        'String',['Hints:', newline, ...
                  '1) Choose an existing TXT via Browse, or click (1) to generate a random series.', newline, ...
                  '2) Set max_dim / max_delay / stride; click (2) to run Takens → VR and plot PD.', newline, ...
                  '3) Configure SW1PerS params; click "SW1PerS Score" to compute.'], ...
        'BackgroundColor','w', ...
        'HorizontalAlignment','left', ...
        'FontName',defaultFont, ...
        'FontSize',12, ...
        'Position',[550 100 400 250]); %#ok<NASGU>

    % =================== Callbacks ===================
    function onBrowse(~,~)
        [fn,fp] = uigetfile('*.txt','Choose TXT file');
        if isequal(fn,0), return; end
        edtPath.String = fullfile(fp,fn);
    end

    % ========= (1) Generate random time series =========
    function onGenerateTS(~,~)
        try
            fname = char(edtPath.String);
            n     = str2double(edtN.String);
            tmin  = str2double(edtMin.String);
            tmax  = str2double(edtMax.String);

            if any(isnan([n tmin tmax]))
                error('NaN found in parameters.');
            end
            if tmax <= tmin
                error('Require t_max > t_min.');
            end

            outdir = fileparts(fname); if isempty(outdir), outdir = "."; end
            if exist(outdir,'dir')==0, mkdir(outdir); end
            delete(fullfile(outdir,'*.txt'));

            figs_before = findall(0,'Type','figure');

            create_random_time_series(fname, n, 1, [tmin, tmax], ...
                'allow_trending',  (chkAllowTrend.Value==1), ...
                'use_meanshift',   (chkMeanShift.Value==1));

            % Close any extra figs created by generator
            figs_after = findall(0,'Type','figure');
            new_figs = setdiff(figs_after, figs_before);
            new_figs(new_figs==f) = [];
            delete(new_figs);

            % Reformat TXT to 2 lines: time row, then "1 \t series"
            A = readmatrix(fname);
            if isempty(A), error('Generated file is empty: %s', fname); end
            if isvector(A), y = A(:); else, y = A(:,1); end
            t = linspace(tmin, tmax, numel(y)).';

            fid = fopen(fname,'w');
            fprintf(fid,'%.16g', t(1));
            for i = 2:numel(t), fprintf(fid,'\t%.16g', t(i)); end
            fprintf(fid,'\n');
            fprintf(fid,'1');
            for i = 1:numel(y), fprintf(fid,'\t%.16g', y(i)); end
            fprintf(fid,'\n');
            fclose(fid);

            cla(axTS); plot(axTS, t, y, '.-'); grid(axTS,'on');
            title(axTS, sprintf('Random Series (n=%d, t∈[%.3g, %.3g])', numel(y), tmin, tmax));
            xlabel(axTS,'t'); ylabel(axTS,'y');

            txtInfo.String = sprintf('Status: generated → %s', fname);
        catch ME
            txtInfo.String = ['Status: generation failed - ' ME.message];
            errordlg(ME.message,'Generation failed');
        end
    end

    % ========= (2) Compute PD (Takens → VR via ripser) =========
    function onComputePD(~,~)
        try
            fname     = char(edtPath.String);
            max_dim   = str2double(edtMaxDim.String);
            max_delay = str2double(edtMaxDelay.String);
            stride    = str2double(edtStride.String);
            coeff     = 2;    % ripser build assumed Z/2
            vrmaxdim  = 2;
            capInf    = (chkCapInf.Value==1);

            if any(isnan([max_dim max_delay stride coeff vrmaxdim]))
                error('Invalid numeric input.');
            end

            y = readmatrix(fname); %#ok<NASGU>
            if isempty(y), error('TXT is empty or unreadable: %s', fname); end

            y = readmatrix(fname);
            y = y(:);

            [~, tau_opt, d_opt] = optimal_para_search(y, max_dim, max_delay, stride);

            cla(axPD);
            [~, ~] = vrph_from_sliding_windows(fname, d_opt, tau_opt, ...
                'Stride',      stride, ...
                'coeff',       2, ...
                'maxdim',      2, ...
                'CapInf',      capInf, ...
                'Axes',        axPD);

            txtInfo.String = sprintf('Status: done. d*=%d, tau*=%d, stride=%d, coeff=%d(actual=2), maxdim=%d', ...
                                      d_opt, tau_opt, stride, coeff, vrmaxdim);
            drawnow;
        catch ME
            txtInfo.String = ['Status: error - ' ME.message];
            errordlg(ME.getReport('basic','hyperlinks','off'), 'PD computation error');
        end
    end

    % ========= (3) SW1PerS =========
    function onSw1PerSBtn(~,~)
        try
            txtPath = strtrim(edtPath.String);
            if isempty(txtPath) || exist(txtPath,'file')~=2
                error('TXT not found: %s', txtPath);
            end
            y = read_series_from_txt(txtPath, 1);
            y = double(y(:));

            stride      = str2double(edtStride.String); %#ok<NASGU>
            coeff       = 2; %#ok<NASGU>
            maxdim      = 2; %#ok<NASGU>
            featureType = str2double(edtFeatType.String);
            numCycles   = str2double(edtNumCycles.String);
            numPoints   = str2double(edtNumPoints.String);
            movavgWin   = str2double(edtMAWin.String);

            allowTrend  = logical(chkAllowTrend.Value);
            useMS       = logical(chkMeanShift.Value);
            msEps       = str2double(edtMSEps.String);

            if logical(chkExpAvg.Value)
                expavg = str2double(strtrim(edtExpAlpha.String));
                if isnan(expavg), expavg = []; end
            else
                expavg = [];
            end

            if logical(chkMAvg.Value)
                movingavg = movavgWin;
            else
                movingavg = [];
            end

            if logical(chkSSpline.Value)
                smoothspline = str2double(strtrim(edtSSigma.String));
                if isnan(smoothspline), smoothspline = []; end
            else
                smoothspline = [];
            end

            score = SW1PerS_v1( ...
                y, ...
                numCycles, ...
                featureType, ...
                numPoints, ...
                allowTrend, ...
                useMS, ...
                msEps, ...
                ~isempty(expavg), ...
                expavg, ...
                ~isempty(movingavg), ...
                movingavg, ...
                ~isempty(smoothspline), ...
                smoothspline ...
            );

            series_id = '1';
            per_val   = numCycles;
            set(uit, 'Data', {series_id, score, per_val});
            txtInfo.String = sprintf('Status: SW1PerS done, score = %.4f', score);
        catch ME
            txtInfo.String = ['Status: SW1PerS failed - ' ME.message];
            errordlg(ME.message,'SW1PerS failed','modal');
        end
    end
end

% ======== Helper: read the idx-th series from TXT (skip 1st time row) ========
function y = read_series_from_txt(fname, idx)
    if ~exist(fname,'file')
        error('TXT not found: %s', fname);
    end
    fid = fopen(fname,'r');
    header = fgetl(fid); %#ok<NASGU>
    lines = textscan(fid,'%s','delimiter','\n','whitespace','');
    fclose(fid);
    lines = lines{1};
    if idx < 1 || idx > numel(lines)
        error('SeriesIndex=%d out of range (total %d lines).', idx, numel(lines));
    end
    parts = split(string(lines{idx}), sprintf('\t'));
    if numel(parts) < 2
        error('Malformed line: %s', lines{idx});
    end
    vals = str2double(parts(2:end));
    if any(~isfinite(vals)), error('Series contains NaN/Inf'); end
    y = vals(:);
end
