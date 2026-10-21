function START
% START  Launcher menu for all analyses in this project.
%
%   Usage:
%     1) Open this file in MATLAB and press Run, OR type
%            START
%        in the command window and press Enter.
%     2) A menu opens. Press "Run" on a card - figure windows pop up and
%        the text output appears at the bottom of the menu and in the
%        MATLAB command window.
%
%   The menu adds all project folders to the MATLAB path by itself.

    rootDir = fileparts(mfilename('fullpath'));

    % ---------- ANALYSIS LIST ----------
    A = struct('name',{},'folder',{},'file',{},'desc',{},'cat',{});
    A(end+1) = item('DIC - coarse correlation', '6_DIC', 'dic_coarse.m', ...
        'Fast integer-pixel DIC (CPU): displacement and strain fields, reference vs. last frame of the series. ~0.5-2 min.', 'Experiment');
    A(end+1) = item('DIC - subpixel engine', '6_DIC', 'dic_icgn.m', ...
        'Subpixel DIC engine (IC-GN + ZNSSD + splines + PLS strains) with a built-in self-test; processes the full loading series. ~1-3 min.', 'Experiment');
    A(end+1) = item('DIC - validation vs DICe', '9_Validation_DICe', 'validate_dice.m', ...
        'One-button validation: compares dic_icgn.m with DICe (Sandia) official gold results on the identical images. Self-contained.', 'Experiment');
    A(end+1) = item('DIC - Challenge 2.0 benchmark', '10_Challenge2', 'challenge2_benchmark.m', ...
        'Noise floor, spatial resolution and MEI on the DIC Challenge 2.0 Star images, vs published codes (batched IC-GN, all CPU cores). Tens of minutes; see 10_Challenge2/README.md.', 'Experiment');

    % add project folders to the MATLAB path
    for k = 1:numel(A)
        p = fullfile(rootDir, A(k).folder);
        if isfolder(p), addpath(p); end
    end

    % ---------- PALETTE (dark theme) ----------
    col.bg      = [0.106 0.110 0.165];
    col.panel   = [0.157 0.165 0.235];
    col.card    = [0.196 0.204 0.290];
    col.accent  = [0.478 0.635 0.965];
    col.txt     = [0.871 0.894 0.992];
    col.sub     = [0.620 0.660 0.800];
    col.console = [0.070 0.072 0.110];
    catColor = containers.Map( ...
        {'Plates','Strength','Dynamics','Experiment','Geometry'}, ...
        {[0.31 0.55 0.93],[0.91 0.50 0.36],[0.40 0.73 0.52],[0.66 0.47 0.90],[0.91 0.73 0.36]});

    % ---------- WINDOW ----------
    fig = uifigure('Name','Engineering Analysis Lab - MATLAB', ...
        'Position',[80 60 1120 760],'Color',col.bg,'HandleVisibility','off');

    root = uigridlayout(fig,[3 1]);
    root.RowHeight   = {96,'1x',180};
    root.ColumnWidth = {'1x'};
    root.Padding     = [18 18 18 18];
    root.RowSpacing  = 14;
    root.BackgroundColor = col.bg;

    % ----- HEADER -----
    head = uipanel(root,'BackgroundColor',col.panel,'BorderType','none');
    hgl = uigridlayout(head,[2 2]);
    hgl.RowHeight = {'1x','1x'};
    hgl.ColumnWidth = {'1x',170};
    hgl.Padding = [22 12 18 12];
    hgl.ColumnSpacing = 10;
    hgl.BackgroundColor = col.panel;

    lblT = uilabel(hgl,'Text','Engineering Analysis Lab', ...
        'FontSize',24,'FontWeight','bold','FontColor',col.txt);
    lblT.Layout.Row = 1; lblT.Layout.Column = 1;
    lblS = uilabel(hgl,'Text','Pick an analysis - figures open in separate windows, text output appears below.', ...
        'FontSize',13,'FontColor',col.sub);
    lblS.Layout.Row = 2; lblS.Layout.Column = 1;

    bFolder = uibutton(hgl,'Text',' Project folder','FontSize',13,'FontWeight','bold', ...
        'BackgroundColor',col.card,'FontColor',col.txt, ...
        'ButtonPushedFcn',@(~,~) openFolder(rootDir));
    bFolder.Layout.Row = 1; bFolder.Layout.Column = 2;

    % ----- ANALYSIS CARDS -----
    nRows = ceil(numel(A)/2);
    center = uipanel(root,'BackgroundColor',col.bg,'BorderType','none');
    cgl = uigridlayout(center,[nRows 2]);
    cgl.Scrollable = 'on';
    cgl.BackgroundColor = col.bg;
    cgl.Padding = [0 0 0 0];
    cgl.RowSpacing = 12; cgl.ColumnSpacing = 12;
    cgl.RowHeight = repmat({104},1,nRows);
    cgl.ColumnWidth = {'1x','1x'};

    runBtns = gobjects(1, numel(A));
    for k = 1:numel(A)
        c = catColor(A(k).cat);
        card = uipanel(cgl,'BackgroundColor',col.card,'BorderType','none');
        card.Layout.Row = ceil(k/2);
        card.Layout.Column = 2 - mod(k,2);   % odd k -> 1, even k -> 2

        gl = uigridlayout(card,[2 2]);
        gl.RowHeight   = {'fit','1x'};
        gl.ColumnWidth = {'1x',150};
        gl.Padding = [16 12 16 12];
        gl.RowSpacing = 6;
        gl.BackgroundColor = col.card;

        nm = uilabel(gl,'Text',A(k).name,'FontSize',16,'FontWeight','bold','FontColor',col.txt);
        nm.Layout.Row = 1; nm.Layout.Column = 1;

        ds = uilabel(gl,'Text',A(k).desc,'FontSize',11.5,'FontColor',col.sub,'WordWrap','on');
        ds.Layout.Row = 2; ds.Layout.Column = [1 2];

        btn = uibutton(gl,'Text','Run','FontSize',14,'FontWeight','bold', ...
            'FontColor',[1 1 1],'BackgroundColor',c, ...
            'Tooltip',['Category: ' A(k).cat '   |   file: ' A(k).folder '/' A(k).file], ...
            'ButtonPushedFcn',@(~,~) onRun(fig, A(k)));
        btn.Layout.Row = 1; btn.Layout.Column = 2;
        runBtns(k) = btn;
    end

    % ----- BOTTOM PANEL: status + console -----
    bottom = uipanel(root,'BackgroundColor',col.panel,'BorderType','none');
    bgl = uigridlayout(bottom,[2 1]);
    bgl.RowHeight = {'fit','1x'};
    bgl.ColumnWidth = {'1x'};
    bgl.Padding = [16 10 16 14];
    bgl.RowSpacing = 8;
    bgl.BackgroundColor = col.panel;

    status = uilabel(bgl,'Text','Ready. Pick an analysis above.', ...
        'FontSize',14,'FontWeight','bold','FontColor',col.accent);
    status.Layout.Row = 1;

    console = uitextarea(bgl,'Editable','off','FontName',get(groot,'FixedWidthFontName'),'FontSize',12.5, ...
        'BackgroundColor',col.console,'FontColor',[0.80 0.86 0.78], ...
        'Value',{'- text output of the selected analysis appears here -'});
    console.Layout.Row = 2;

    % store handles and context
    S.status  = status;
    S.console = console;
    S.rootDir = rootDir;
    S.runBtns = runBtns;
    fig.UserData = S;
end

% ===================== HELPERS =====================

function it = item(name, folder, file, desc, cat)
    it = struct('name',name,'folder',folder,'file',file,'desc',desc,'cat',cat);
end

function onRun(fig, a)
    S = fig.UserData;
    scriptPath = fullfile(S.rootDir, a.folder, a.file);

    if ~isfile(scriptPath)
        S.status.Text = ['File not found: ' a.folder '/' a.file];
        S.status.FontColor = [0.95 0.45 0.45];
        return;
    end

    S.status.Text = ['Running: ' a.name '  ...  (the window may freeze while computing)'];
    S.status.FontColor = [0.95 0.80 0.45];
    S.console.Value = {'Computing, please wait...'};
    % run buttons are disabled during computation (double-click guard);
    % onCleanup re-enables them on ANY exit, including errors and Ctrl+C
    set(S.runBtns, 'Enable', 'off');
    restoreBtns = onCleanup(@() set(S.runBtns, 'Enable', 'on')); %#ok<NASGU>

    drawnow;

    try
        t0 = tic;
        out = i_runscript(scriptPath);
        elapsed = toc(t0);
        if isempty(strtrim(out))
            out = sprintf('"%s" finished.\nSee the opened figure windows for results.', a.name);
        end
        S.console.Value = cellstr(splitlines(string(out)));
        S.status.Text = sprintf('OK: Done: %s  (%.1f s)', a.name, elapsed);
        S.status.FontColor = [0.55 0.85 0.55];
    catch err
        msg = [{['ERROR while running "' a.name '":']}; {''}; {err.message}];
        if ~isempty(err.stack)
            msg = [msg; {''}; {'Stack:'}; cellstr(string({err.stack.name}'))];
        end
        if contains(err.message,'multicuboid') || contains(err.message,'createpde') ...
                || contains(err.message,'structural')
            msg = [msg; {''}; {'The PDE Toolbox seems to be missing - this analysis requires it.'}];
        end
        S.console.Value = msg;
        S.status.Text = ['Error: ' a.name '  (details in the output box)'];
        S.status.FontColor = [0.95 0.45 0.45];
    end
end

function out = i_runscript(scriptPath)
    % Isolated run: close all / clc / clear inside a script do not affect
    % the menu window (it is created with HandleVisibility='off').
    % rehash refreshes MATLAB's path/file cache in case files were changed
    % outside of MATLAB (sync tools, file managers, etc.).
    rehash;
    out = evalc('run(scriptPath)');
end

function openFolder(rootDir)
    if ispc
        winopen(rootDir);
    elseif ismac
        system(['open "' rootDir '"']);
    else
        system(['xdg-open "' rootDir '"']);
    end
end
