function [GmObj, StateStr, InputStr, OutputStr] = ImpedanceDataModel(ApparatusBus, Type, PowerFlow, Para, Ts, ListBus)
%IMPEDANCEDATAMODEL Build a 2-by-2 admittance model from impedance data.
% The impedance matrix is inverted at each frequency, then all four
% admittance channels are identified together so that they share one pole set.
%
% Types 4000 and 4001 represent black-box GFL and GFM devices, respectively.
% Para.ImpedanceDataType selects the input data domain independently.
%
% ImpedanceDataType 1 (dq) columns:
%   Frequency, Zdd_mag, Zdd_phase, Zdq_mag, Zdq_phase,
%   Zqd_mag, Zqd_phase, Zqq_mag, Zqq_phase
% ImpedanceDataType 2 uses the same layout for Zpp, Zpn, Znp, and Znn and is
% converted from the sequence domain to the dq domain before identification.
%
% Required parameters:
%   Para.ImpedanceDataFile - Path to the impedance-data workbook.
%   Para.ImpedanceDataType - 1: dq domain; 2: sequence domain.
% Optional parameters:
%   Para.FitInitialOrder   - Initial order for directed search.
%   Para.FitOrderRange     - Candidate orders when no initial order is given.

    if ~isfield(Para,'ImpedanceDataType') || isempty(Para.ImpedanceDataType)
        error('SimplusGT:MissingImpedanceDataType', ...
            ['Para.ImpedanceDataType is required (1: dq, 2: sequence). ' ...
             'Types 4000 and 4001 now identify GFL and GFM devices, ' ...
             'not the impedance-data domain. Update the configuration.']);
    end
    validateattributes(Para.ImpedanceDataType,{'numeric'}, ...
        {'real','scalar','finite','integer','>=',1,'<=',2}, ...
        mfilename,'Para.ImpedanceDataType');

    %% Load impedance data
    if ~isfield(Para, 'ImpedanceData') || isempty(Para.ImpedanceData)
        if isfield(Para, 'ImpedanceDataFile') && ~isempty(Para.ImpedanceDataFile)
            impedanceFile = resolveImpedanceFile(Para.ImpedanceDataFile);
            fprintf('Reading impedance data file: %s\n', impedanceFile);
            try
                numData = xlsread(impedanceFile);
                if isempty(numData)
                    error('The impedance data file is empty or has an invalid format.');
                end
                Para.ImpedanceData.Frequency_Hz   = numData(:, 1);
                Para.ImpedanceData.Zdd_mag        = numData(:, 2);
                Para.ImpedanceData.Zdd_phase_deg  = numData(:, 3);
                Para.ImpedanceData.Zdq_mag        = numData(:, 4);
                Para.ImpedanceData.Zdq_phase_deg  = numData(:, 5);
                Para.ImpedanceData.Zqd_mag        = numData(:, 6);
                Para.ImpedanceData.Zqd_phase_deg  = numData(:, 7);
                Para.ImpedanceData.Zqq_mag        = numData(:, 8);
                Para.ImpedanceData.Zqq_phase_deg  = numData(:, 9);
            catch ME
                error('Failed to read the impedance data file: %s', ME.message);
            end
        else
            error('Para.ImpedanceDataFile is required.');
        end
    end

    % Convert sequence-domain data for either device type to the dq domain.
    if Para.ImpedanceDataType == 2
        fprintf('  Sequence impedance: converting pn to dq...\n');

        % The column positions match dq data, but their meanings differ.
        f_vec       = Para.ImpedanceData.Frequency_Hz;
        Zpp_mag     = Para.ImpedanceData.Zdd_mag;
        Zpp_phase   = Para.ImpedanceData.Zdd_phase_deg;
        Zpn_mag     = Para.ImpedanceData.Zdq_mag;
        Zpn_phase   = Para.ImpedanceData.Zdq_phase_deg;
        Znp_mag     = Para.ImpedanceData.Zqd_mag;
        Znp_phase   = Para.ImpedanceData.Zqd_phase_deg;
        Znn_mag     = Para.ImpedanceData.Zqq_mag;
        Znn_phase   = Para.ImpedanceData.Zqq_phase_deg;

        Zpp = Zpp_mag .* exp(1j * deg2rad(Zpp_phase));
        Zpn = Zpn_mag .* exp(1j * deg2rad(Zpn_phase));
        Znp = Znp_mag .* exp(1j * deg2rad(Znp_phase));
        Znn = Znn_mag .* exp(1j * deg2rad(Znn_phase));

        % Inverse sequence-to-dq transformation:
        %   Z_pn(f) = T * Z_dq(f - f0) * inv(T)
        %   Z_dq(f - f0) = inv(T) * Z_pn(f) * T
        %     z_dd = (z_pp + z_nn + z_pn + z_np) / 2
        %     z_qq = (z_pp + z_nn - z_pn - z_np) / 2
        %     z_qd = -j * (z_pp - z_nn + z_pn - z_np) / 2
        %     z_dq = -j * (z_pn - z_np - z_pp + z_nn) / 2
        Zdd = (Zpp + Znn + Zpn + Znp) / 2;
        Zqq = (Zpp + Znn - Zpn - Znp) / 2;
        Zqd = -1j * (Zpp - Znn + Zpn - Znp) / 2;
        Zdq = -1j * (Zpn - Znp - Zpp + Znn) / 2;

        if isfield(Para,'Fbase') && ~isempty(Para.Fbase)
            f0 = Para.Fbase;
        else
            f0 = 50;
        end

        % Keep only samples that remain at positive frequency after shifting.
        idx_keep = f_vec > f0;
        f_vec = f_vec(idx_keep);

        Zdd = Zdd(idx_keep);
        Zdq = Zdq(idx_keep);
        Zqd = Zqd(idx_keep);
        Zqq = Zqq(idx_keep);

        Para.ImpedanceData.Frequency_Hz = f_vec - f0;
        Para.ImpedanceData.Zdd_mag        = abs(Zdd);
        Para.ImpedanceData.Zdd_phase_deg  = rad2deg(angle(Zdd));
        Para.ImpedanceData.Zdq_mag        = abs(Zdq);
        Para.ImpedanceData.Zdq_phase_deg  = rad2deg(angle(Zdq));
        Para.ImpedanceData.Zqd_mag        = abs(Zqd);
        Para.ImpedanceData.Zqd_phase_deg  = rad2deg(angle(Zqd));
        Para.ImpedanceData.Zqq_mag        = abs(Zqq);
        Para.ImpedanceData.Zqq_phase_deg  = rad2deg(angle(Zqq));

        fprintf('  Sequence-to-dq conversion completed.\n');
    end

    %% Convert impedance to admittance at each frequency
    f_Hz = Para.ImpedanceData.Frequency_Hz;
    w    = 2 * pi * f_Hz;
    N_f  = length(w);

    Zdd = Para.ImpedanceData.Zdd_mag .* exp(1j * deg2rad(Para.ImpedanceData.Zdd_phase_deg));
    Zdq = Para.ImpedanceData.Zdq_mag .* exp(1j * deg2rad(Para.ImpedanceData.Zdq_phase_deg));
    Zqd = Para.ImpedanceData.Zqd_mag .* exp(1j * deg2rad(Para.ImpedanceData.Zqd_phase_deg));
    Zqq = Para.ImpedanceData.Zqq_mag .* exp(1j * deg2rad(Para.ImpedanceData.Zqq_phase_deg));

    resp_MIMO = zeros(2, 2, N_f);
    for k = 1:N_f
        Z_k = [Zdd(k), Zdq(k); Zqd(k), Zqq(k)];
        Y_k = Z_k \ eye(2);
        resp_MIMO(:, :, k) = Y_k;
    end
    frdata_Y_MIMO = frd(resp_MIMO, w);

    %% Identify a MIMO state-space model
    fprintf('Starting MIMO state-space identification...\n');
    opt = ssestOptions('Display', 'off');

    max_safe_order = floor(N_f / 2);

    SEARCH_STEP       = 2;
    PLATEAU_THRESHOLD = 99.5;

    % Use a directed search when an initial order is available.
    if isfield(Para, 'FitInitialOrder') && ~isempty(Para.FitInitialOrder) ...
                                       && Para.FitInitialOrder > 0

        nx_initial = Para.FitInitialOrder;
        if nx_initial > max_safe_order
            nx_initial = max_safe_order;
        end
        if nx_initial < 2
            nx_initial = 2;
        end

        fprintf('  Directed search: initial=%d, step=%d, threshold=%.1f%%\n', ...
                nx_initial, SEARCH_STEP, PLATEAU_THRESHOLD);

        try
            sys_initial = ssest(frdata_Y_MIMO, nx_initial, opt);
            fit_initial = mean(sys_initial.Report.Fit.FitPercent(:));
            fprintf('    Initial order %2d: fit=%.2f%%\n', nx_initial, fit_initial);
        catch ME
            error('Identification failed at initial order %d: %s', ...
                nx_initial, ME.message);
        end

        best_order = nx_initial;
        best_fit   = fit_initial;
        best_sys   = sys_initial;

        if fit_initial < PLATEAU_THRESHOLD
            fprintf('    Fit is below the threshold; searching upward...\n');
            nx = best_order;
            while nx + SEARCH_STEP <= max_safe_order
                try
                    sys_tmp = ssest(frdata_Y_MIMO, nx + SEARCH_STEP, opt);
                    fit_tmp = mean(sys_tmp.Report.Fit.FitPercent(:));
                    delta = fit_tmp - best_fit;
                    fprintf('      %2d: %.2f%% (delta=%.2f%%)', ...
                        nx + SEARCH_STEP, fit_tmp, delta);

                    nx = nx + SEARCH_STEP;

                    if fit_tmp >= PLATEAU_THRESHOLD
                        best_order = nx;
                        best_fit   = fit_tmp;
                        best_sys   = sys_tmp;
                        fprintf(' [threshold reached]\n');
                        break;
                    elseif fit_tmp > best_fit
                        best_order = nx;
                        best_fit   = fit_tmp;
                        best_sys   = sys_tmp;
                        fprintf(' [improved]\n');
                    else
                        fprintf('\n');
                    end
                catch ME
                    fprintf(' [failed] %s\n', ME.message);
                    break;
                end
            end
        end

        fprintf('  Selected order %d with %.2f%% fit.\n', best_order, best_fit);
        sys_best = best_sys;

    else
        % Scan candidate orders when no initial order is provided.
        if isfield(Para, 'FitOrderRange') && ~isempty(Para.FitOrderRange)
            sys_order_range = Para.FitOrderRange;
        else
            sys_order_range = [4:4:20, 28:8:48];
        end
        sys_order_range(sys_order_range > max_safe_order) = [];
        sys_order_range(sys_order_range < 2) = [];
        if isempty(sys_order_range)
            sys_order_range = 4:4:max_safe_order;
        end

        n_orders  = length(sys_order_range);
        fit_vals  = -inf(n_orders, 1);
        sys_cells = cell(n_orders, 1);

        fprintf('  Automatic scan: %d candidates [%s]\n', n_orders, ...
                strjoin(arrayfun(@num2str, sys_order_range, 'UniformOutput', false), ', '));

        for idx = 1:n_orders
            nx = sys_order_range(idx);
            try
                sys_tmp = ssest(frdata_Y_MIMO, nx, opt);
                fit_vals(idx) = mean(sys_tmp.Report.Fit.FitPercent(:));
                sys_cells{idx} = sys_tmp;
                fprintf('    %2d: %.2f%%\n', nx, fit_vals(idx));
            catch ME
                fprintf('    %2d: failed (%s)\n', nx, ME.message);
            end
        end

        valid_idx = find(fit_vals > -inf);
        if isempty(valid_idx)
            error('ssest identification failed for every candidate order.');
        end

        [max_fit, ~] = max(fit_vals(valid_idx));
        plateau_idx = valid_idx(fit_vals(valid_idx) >= max_fit - 0.1);
        selected_idx = plateau_idx(1);

        best_order = sys_order_range(selected_idx);
        best_fit   = fit_vals(selected_idx);

        if best_fit < 99.5
            warning(['The best fit is only %.2f%%. The data quality or ' ...
                'frequency range may be insufficient.'], best_fit);
        end

        fprintf('  Selected order %d with %.2f%% fit.\n', best_order, best_fit);
        sys_best = sys_cells{selected_idx};
    end

    % Refit once at the selected order.
    try
        sys_fit = ssest(frdata_Y_MIMO, best_order, opt);
    catch
        warning('Refitting failed; using the cached model.');
        sys_fit = sys_best;
    end

    GmObj = ss(sys_fit);
    GmObj = minreal(GmObj);

    %% Assign signal and state names
    InputStr  = {'v_d', 'v_q'};
    OutputStr = {'i_d', 'i_q'};
    try
        GmObj.InputName  = InputStr;
        GmObj.OutputName = OutputStr;
    catch
    end

    n = size(GmObj.A, 1);
    StateStr = arrayfun(@(k) sprintf('x_%d', k), 1:n, 'UniformOutput', false);
    try
        GmObj.StateName = StateStr;
    catch
    end

    fprintf('  Final model: order %d after minreal.\n', n);
end

function resolvedFile = resolveImpedanceFile(configuredFile)
% First honor the configured absolute/relative path.  If it no longer
% exists (for example Excel2Json expanded a bare name against the project
% root), search the project tree by the same file name.
    configuredFile = char(configuredFile);
    if exist(configuredFile, 'file') == 2
        resolvedFile = configuredFile;
        return;
    end

    [~,baseName,extension] = fileparts(configuredFile);
    targetName = [baseName,extension];
    classDir = fileparts(mfilename('fullpath'));
    packageDir = fileparts(classDir);
    projectRoot = fileparts(packageDir);
    matches = dir(fullfile(projectRoot,'**',targetName));
    matches = matches(~[matches.isdir]);

    if isempty(matches)
        error(['The impedance data file does not exist: %s. ' ...
            'No file named %s was found in the project tree.'], ...
            configuredFile,targetName);
    end

    matchPaths = arrayfun(@(item) fullfile(item.folder,item.name), ...
        matches,'UniformOutput',false);
    matchPaths = unique(matchPaths,'stable');
    if numel(matchPaths) > 1
        error(['Multiple impedance data files named %s were found. ' ...
            'Specify an unambiguous path in the configuration:\n%s'], ...
            targetName,strjoin(matchPaths,newline));
    end

    resolvedFile = matchPaths{1};
    fprintf(['The configured path was not found. Using the impedance data ' ...
        'file located in the project tree: %s\n'],resolvedFile);
end
