classdef PostProcessor
    % PostProcessor
    % Export hydrostatic and GZ results and generate publication-style plots.
    %
    % Main improvements:
    %   - clearer axis labels with proper physical meaning
    %   - more polished and consistent plotting style
    %   - improved legends, annotations, and zero-reference handling
    %   - safer extraction of derived hydrostatic quantities
    %   - richer hydrostatics figure with better visual hierarchy

    properties
        cfg struct
    end

    methods
        function obj = PostProcessor(cfg)
            if nargin < 1 || isempty(cfg)
                cfg = struct();
            end
            obj.cfg = cfg;
        end

        function generate(obj, outputTag, results, mesh, inputData)
            COND = results.COND;

            obj.ensureOutputFolder();

            obj.writeHydroDat(outputTag, COND, inputData.VCG);
            obj.writeHydroTxt(outputTag, COND, inputData.VCG);

            matFile = fullfile('Output', ['Hydro_', outputTag]);
            cfg = obj.cfg;
            save(matFile, 'COND', 'mesh', 'inputData', 'results', 'cfg');
            fprintf('  Saved:   %s.mat\n', matFile);

            fprintf('  Generating hydrostatic plots ...\n');
            obj.plotHydrostaticsBasic(outputTag, COND, inputData.VCG);
            fprintf('  Hydrostatic plots done.\n');

            if ~isempty(results.GZdata)
                obj.writeGZDat(outputTag, results.GZdata);
                obj.writeGZTxt(outputTag, results.GZdata);
                obj.plotGZ(outputTag, results.GZdata);
            end
        end
    end

    methods (Access = private)
        function ensureOutputFolder(~)
            if exist('Output', 'dir') ~= 7
                mkdir('Output');
            end
        end

        function writeHydroDat(obj, outputTag, COND, vcg)
            datFile = fullfile('Output', ['Hydrostatics_', outputTag, '.dat']);
            fid = fopen(datFile, 'w');
            if fid == -1
                warning('PostProcessor:FileOpen', 'Could not open %s for writing.', datFile);
                return;
            end
            cleaner = onCleanup(@() fclose(fid));

            fprintf(fid, '# T[m]\t V[m3]\t XB[m]\t YB[m]\t ZB[m]\t Ix[m4]\t Iy[m4]\t BMt[m]\t BMl[m]\t KMt[m]\t GMt[m]\t GMl[m]\t MCT[tm/cm]\t Af[m2]\t Lwl[m]\t Bwl[m]\t Ams[m2]\t Ws[m2]\t XF[m]\n');

            for k = 1:numel(COND)
                c = COND(k);
                if ~isfinite(c.V)
                    continue;
                end

                hs = obj.computeDerivedHydro(c, vcg);

                fprintf(fid, ...
                    '%6.3f\t%10.3f\t%8.4f\t%8.4f\t%8.4f\t%12.3f\t%12.3f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.3f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\n', ...
                    c.T, c.V, c.B(1), c.B(2), c.B(3), c.I(1), c.I(2), ...
                    hs.BMt, hs.BMl, hs.KMt, hs.GMt, hs.GMl, hs.MCT, ...
                    c.Af, c.Lwl, c.Bwl, c.Ams, c.Ws, c.F(1));
            end

            fprintf('  Written: %s\n', datFile);
        end

        function writeHydroTxt(obj, outputTag, COND, vcg)
            txtFile = fullfile('Output', ['Hydrostatics_', outputTag, '.txt']);
            fid = fopen(txtFile, 'w');
            if fid == -1
                warning('PostProcessor:FileOpen', 'Could not open %s for writing.', txtFile);
                return;
            end
            cleaner = onCleanup(@() fclose(fid));

            sep = repmat('=', 1, 243);
            fprintf(fid, '%s\n', sep);
            fprintf(fid, '  NAUTILUS - Hydrostatic Results (Mesh)\n');
            fprintf(fid, '  Hull: %s    |   Drafts: %d   |   Density: %.0f kg/m3\n', outputTag, numel(COND), obj.cfg.rho);
            fprintf(fid, '%s\n\n', sep);

            hdr = {'T[m]','V[m3]','Delta[t]','XB[m]','YB[m]','KB[m]', ...
                'Ix[m4]','Iy[m4]','BMt[m]','BMl[m]','KMt[m]','GMt[m]','GMl[m]','MCT[tm/cm]', ...
                'Af[m2]','TPC[t/cm]','XF[m]','Lwl[m]','Bwl[m]','Ws[m2]','Cb','Cw','Cm','Cp'};
            fmt_h = '%-8s  %-10s  %-9s  %-8s  %-8s  %-8s  %-12s  %-12s  %-8s  %-9s  %-8s  %-8s  %-9s  %-10s  %-8s  %-9s  %-8s  %-8s  %-8s  %-8s  %-6s  %-6s  %-6s  %-6s\n';
            fprintf(fid, fmt_h, hdr{:});
            fprintf(fid, '%s\n', repmat('-', 1, 243));

            fmt_d = '%-8.3f  %-10.3f  %-9.2f  %-8.4f  %-8.4f  %-8.4f  %-12.3f  %-12.3f  %-8.4f  %-9.4f  %-8.4f  %-8.4f  %-9.4f  %-10.4f  %-8.3f  %-9.4f  %-8.4f  %-8.4f  %-8.4f  %-8.3f  %-6.4f  %-6.4f  %-6.4f  %-6.4f\n';

            for k = 1:numel(COND)
                c = COND(k);
                if ~isfinite(c.V)
                    continue;
                end

                hs = obj.computeDerivedHydro(c, vcg);

                fprintf(fid, fmt_d, ...
                    c.T, c.V, hs.Delta, c.B(1), c.B(2), c.B(3), c.I(1), c.I(2), ...
                    hs.BMt, hs.BMl, hs.KMt, hs.GMt, hs.GMl, hs.MCT, ...
                    c.Af, hs.TPC, c.F(1), c.Lwl, c.Bwl, c.Ws, hs.Cb, hs.Cw, hs.Cm, hs.Cp);
            end

            fprintf('  Written: %s\n', txtFile);
        end

        function writeGZDat(~, outputTag, GZdata)
            datFile = fullfile('Output', ['GZ_', outputTag, '.dat']);
            fid = fopen(datFile, 'w');
            if fid == -1
                warning('PostProcessor:FileOpen', 'Could not open %s for writing.', datFile);
                return;
            end
            cleaner = onCleanup(@() fclose(fid));

            fprintf(fid, '# heel[deg]\t trim[deg]\t draft[m]\t GZ[m]\t V[m3]\t converged\n');
            for k = 1:numel(GZdata)
                g = GZdata(k);
                fprintf(fid, '%8.3f\t%8.4f\t%8.4f\t%8.5f\t%12.4f\t%d\n', ...
                    g.heel, g.trim, g.draft, g.GZ, g.V, g.conv);
            end

            fprintf('  Written: %s\n', datFile);
        end

        function writeGZTxt(~, outputTag, GZdata)
            txtFile = fullfile('Output', ['GZ_', outputTag, '.txt']);
            fid = fopen(txtFile, 'w');
            if fid == -1
                warning('PostProcessor:FileOpen', 'Could not open %s for writing.', txtFile);
                return;
            end
            cleaner = onCleanup(@() fclose(fid));

            fprintf(fid, 'NAUTILUS - GZ Curve Results\n');
            fprintf(fid, '%s\n', repmat('=', 1, 70));
            fprintf(fid, '%-10s %-10s %-10s %-12s %-12s %-10s\n', ...
                'Heel[deg]', 'Trim[deg]', 'Draft[m]', 'GZ[m]', 'Volume[m3]', 'Status');
            fprintf(fid, '%s\n', repmat('-', 1, 70));

            for k = 1:numel(GZdata)
                g = GZdata(k);
                status = 'FAIL';
                if g.conv
                    status = 'OK';
                end
                fprintf(fid, '%-10.2f %-10.4f %-10.4f %-12.5f %-12.4f %-10s\n', ...
                    g.heel, g.trim, g.draft, g.GZ, g.V, status);
            end

            heelAll = [GZdata.heel];
            gzAll   = [GZdata.GZ];
            valid   = isfinite(heelAll) & isfinite(gzAll);

            if any(valid)
                heelV = heelAll(valid);
                gzV   = gzAll(valid);

                [gzMax, iMax] = max(gzV);
                heelMax = heelV(iMax);

                iZero = find(gzV < 0, 1, 'first');
                if ~isempty(iZero) && iZero > 1
                    heelZero = interp1(gzV(iZero-1:iZero), heelV(iZero-1:iZero), 0, 'linear', 'extrap');
                else
                    heelZero = NaN;
                end

                gzArea = trapz(deg2rad(heelV), gzV);

                fprintf(fid, '\n%s\n', repmat('-', 1, 70));
                fprintf(fid, 'Stability Summary\n');
                fprintf(fid, '  Max GZ:   %.4f m at %.1f deg\n', gzMax, heelMax);
                if ~isnan(heelZero)
                    fprintf(fid, '  Range:    0 to %.1f deg\n', heelZero);
                end
                fprintf(fid, '  GZ area:  %.4f m*rad\n', gzArea);
            end

            fprintf('  Written: %s\n', txtFile);
        end

        function plotGZ(obj, outputTag, GZdata)
            heel = [GZdata.heel];
            gz   = [GZdata.GZ];
            ok   = isfinite(heel) & isfinite(gz);

            if ~any(ok)
                return;
            end

            heel = heel(ok);
            gz   = gz(ok);

            C = obj.publicationColors();
            fig = obj.createPublicationFigure([160 120 1020 650]);
            ax = axes('Parent', fig);
            hold(ax, 'on');

            obj.drawZeroLine(ax, C);

            plot(ax, heel, gz, '-o', ...
                'LineWidth', 2.8, ...
                'Color', C.blue, ...
                'MarkerSize', 6.8, ...
                'MarkerFaceColor', [1 1 1], ...
                'MarkerEdgeColor', C.blue);

            % light area fill under positive GZ
            gzFill = gz;
            gzFill(gzFill < 0) = 0;
            patch(ax, [heel, fliplr(heel)], [gzFill, zeros(size(gzFill))], C.blue, ...
                'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');

            [gzPeak, iPeak] = max(gz);
            heelPeak = heel(iPeak);

            plot(ax, heelPeak, gzPeak, 'p', ...
                'MarkerSize', 13, ...
                'LineWidth', 1.3, ...
                'MarkerFaceColor', C.gold, ...
                'MarkerEdgeColor', C.dark, ...
                'DisplayName', '$GZ_{\max}$');

            text(ax, heelPeak, gzPeak, ...
                sprintf('$\\;GZ_{\\max}=%.3f\\,\\mathrm{m}\\ @\\ %.1f^\\circ$', gzPeak, heelPeak), ...
                'Interpreter', 'latex', ...
                'FontSize', 12, ...
                'Color', C.dark, ...
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'left');

            obj.styleAxes(ax);
            xlim(ax, [min(heel), max(heel)]);

            % Better labels
            xlabel(ax, '$\phi\,[^\circ]$', 'Interpreter', 'latex');
            ylabel(ax, '$GZ\;[\mathrm{m}]$', 'Interpreter', 'latex');

            obj.applyAxesTitle(ax, 'Righting Lever Curve');
            legend(ax, {'$GZ(\phi)$', '$GZ_{\max}$'}, ...
                'Interpreter', 'latex', ...
                'Location', 'northwest', ...
                'Box', 'off');

            if obj.cfg.savePlots
                obj.exportFigure(fig, fullfile('Output', ['GZ_', outputTag]));
            end
        end

        function plotHydrostaticsBasic(obj, outputTag, COND, vcg)
            % Compute all hydrostatic arrays once
            HS = obj.extractHydroArrays(COND, vcg);

            valid = isfinite(HS.T) & isfinite(HS.V) & HS.V > 0;
            if ~any(valid)
                warning('PostProcessor:NoValidData', 'No valid hydrostatic data available for plotting.');
                return;
            end

            C = obj.publicationColors();
            fig = obj.createPublicationFigure([90 70 1320 860]);
            tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

            % -----------------------------------------------------------------
            % 1) Displacement
            % -----------------------------------------------------------------
            ax = nexttile(tl, 1);
            hold(ax, 'on');
            plot(ax, HS.T(valid), HS.V(valid), '-o', ...
                'LineWidth', 2.5, ...
                'Color', C.blue, ...
                'MarkerSize', 6.0, ...
                'MarkerFaceColor', 'white', ...
                'MarkerEdgeColor', C.blue);
            obj.fillUnderCurve(ax, HS.T(valid), HS.V(valid), C.blue, 0.08);
            obj.styleAxes(ax);
            xlabel(ax, '$T\;[\mathrm{m}]$', 'Interpreter', 'latex');
            ylabel(ax, '$\nabla\;[\mathrm{m}^3]$', 'Interpreter', 'latex');
            obj.applyAxesTitle(ax, 'Displacement Volume');

            % -----------------------------------------------------------------
            % 2) Stability heights
            % -----------------------------------------------------------------
            ax = nexttile(tl, 2);
            validStab = valid & isfinite(HS.KB) & isfinite(HS.KMt) & isfinite(HS.GMt);
            hold(ax, 'on');
            if any(validStab)
                plot(ax, HS.T(validStab), HS.KB(validStab), '-o', ...
                    'LineWidth', 2.1, 'Color', C.green, ...
                    'MarkerSize', 5.5, 'MarkerFaceColor', 'white', ...
                    'DisplayName', '$KB$');

                plot(ax, HS.T(validStab), HS.KMt(validStab), '-s', ...
                    'LineWidth', 2.1, 'Color', C.orange, ...
                    'MarkerSize', 5.5, 'MarkerFaceColor', 'white', ...
                    'DisplayName', '$KM_t$');

                plot(ax, HS.T(validStab), HS.GMt(validStab), '-^', ...
                    'LineWidth', 2.1, 'Color', C.purple, ...
                    'MarkerSize', 5.5, 'MarkerFaceColor', 'white', ...
                    'DisplayName', '$GM_t$');

                obj.drawZeroLine(ax, C);
                legend(ax, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');
            end
            obj.styleAxes(ax);
            xlabel(ax, '$T\;[\mathrm{m}]$', 'Interpreter', 'latex');
            ylabel(ax, '$\mathrm{Height}\;[\mathrm{m}]$', 'Interpreter', 'latex');
            obj.applyAxesTitle(ax, 'Transverse Stability Heights');

            % -----------------------------------------------------------------
            % 3) Waterplane area and wetted surface
            % -----------------------------------------------------------------
            ax = nexttile(tl, 3);
            validArea = valid & isfinite(HS.Af) & isfinite(HS.Ws);
            hold(ax, 'on');
            if any(validArea)
                plot(ax, HS.T(validArea), HS.Af(validArea), '-o', ...
                    'LineWidth', 2.2, ...
                    'Color', C.teal, ...
                    'MarkerSize', 5.8, ...
                    'MarkerFaceColor', 'white', ...
                    'DisplayName', '$A_{\mathrm{WP}}$');

                plot(ax, HS.T(validArea), HS.Ws(validArea), '-s', ...
                    'LineWidth', 2.2, ...
                    'Color', C.dark, ...
                    'MarkerSize', 5.8, ...
                    'MarkerFaceColor', 'white', ...
                    'DisplayName', '$S_{\mathrm{wet}}$');

                legend(ax, 'Interpreter', 'latex', 'Location', 'northwest', 'Box', 'off');
            end
            obj.styleAxes(ax);
            xlabel(ax, '$T\;[\mathrm{m}]$', 'Interpreter', 'latex');
            ylabel(ax, '$\mathrm{Area}\;[\mathrm{m}^2]$', 'Interpreter', 'latex');
            obj.applyAxesTitle(ax, 'Waterplane and Wetted Areas');

            % -----------------------------------------------------------------
            % 4) Waterline dimensions
            % -----------------------------------------------------------------
            ax = nexttile(tl, 4);
            validDims = valid & isfinite(HS.Lwl) & isfinite(HS.Bwl);
            hold(ax, 'on');
            if any(validDims)
                plot(ax, HS.T(validDims), HS.Lwl(validDims), '-o', ...
                    'LineWidth', 2.2, ...
                    'Color', C.gold, ...
                    'MarkerSize', 5.8, ...
                    'MarkerFaceColor', 'white', ...
                    'DisplayName', '$L_{\mathrm{WL}}$');

                plot(ax, HS.T(validDims), HS.Bwl(validDims), '-s', ...
                    'LineWidth', 2.2, ...
                    'Color', C.blue2, ...
                    'MarkerSize', 5.8, ...
                    'MarkerFaceColor', 'white', ...
                    'DisplayName', '$B_{\mathrm{WL}}$');

                legend(ax, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');
            end
            obj.styleAxes(ax);
            xlabel(ax, '$T\;[\mathrm{m}]$', 'Interpreter', 'latex');
            ylabel(ax, '$\mathrm{Length}\;[\mathrm{m}]$', 'Interpreter', 'latex');
            obj.applyAxesTitle(ax, 'Waterline Dimensions');

            % -----------------------------------------------------------------
            % 5) Hydrostatic coefficients
            % -----------------------------------------------------------------
            ax = nexttile(tl, 5);
            validCoeff = valid & isfinite(HS.Cb) & isfinite(HS.Cw) & isfinite(HS.Cm) & isfinite(HS.Cp);
            hold(ax, 'on');
            if any(validCoeff)
                plot(ax, HS.T(validCoeff), HS.Cb(validCoeff), '-o', ...
                    'LineWidth', 2.0, 'Color', C.blue, ...
                    'MarkerSize', 5.2, 'MarkerFaceColor', 'white', ...
                    'DisplayName', '$C_B$');

                plot(ax, HS.T(validCoeff), HS.Cw(validCoeff), '-s', ...
                    'LineWidth', 2.0, 'Color', C.green, ...
                    'MarkerSize', 5.2, 'MarkerFaceColor', 'white', ...
                    'DisplayName', '$C_W$');

                plot(ax, HS.T(validCoeff), HS.Cm(validCoeff), '-^', ...
                    'LineWidth', 2.0, 'Color', C.orange, ...
                    'MarkerSize', 5.2, 'MarkerFaceColor', 'white', ...
                    'DisplayName', '$C_M$');

                plot(ax, HS.T(validCoeff), HS.Cp(validCoeff), '-d', ...
                    'LineWidth', 2.0, 'Color', C.purple, ...
                    'MarkerSize', 5.2, 'MarkerFaceColor', 'white', ...
                    'DisplayName', '$C_P$');

                ylim(ax, [0, max(1.05, 1.05*max([HS.Cb(validCoeff), HS.Cw(validCoeff), HS.Cm(validCoeff), HS.Cp(validCoeff)]))]);
                legend(ax, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');
            end
            obj.styleAxes(ax);
            xlabel(ax, '$T\;[\mathrm{m}]$', 'Interpreter', 'latex');
            ylabel(ax, '$\mathrm{Coefficient}\;[-]$', 'Interpreter', 'latex');
            obj.applyAxesTitle(ax, 'Hydrostatic Coefficients');

            % -----------------------------------------------------------------
            % 6) TPC and MCT
            % -----------------------------------------------------------------
            ax = nexttile(tl, 6);
            validCap = valid & isfinite(HS.TPC) & isfinite(HS.MCT);
            hold(ax, 'on');
            if any(validCap)
                yyaxis(ax, 'left');
                plot(ax, HS.T(validCap), HS.TPC(validCap), '-o', ...
                    'LineWidth', 2.2, ...
                    'Color', C.teal, ...
                    'MarkerSize', 5.8, ...
                    'MarkerFaceColor', 'white');
                ylabel(ax, '$TPC\;[\mathrm{t/cm}]$', 'Interpreter', 'latex');

                yyaxis(ax, 'right');
                plot(ax, HS.T(validCap), HS.MCT(validCap), '-s', ...
                    'LineWidth', 2.2, ...
                    'Color', C.red, ...
                    'MarkerSize', 5.8, ...
                    'MarkerFaceColor', 'white');
                ylabel(ax, '$MCT_{1\mathrm{cm}}\;[\mathrm{t\,m/cm}]$', 'Interpreter', 'latex');

                yyaxis(ax, 'left');
            end
            obj.styleAxes(ax);
            xlabel(ax, '$T\;[\mathrm{m}]$', 'Interpreter', 'latex');
            obj.applyAxesTitle(ax, 'Hydrostatic Stiffness Indicators');

            obj.applyLayoutTitle(tl, 'Hydrostatic Curves');

            if obj.cfg.savePlots
                obj.exportFigure(fig, fullfile('Output', ['Hydrostatics_', outputTag]));
            end
        end

        function HS = extractHydroArrays(obj, COND, vcg)
            n = numel(COND);

            HS = struct();
            HS.T   = nan(1,n);
            HS.V   = nan(1,n);
            HS.KB  = nan(1,n);
            HS.BMt = nan(1,n);
            HS.BMl = nan(1,n);
            HS.KMt = nan(1,n);
            HS.GMt = nan(1,n);
            HS.GMl = nan(1,n);
            HS.Af  = nan(1,n);
            HS.Lwl = nan(1,n);
            HS.Bwl = nan(1,n);
            HS.Ams = nan(1,n);
            HS.Ws  = nan(1,n);
            HS.TPC = nan(1,n);
            HS.MCT = nan(1,n);
            HS.Cb  = nan(1,n);
            HS.Cw  = nan(1,n);
            HS.Cm  = nan(1,n);
            HS.Cp  = nan(1,n);

            for i = 1:n
                c = COND(i);

                HS.T(i)   = c.T;
                HS.V(i)   = c.V;
                HS.KB(i)  = c.B(3);
                HS.Af(i)  = c.Af;
                HS.Lwl(i) = c.Lwl;
                HS.Bwl(i) = c.Bwl;
                HS.Ams(i) = c.Ams;
                HS.Ws(i)  = c.Ws;

                if ~isfinite(c.V) || c.V <= 0
                    continue;
                end

                hs = obj.computeDerivedHydro(c, vcg);
                HS.BMt(i) = hs.BMt;
                HS.BMl(i) = hs.BMl;
                HS.KMt(i) = hs.KMt;
                HS.GMt(i) = hs.GMt;
                HS.GMl(i) = hs.GMl;
                HS.TPC(i) = hs.TPC;
                HS.MCT(i) = hs.MCT;
                HS.Cb(i)  = hs.Cb;
                HS.Cw(i)  = hs.Cw;
                HS.Cm(i)  = hs.Cm;
                HS.Cp(i)  = hs.Cp;
            end
        end

        function hs = computeDerivedHydro(obj, c, vcg)
            tiny = 1e-9;

            V = c.V;
            KB = c.B(3);
            BMt = c.I(1) / max(V, tiny);
            BMl = c.I(2) / max(V, tiny);
            KMt = KB + BMt;
            GMt = KMt - vcg;
            GMl = KB + BMl - vcg;
            Delta = obj.cfg.rho * V / 1000;
            TPC = c.Af * obj.cfg.rho / (1000 * 100);
            MCT = Delta * GMl / max(100 * c.Lwl, tiny);

            denomCb = max(c.Lwl * c.Bwl * c.T, tiny);
            denomCw = max(c.Lwl * c.Bwl, tiny);
            denomCm = max(c.Bwl * c.T, tiny);
            denomCp = max(c.Ams * c.Lwl, tiny);

            Cb = obj.clamp01(V / denomCb);
            Cw = obj.clamp01(c.Af / denomCw);
            Cm = obj.clamp01(c.Ams / denomCm);
            Cp = obj.clamp01(V / denomCp);

            hs = struct();
            hs.V = V;
            hs.KB = KB;
            hs.BMt = BMt;
            hs.BMl = BMl;
            hs.KMt = KMt;
            hs.GMt = GMt;
            hs.GMl = GMl;
            hs.Delta = Delta;
            hs.TPC = TPC;
            hs.MCT = MCT;
            hs.Cb = Cb;
            hs.Cw = Cw;
            hs.Cm = Cm;
            hs.Cp = Cp;
        end

        function y = clamp01(~, x)
            y = max(0, min(1.01, x));
        end

        function drawZeroLine(~, ax, C)
            xl = xlim(ax);
            plot(ax, xl, [0 0], '--', ...
                'Color', C.red, ...
                'LineWidth', 1.1, ...
                'HandleVisibility', 'off');
        end

        function fillUnderCurve(~, ax, x, y, color, alphaVal)
            if numel(x) < 2 || numel(y) < 2
                return;
            end
            patch(ax, [x, fliplr(x)], [y, zeros(size(y))], color, ...
                'FaceAlpha', alphaVal, ...
                'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end

        function hideAxesToolbars(~, fig)
            axs = findall(fig, 'Type', 'axes');
            for i = 1:numel(axs)
                ax = axs(i);
                if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar) && isprop(ax.Toolbar, 'Visible')
                    ax.Toolbar.Visible = 'off';
                end
            end
        end

        function fig = createPublicationFigure(obj, position)
            if nargin < 2
                position = [120 120 1000 700];
            end
            fig = figure( ...
                'Color', 'w', ...
                'Position', position, ...
                'Renderer', 'painters', ...
                'InvertHardcopy', 'off');
            obj.hideAxesToolbars(fig);
        end

        function styleAxes(~, ax)
            ax.FontName = 'Times New Roman';
            ax.FontSize = 12.5;
            ax.LineWidth = 1.05;
            ax.TickDir = 'out';
            ax.TickLength = [0.014 0.014];
            ax.Box = 'on';
            ax.Layer = 'top';
            ax.GridAlpha = 0.14;
            ax.MinorGridAlpha = 0.08;
            ax.XMinorGrid = 'on';
            ax.YMinorGrid = 'on';
            ax.XGrid = 'on';
            ax.YGrid = 'on';
            ax.TitleFontWeight = 'normal';
        end

        function exportFigure(obj, fig, basePath)
            drawnow;
            obj.hideAxesToolbars(fig);
            set(fig, 'PaperPositionMode', 'auto');

            print(fig, [basePath, '.png'], '-dpng', '-r350');
            fprintf('  Saved:   %s.png\n', basePath);

            try
                print(fig, [basePath, '.pdf'], '-dpdf', '-painters');
                fprintf('  Saved:   %s.pdf\n', basePath);
            catch ME
                warning('PostProcessor:PDFExportFailed', ...
                    'Could not export PDF for %s: %s', basePath, ME.message);
            end
        end

        function applyAxesTitle(obj, target, heading)
            title(target, obj.latexHeading(heading), ...
                'Interpreter', 'latex', ...
                'FontSize', 14);
        end

        function applyLayoutTitle(obj, target, heading)
            title(target, obj.latexHeading(heading), ...
                'Interpreter', 'latex', ...
                'FontSize', 18);
        end

        function txt = latexHeading(obj, txt)
            txt = ['$\mathrm{' obj.latexText(txt) '}$'];
        end

        function txt = latexText(obj, txt)
            txt = obj.escapeLatex(txt);
            txt = strrep(txt, '-', '{-}');
            txt = strrep(txt, ' ', '\ ');
        end

        function txt = escapeLatex(~, txt)
            txt = strrep(txt, '\', '\textbackslash ');
            txt = strrep(txt, '_', '\_');
            txt = strrep(txt, '%', '\%');
            txt = strrep(txt, '#', '\#');
            txt = strrep(txt, '&', '\&');
            txt = strrep(txt, '{', '\{');
            txt = strrep(txt, '}', '\}');
            txt = strrep(txt, '^', '\^{}');
        end

        function C = publicationColors(~)
            C = struct();
            C.blue   = [0.07 0.32 0.62];
            C.blue2  = [0.23 0.53 0.82];
            C.green  = [0.13 0.55 0.23];
            C.orange = [0.82 0.39 0.06];
            C.purple = [0.43 0.21 0.62];
            C.teal   = [0.02 0.52 0.58];
            C.gold   = [0.78 0.58 0.12];
            C.red    = [0.70 0.16 0.16];
            C.dark   = [0.18 0.18 0.18];
        end
    end
end