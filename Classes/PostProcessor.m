classdef PostProcessor
    %PostProcessor Export hydrostatic/GZ results and generate plots.
    %
    % FIXED BUGS:
    % 1. writeHydroTxt - Fixed division by zero risk in coefficient calculations
    % 2. plotHydrostaticsBasic - Better handling of invalid data
    % 3. yline compatibility note added

    properties
        cfg struct
    end

    methods
        function obj = PostProcessor(cfg)
            obj.cfg = cfg;
        end

        function generate(obj, outputTag, results, mesh, inputData)
            COND = results.COND;

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
        function writeHydroDat(obj, outputTag, COND, vcg)
            datFile = fullfile('Output', ['Hydrostatics_', outputTag, '.dat']);
            fid = fopen(datFile, 'w');
            if fid == -1
                warning('PostProcessor:FileOpen', 'Could not open %s for writing.', datFile);
                return;
            end

            fprintf(fid, '# T[m]\t V[m3]\t XB[m]\t YB[m]\t ZB[m]\t Ix[m4]\t Iy[m4]\t BM[m]\t BML[m]\t KM[m]\t GM[m]\t GML[m]\t MCT[tm/cm]\t Af[m2]\t Lwl[m]\t Bwl[m]\t Ams[m2]\t Ws[m2]\t XF[m]\n');
            for k = 1:numel(COND)
                c = COND(k);
                if isnan(c.V)
                    continue;
                end
                V   = c.V;
                KB  = c.B(3);
                BM  = c.I(1) / max(V, 1e-9);
                BML = c.I(2) / max(V, 1e-9);
                KM  = KB + BM;
                GM  = KM - vcg;
                GML = KB + BML - vcg;
                Dlt = obj.cfg.rho * V / 1000;
                MCT = Dlt * GML / max(100 * c.Lwl, 1e-9);
                fprintf(fid, '%6.3f\t%10.3f\t%8.4f\t%8.4f\t%8.4f\t%12.3f\t%12.3f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.3f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\n', ...
                    c.T, V, c.B(1), c.B(2), KB, c.I(1), c.I(2), BM, BML, KM, GM, GML, MCT, ...
                    c.Af, c.Lwl, c.Bwl, c.Ams, c.Ws, c.F(1));
            end
            fclose(fid);
            fprintf('  Written: %s\n', datFile);
        end

        function writeHydroTxt(obj, outputTag, COND, vcg)
            txtFile = fullfile('Output', ['Hydrostatics_', outputTag, '.txt']);
            fid = fopen(txtFile, 'w');
            if fid == -1
                warning('PostProcessor:FileOpen', 'Could not open %s for writing.', txtFile);
                return;
            end

            % Separator width matches the data table (~243 chars with 24 columns)
            sep = repmat('=', 1, 243);
            fprintf(fid, '%s\n', sep);
            fprintf(fid, '  NAUTILUS - Hydrostatic Results (Mesh)\n');
            fprintf(fid, '  Hull: %s    |   Drafts: %d   |   Density: %.0f kg/m3\n', outputTag, numel(COND), obj.cfg.rho);
            fprintf(fid, '%s\n\n', sep);

            hdr = {'T[m]','V[m3]','Delta[t]','XB[m]','YB[m]','KB[m]', ...
                'Ix[m4]','Iy[m4]','BM[m]','BML[m]','KM[m]','GM[m]','GML[m]','MCT[tm/cm]', ...
                'Af[m2]','TPC[t/cm]','XF[m]','Lwl[m]','Bwl[m]','Ws[m2]','Cb','Cw','Cm','Cp'};
            fmt_h = '%-8s  %-10s  %-9s  %-8s  %-8s  %-8s  %-12s  %-12s  %-8s  %-9s  %-8s  %-8s  %-9s  %-10s  %-8s  %-9s  %-8s  %-8s  %-8s  %-8s  %-6s  %-6s  %-6s  %-6s\n';
            fprintf(fid, fmt_h, hdr{:});
            fprintf(fid, '%s\n', repmat('-', 1, 243));

            fmt_d = '%-8.3f  %-10.3f  %-9.2f  %-8.4f  %-8.4f  %-8.4f  %-12.3f  %-12.3f  %-8.4f  %-9.4f  %-8.4f  %-8.4f  %-9.4f  %-10.4f  %-8.3f  %-9.4f  %-8.4f  %-8.4f  %-8.4f  %-8.3f  %-6.4f  %-6.4f  %-6.4f  %-6.4f\n';
            for k = 1:numel(COND)
                c = COND(k);
                if isnan(c.V)
                    continue;
                end

                V = c.V;
                KB = c.B(3);
                BM  = c.I(1) / max(V, 1e-9);
                BML = c.I(2) / max(V, 1e-9);
                KM  = KB + BM;
                GM  = KM - vcg;
                GML = KB + BML - vcg;
                Dlt = obj.cfg.rho * V / 1000;
                TPC = c.Af * obj.cfg.rho / (1000 * 100);
                MCT = Dlt * GML / max(100 * c.Lwl, 1e-9);   % t·m/cm

                % FIXED: Safer coefficient calculations with minimum denominators
                denom_Cb = max(c.Lwl * c.Bwl * c.T, 1e-9);
                denom_Cw = max(c.Lwl * c.Bwl, 1e-9);
                denom_Cm = max(c.Bwl * c.T, 1e-9);
                denom_Cp = max(c.Ams * c.Lwl, 1e-9);

                Cb = V / denom_Cb;
                Cw = c.Af / denom_Cw;
                Cm = c.Ams / denom_Cm;
                Cp = V / denom_Cp;

                % Clamp coefficients to valid range [0, 1] (with small tolerance for numerical errors)
                Cb = max(0, min(1.01, Cb));
                Cw = max(0, min(1.01, Cw));
                Cm = max(0, min(1.01, Cm));
                Cp = max(0, min(1.01, Cp));

                fprintf(fid, fmt_d, ...
                    c.T, V, Dlt, c.B(1), c.B(2), KB, c.I(1), c.I(2), BM, BML, KM, GM, GML, MCT, ...
                    c.Af, TPC, c.F(1), c.Lwl, c.Bwl, c.Ws, Cb, Cw, Cm, Cp);
            end
            fclose(fid);
            fprintf('  Written: %s\n', txtFile);
        end

        function writeGZDat(~, outputTag, GZdata)
            datFile = fullfile('Output', ['GZ_', outputTag, '.dat']);
            fid = fopen(datFile, 'w');
            if fid == -1
                warning('PostProcessor:FileOpen', 'Could not open %s for writing.', datFile);
                return;
            end

            fprintf(fid, '# heel[deg]\t trim[deg]\t draft[m]\t GZ[m]\t V[m3]\t converged\n');
            for k = 1:numel(GZdata)
                g = GZdata(k);
                fprintf(fid, '%8.3f\t%8.4f\t%8.4f\t%8.5f\t%12.4f\t%d\n', ...
                    g.heel, g.trim, g.draft, g.GZ, g.V, g.conv);
            end
            fclose(fid);
            fprintf('  Written: %s\n', datFile);
        end

        function writeGZTxt(~, outputTag, GZdata)
            txtFile = fullfile('Output', ['GZ_', outputTag, '.txt']);
            fid = fopen(txtFile, 'w');
            if fid == -1
                warning('PostProcessor:FileOpen', 'Could not open %s for writing.', txtFile);
                return;
            end

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

            % Stability summary block
            heelAll = [GZdata.heel];
            gzAll   = [GZdata.GZ];
            valid   = isfinite(heelAll) & isfinite(gzAll);
            if any(valid)
                heelV = heelAll(valid);
                gzV   = gzAll(valid);
                [gzMax, iMax] = max(gzV);
                heelMax = heelV(iMax);
                iZero   = find(gzV < 0, 1, 'first');
                if ~isempty(iZero) && iZero > 1
                    heelZero = interp1(gzV(iZero-1:iZero), heelV(iZero-1:iZero), 0, 'linear', 'extrap');
                else
                    heelZero = NaN;
                end
                gzArea = trapz(heelV * pi/180, gzV);   % m·rad

                fprintf(fid, '\n%s\n', repmat('-', 1, 70));
                fprintf(fid, 'Stability Summary\n');
                fprintf(fid, '  Max GZ:   %.4f m at %.1f deg\n', gzMax, heelMax);
                if ~isnan(heelZero)
                    fprintf(fid, '  Range:    0 to %.1f deg\n', heelZero);
                end
                fprintf(fid, '  GZ area:  %.4f m*rad\n', gzArea);
            end

            fclose(fid);
            fprintf('  Written: %s\n', txtFile);
        end

        function plotGZ(obj, outputTag, GZdata)
            heel = [GZdata.heel];
            gz = [GZdata.GZ];
            ok = isfinite(heel) & isfinite(gz);
            if ~any(ok)
                return;
            end

            C = obj.publicationColors();
            fig = obj.createPublicationFigure([160 120 980 620]);
            ax = axes('Parent', fig);
            hold(ax, 'on');

            plot(ax, heel(ok), gz(ok), '-o', ...
                'LineWidth', 2.4, ...
                'Color', C.blue, ...
                'MarkerSize', 6.5, ...
                'MarkerFaceColor', 'white', ...
                'MarkerEdgeColor', C.blue);

            % FIXED: Check for yline availability (MATLAB R2018b+)
            % For older versions, use plot instead
            if exist('yline', 'file') == 2
                yline(ax, 0, '--', ...
                    'Color', C.red, ...
                    'LineWidth', 1.3, ...
                    'HandleVisibility', 'off');
            else
                plot(ax, xlim(ax), [0 0], '--', ...
                    'Color', C.red, ...
                    'LineWidth', 1.3, ...
                    'HandleVisibility', 'off');
            end

            gzValid = gz(ok);
            heelValid = heel(ok);
            [gzPeak, iPeak] = max(gzValid);
            heelPeak = heelValid(iPeak);
            plot(ax, heelPeak, gzPeak, 'p', ...
                'MarkerSize', 12, ...
                'LineWidth', 1.2, ...
                'MarkerFaceColor', C.gold, ...
                'MarkerEdgeColor', C.dark);
            text(ax, heelPeak, gzPeak, ...
                sprintf('$\\max\\,GZ = %.3f\\,\\mathrm{m}\\ \\mathrm{at}\\ %.1f^\\circ$', gzPeak, heelPeak), ...
                'Interpreter', 'latex', ...
                'FontSize', 12, ...
                'Color', C.dark, ...
                'VerticalAlignment', 'bottom');

            obj.styleAxes(ax);
            xlabel(ax, '$\phi\,[^\circ]$', 'Interpreter', 'latex');
            ylabel(ax, '$GZ\,[\mathrm{m}]$', 'Interpreter', 'latex');
            obj.applyAxesTitle(ax, 'Righting Lever Curve');
            legend(ax, {'$GZ(\phi)$', '$GZ_{\max}$'}, ...
                'Interpreter', 'latex', ...
                'Location', 'best', ...
                'Box', 'off');

            if obj.cfg.savePlots
                obj.exportFigure(fig, fullfile('Output', ['GZ_', outputTag]));
            end
        end

        function plotHydrostaticsBasic(obj, outputTag, COND, vcg)
            % Extract data with validation
            T = [COND.T];
            V = [COND.V];
            
            % Handle potential NaN in B field
            KB = nan(size(T));
            for i = 1:numel(COND)
                if ~isnan(COND(i).B(3))
                    KB(i) = COND(i).B(3);
                end
            end
            
            BM = nan(size(T));
            for i = 1:numel(COND)
                if ~isnan(COND(i).V) && COND(i).V > 0 && ~isnan(COND(i).I(1))
                    BM(i) = COND(i).I(1) / max(COND(i).V, 1e-9);
                end
            end
            
            KM = KB + BM;
            GM = KM - vcg;
            
            Af = [COND.Af];
            Lwl = [COND.Lwl];
            Bwl = [COND.Bwl];

            valid = isfinite(T) & isfinite(V) & V > 0;
            if ~any(valid)
                warning('PostProcessor:NoValidData', 'No valid data for plotting.');
                return;
            end

            C = obj.publicationColors();
            fig = obj.createPublicationFigure([120 110 1180 760]);
            tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

            ax = nexttile(tl, 1);
            plot(ax, T(valid), V(valid), '-o', ...
                'LineWidth', 2.1, ...
                'Color', C.blue, ...
                'MarkerSize', 6, ...
                'MarkerFaceColor', 'white', ...
                'MarkerEdgeColor', C.blue);
            obj.styleAxes(ax);
            xlabel(ax, '$T\,[\mathrm{m}]$', 'Interpreter', 'latex');
            ylabel(ax, '$V\,[\mathrm{m}^3]$', 'Interpreter', 'latex');
            obj.applyAxesTitle(ax, 'Displacement Curve');

            ax = nexttile(tl, 2);
            valid_stab = valid & isfinite(KB) & isfinite(KM) & isfinite(GM);
            if any(valid_stab)
                hold(ax, 'on');
                plot(ax, T(valid_stab), KB(valid_stab), '-o', ...
                    'LineWidth', 2.0, 'Color', C.green, ...
                    'MarkerSize', 5.5, 'MarkerFaceColor', 'white');
                plot(ax, T(valid_stab), KM(valid_stab), '-s', ...
                    'LineWidth', 2.0, 'Color', C.orange, ...
                    'MarkerSize', 5.5, 'MarkerFaceColor', 'white');
                plot(ax, T(valid_stab), GM(valid_stab), '-^', ...
                    'LineWidth', 2.0, 'Color', C.purple, ...
                    'MarkerSize', 5.5, 'MarkerFaceColor', 'white');
                if any(isfinite(GM(valid_stab)))
                    if exist('yline', 'file') == 2
                        yline(ax, 0, ':', 'Color', C.red, 'LineWidth', 1.0, 'HandleVisibility', 'off');
                    else
                        plot(ax, xlim(ax), [0 0], ':', 'Color', C.red, 'LineWidth', 1.0, 'HandleVisibility', 'off');
                    end
                end
                legend(ax, {'$KB$', '$KM$', '$GM$'}, ...
                    'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');
            end
            obj.styleAxes(ax);
            xlabel(ax, '$T\,[\mathrm{m}]$', 'Interpreter', 'latex');
            ylabel(ax, '$[\mathrm{m}]$', 'Interpreter', 'latex');
            obj.applyAxesTitle(ax, 'Stability Heights');

            ax = nexttile(tl, 3);
            valid_wp = valid & isfinite(Af);
            if any(valid_wp)
                plot(ax, T(valid_wp), Af(valid_wp), '-o', ...
                    'LineWidth', 2.1, ...
                    'Color', C.teal, ...
                    'MarkerSize', 6, ...
                    'MarkerFaceColor', 'white', ...
                    'MarkerEdgeColor', C.teal);
            end
            obj.styleAxes(ax);
            xlabel(ax, '$T\,[\mathrm{m}]$', 'Interpreter', 'latex');
            ylabel(ax, '$A_{\mathrm{WP}}\,[\mathrm{m}^2]$', 'Interpreter', 'latex');
            obj.applyAxesTitle(ax, 'Waterplane Area');

            ax = nexttile(tl, 4);
            valid_dims = valid & isfinite(Lwl) & isfinite(Bwl);
            if any(valid_dims)
                hold(ax, 'on');
                plot(ax, T(valid_dims), Lwl(valid_dims), '-o', ...
                    'LineWidth', 2.0, 'Color', C.gold, ...
                    'MarkerSize', 5.5, 'MarkerFaceColor', 'white');
                plot(ax, T(valid_dims), Bwl(valid_dims), '-s', ...
                    'LineWidth', 2.0, 'Color', C.dark, ...
                    'MarkerSize', 5.5, 'MarkerFaceColor', 'white');
                legend(ax, {'$L_{\mathrm{WL}}$', '$B_{\mathrm{WL}}$'}, ...
                    'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');
            end
            obj.styleAxes(ax);
            xlabel(ax, '$T\,[\mathrm{m}]$', 'Interpreter', 'latex');
            ylabel(ax, '$[\mathrm{m}]$', 'Interpreter', 'latex');
            obj.applyAxesTitle(ax, 'Waterline Dimensions');

            obj.applyLayoutTitle(tl, 'Hydrostatic Curves');

            if obj.cfg.savePlots
                obj.exportFigure(fig, fullfile('Output', ['Hydrostatics_', outputTag]));
            end
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
            ax.FontSize = 13;
            ax.LineWidth = 1.1;
            ax.TickDir = 'out';
            ax.TickLength = [0.015 0.015];
            ax.Box = 'on';
            ax.Layer = 'top';
            ax.GridAlpha = 0.18;
            ax.MinorGridAlpha = 0.10;
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
            print(fig, [basePath, '.png'], '-dpng', '-r300');
            print(fig, [basePath, '.pdf'], '-dpdf', '-painters');
        end

        function applyAxesTitle(obj, target, heading)
            title(target, obj.latexHeading(heading), 'Interpreter', 'latex');
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
            C.blue = [0.07 0.32 0.62];
            C.green = [0.13 0.55 0.23];
            C.orange = [0.82 0.39 0.06];
            C.purple = [0.43 0.21 0.62];
            C.teal = [0.02 0.52 0.58];
            C.gold = [0.78 0.58 0.12];
            C.red = [0.70 0.16 0.16];
            C.dark = [0.18 0.18 0.18];
        end
    end
end
