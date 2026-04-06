classdef GeometryChecker
    %GeometryChecker Run mesh sanity checks for hydrostatic analysis.
    %
    % FIXED BUGS:
    % 1. Improved face index validation with tolerance
    % 2. Better degenerate face detection
    % 3. More informative error messages

    methods
        function [mesh, report] = check(~, mesh, options)
            if nargin < 3 || isempty(options)
                options = struct();
            end
            if ~isfield(options, 'areaTolerance'); options.areaTolerance = 1e-12; end
            if ~isfield(options, 'removeDegenerateFaces'); options.removeDegenerateFaces = true; end
            if ~isfield(options, 'requireWatertight'); options.requireWatertight = false; end

            report = struct();
            report.isValid = true;
            report.messages = {};

            if ~isstruct(mesh) || ~isfield(mesh, 'vertices') || ~isfield(mesh, 'faces')
                error('GeometryChecker:InvalidMesh', ...
                    'Mesh must be a struct with vertices and faces fields.');
            end

            V = mesh.vertices;
            F = mesh.faces;

            validateattributes(V, {'numeric'}, {'2d', 'ncols', 3, 'finite', 'real'});
            validateattributes(F, {'numeric'}, {'2d', 'ncols', 3, 'real', 'positive'});

            % FIXED: Use rounding with tolerance for face indices
            F = round(F);
            nV = size(V, 1);
            
            % Check for out-of-range indices
            if any(F(:) < 1) || any(F(:) > nV)
                invalidFaces = find(any(F < 1, 2) | any(F > nV, 2));
                error('GeometryChecker:FaceIndexOutOfRange', ...
                    'Mesh face indices exceed vertex array bounds. Found %d invalid faces (e.g., face %d).', ...
                    numel(invalidFaces), invalidFaces(1));
            end

            % Compute face areas
            a = V(F(:,1), :);
            b = V(F(:,2), :);
            c = V(F(:,3), :);
            faceArea = 0.5 * vecnorm(cross(b - a, c - a, 2), 2, 2);
            degMask = faceArea <= options.areaTolerance;

            if any(degMask)
                nDeg = nnz(degMask);
                if options.removeDegenerateFaces
                    F = F(~degMask, :);
                    report.messages{end+1} = sprintf('Removed %d degenerate faces.', nDeg); %#ok<AGROW>
                else
                    report.isValid = false;
                    report.messages{end+1} = sprintf('Found %d degenerate faces.', nDeg); %#ok<AGROW>
                end
            end

            % Edge analysis for manifold check
            edges = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
            edgesSorted = sort(edges, 2);
            [uniqueEdges, ~, edgeIds] = unique(edgesSorted, 'rows');
            edgeCount = accumarray(edgeIds, 1);
            boundaryMask = edgeCount == 1;
            nonManifoldMask = edgeCount > 2;
            nBoundaryEdges = nnz(boundaryMask);
            nNonManifoldEdges = nnz(nonManifoldMask);

            if nBoundaryEdges > 0
                report.messages{end+1} = sprintf('Mesh has %d boundary edges (open surface).', nBoundaryEdges); %#ok<AGROW>
                if options.requireWatertight
                    report.isValid = false;
                end
            end

            if nNonManifoldEdges > 0
                report.messages{end+1} = sprintf('Mesh has %d non-manifold edges.', nNonManifoldEdges); %#ok<AGROW>
                report.isValid = false;
            end

            % Volume check
            signedVolume = GeometryChecker.computeSignedVolume(V, F);
            if abs(signedVolume) < 1e-10
                report.messages{end+1} = 'Signed volume is near zero; geometry may be invalid.'; %#ok<AGROW>
                report.isValid = false;
            end

            % Bounding box check
            boundsMin = min(V, [], 1);
            boundsMax = max(V, [], 1);
            span = boundsMax - boundsMin;
            if any(span <= 1e-9)
                report.messages{end+1} = 'Mesh bounding box is collapsed in one or more directions.'; %#ok<AGROW>
                report.isValid = false;
            end

            % Update mesh with cleaned faces
            mesh.faces = F;
            mesh.nVertices = size(V, 1);
            mesh.nFaces = size(F, 1);

            report.nVertices = mesh.nVertices;
            report.nFaces = mesh.nFaces;
            report.boundaryEdges = nBoundaryEdges;
            report.nonManifoldEdges = nNonManifoldEdges;
            report.boundaryEdgeList = uniqueEdges(boundaryMask, :);
            report.nonManifoldEdgeList = uniqueEdges(nonManifoldMask, :);
            report.signedVolume = signedVolume;
            report.boundingBoxMin = boundsMin;
            report.boundingBoxMax = boundsMax;
        end

        function writeReportToFile(~, outputTag, report)
            txtFile = fullfile('Output', ['MeshSanity_', outputTag, '.txt']);
            fid = fopen(txtFile, 'w');
            if fid == -1
                warning('GeometryChecker:FileOpen', ...
                    'Could not open %s for writing.', txtFile);
                return;
            end

            span = report.boundingBoxMax - report.boundingBoxMin;
            fprintf(fid, 'NAUTILUS - Mesh Sanity Report\n');
            fprintf(fid, '%s\n', repmat('=', 1, 72));
            fprintf(fid, 'Mesh: %s\n', outputTag);
            fprintf(fid, 'Status: %s\n', GeometryChecker.statusLabel(report));
            fprintf(fid, '\nGeometry\n');
            fprintf(fid, '  Vertices:          %d\n', report.nVertices);
            fprintf(fid, '  Faces:             %d\n', report.nFaces);
            fprintf(fid, '  |Signed volume|:   %.6f m^3\n', abs(report.signedVolume));
            fprintf(fid, '  Boundary edges:    %d\n', report.boundaryEdges);
            fprintf(fid, '  Non-manifold:      %d\n', report.nonManifoldEdges);
            fprintf(fid, '\nBounding box\n');
            fprintf(fid, '  xmin/xmax:         %.6f / %.6f m\n', report.boundingBoxMin(1), report.boundingBoxMax(1));
            fprintf(fid, '  ymin/ymax:         %.6f / %.6f m\n', report.boundingBoxMin(2), report.boundingBoxMax(2));
            fprintf(fid, '  zmin/zmax:         %.6f / %.6f m\n', report.boundingBoxMin(3), report.boundingBoxMax(3));
            fprintf(fid, '  Span [L,B,D]:      %.6f  %.6f  %.6f m\n', span(1), span(2), span(3));

            if ~isempty(report.messages)
                fprintf(fid, '\nMessages\n');
                for k = 1:numel(report.messages)
                    fprintf(fid, '  - %s\n', report.messages{k});
                end
            end

            fclose(fid);
            fprintf('  Written: %s\n', txtFile);
        end

        function visualizeMesh(~, outputTag, mesh, report, options)
            if nargin < 5 || isempty(options)
                options = struct();
            end
            if ~isfield(options, 'visualizeInput'); options.visualizeInput = true; end
            if ~options.visualizeInput
                return;
            end
            if ~isfield(options, 'saveFigure'); options.saveFigure = true; end
            if ~isfield(options, 'showEdges'); options.showEdges = false; end
            if ~isfield(options, 'highlightProblemEdges'); options.highlightProblemEdges = true; end
            if ~isfield(options, 'faceAlpha'); options.faceAlpha = 0.96; end

            V = mesh.vertices;
            F = mesh.faces;
            span = report.boundingBoxMax - report.boundingBoxMin;
            fig = figure( ...
                'Color', 'w', ...
                'Position', [80 80 1360 860], ...
                'Renderer', 'painters', ...
                'InvertHardcopy', 'off');
            tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

            ax = nexttile(tl, 1, [2 1]);
            GeometryChecker.plotMesh(ax, V, F, report, options);
            view(ax, 34, 24);
            xlabel(ax, '$x\,[\mathrm{m}]$', 'Interpreter', 'latex');
            ylabel(ax, '$y\,[\mathrm{m}]$', 'Interpreter', 'latex');
            zlabel(ax, '$z\,[\mathrm{m}]$', 'Interpreter', 'latex');
            title(ax, GeometryChecker.latexHeading('Isometric Mesh View'), ...
                'Interpreter', 'latex');
            if exist('camlight', 'file') == 2
                camlight(ax, 'headlight');
                lighting(ax, 'gouraud');
            end

            ax = nexttile(tl, 2);
            GeometryChecker.plotMesh(ax, V, F, report, options);
            view(ax, 0, 90);
            xlabel(ax, '$x\,[\mathrm{m}]$', 'Interpreter', 'latex');
            ylabel(ax, '$y\,[\mathrm{m}]$', 'Interpreter', 'latex');
            zlabel(ax, '');
            title(ax, GeometryChecker.latexHeading('Plan View'), ...
                'Interpreter', 'latex');

            ax = nexttile(tl, 4);
            axis(ax, 'off');
            status = GeometryChecker.statusLabel(report);
            lines = { ...
                ['Status: ' status], ...
                '', ...
                sprintf('N_V = %d', report.nVertices), ...
                sprintf('N_F = %d', report.nFaces), ...
                sprintf('|V_s| = %.6f m^3', abs(report.signedVolume)), ...
                sprintf('N_boundary = %d', report.boundaryEdges), ...
                sprintf('N_non_manifold = %d', report.nonManifoldEdges), ...
                '', ...
                sprintf('Delta x = %.4f m', span(1)), ...
                sprintf('Delta y = %.4f m', span(2)), ...
                sprintf('Delta z = %.4f m', span(3))};
            if options.highlightProblemEdges
                lines{end+1} = '';
                lines{end+1} = 'Edge overlays';
                lines{end+1} = 'Boundary edges are highlighted in red';
                lines{end+1} = 'Non-manifold edges are highlighted in amber';
            end
            if ~isempty(report.messages)
                lines{end+1} = '';
                lines{end+1} = 'Notes';
                for k = 1:numel(report.messages)
                    lines{end+1} = ['- ' report.messages{k}]; %#ok<AGROW>
                end
            end
            text(ax, 0.02, 0.98, strjoin(lines, newline), ...
                'Units', 'normalized', ...
                'Interpreter', 'none', ...
                'FontName', 'Times New Roman', ...
                'FontSize', 12, ...
                'VerticalAlignment', 'top', ...
                'Color', [0.15 0.15 0.15]);

            title(tl, GeometryChecker.latexHeading('Geometry Gate Before Hydrostatics'), ...
                'Interpreter', 'latex', ...
                'FontSize', 18);

            GeometryChecker.hideAxesToolbars(fig);
            if options.saveFigure
                set(fig, 'PaperPositionMode', 'auto');
                outBase = fullfile('Output', ['MeshOverview_', outputTag]);
                print(fig, [outBase, '.png'], '-dpng', '-r300');
                try
                    print(fig, [outBase, '.pdf'], '-dpdf', '-painters');
                    fprintf('  Saved:   %s.pdf\n', outBase);
                catch
                    warning('GeometryChecker:PDFExportFailed', ...
                        'PDF export failed: %s. PNG saved successfully.', lasterr);
                end
                fprintf('  Saved:   %s.png\n', outBase);
            end
        end

        function printReport(~, hullName, report)
            fprintf('  Geometry checks for %s\n', hullName);
            fprintf('    Vertices: %d\n', report.nVertices);
            fprintf('    Faces:    %d\n', report.nFaces);
            fprintf('    |Volume|: %.6f m^3\n', abs(report.signedVolume));
            fprintf('    Boundary edges: %d\n', report.boundaryEdges);
            fprintf('    Non-manifold edges: %d\n', report.nonManifoldEdges);
            if isempty(report.messages)
                fprintf('    Checks: PASS\n');
            else
                for k = 1:numel(report.messages)
                    fprintf('    Note: %s\n', report.messages{k});
                end
                if report.isValid
                    fprintf('    Checks: PASS WITH WARNINGS\n');
                else
                    fprintf('    Checks: FAIL\n');
                end
            end
        end
    end

    methods (Static, Access = private)
        function plotMesh(ax, vertices, faces, report, options)
            edgeColor = 'none';
            if isfield(options, 'showEdges') && options.showEdges
                edgeColor = [0.25 0.25 0.25];
            end

            patch(ax, ...
                'Vertices', vertices, ...
                'Faces', faces, ...
                'FaceColor', [0.24 0.49 0.75], ...
                'FaceAlpha', options.faceAlpha, ...
                'EdgeColor', edgeColor, ...
                'LineWidth', 0.25, ...
                'AmbientStrength', 0.55, ...
                'DiffuseStrength', 0.75, ...
                'SpecularStrength', 0.08);
            hold(ax, 'on');
            if options.highlightProblemEdges
                GeometryChecker.plotEdgeOverlay(ax, vertices, report.boundaryEdgeList, [0.82 0.16 0.16], 1.8);
                GeometryChecker.plotEdgeOverlay(ax, vertices, report.nonManifoldEdgeList, [0.88 0.51 0.16], 2.1);
            end
            axis(ax, 'equal');
            axis(ax, 'tight');
            grid(ax, 'on');
            box(ax, 'on');
            ax.FontName = 'Times New Roman';
            ax.FontSize = 12;
            ax.LineWidth = 1.0;
            ax.TickDir = 'out';
            ax.TickLength = [0.015 0.015];
            ax.Layer = 'top';
            ax.GridAlpha = 0.15;
            ax.MinorGridAlpha = 0.08;
            ax.XMinorGrid = 'on';
            ax.YMinorGrid = 'on';
            ax.ZMinorGrid = 'on';
            view(ax, 3);
        end

        function plotEdgeOverlay(ax, vertices, edgeList, color, lineWidth)
            if isempty(edgeList)
                return;
            end

            xyz1 = vertices(edgeList(:,1), :);
            xyz2 = vertices(edgeList(:,2), :);
            x = [xyz1(:,1), xyz2(:,1), nan(size(edgeList, 1), 1)]';
            y = [xyz1(:,2), xyz2(:,2), nan(size(edgeList, 1), 1)]';
            z = [xyz1(:,3), xyz2(:,3), nan(size(edgeList, 1), 1)]';
            plot3(ax, x(:), y(:), z(:), '-', ...
                'Color', color, ...
                'LineWidth', lineWidth, ...
                'Clipping', 'on');
        end

        function hideAxesToolbars(fig)
            axs = findall(fig, 'Type', 'axes');
            for i = 1:numel(axs)
                ax = axs(i);
                if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar) && isprop(ax.Toolbar, 'Visible')
                    ax.Toolbar.Visible = 'off';
                end
            end
        end

        function txt = statusLabel(report)
            if report.isValid && isempty(report.messages)
                txt = 'PASS';
            elseif report.isValid
                txt = 'PASS WITH WARNINGS';
            else
                txt = 'FAIL';
            end
        end

        function txt = escapeLatex(txt)
            txt = strrep(txt, '\', '\textbackslash ');
            txt = strrep(txt, '_', '\_');
            txt = strrep(txt, '%', '\%');
            txt = strrep(txt, '#', '\#');
            txt = strrep(txt, '&', '\&');
            txt = strrep(txt, '{', '\{');
            txt = strrep(txt, '}', '\}');
            txt = strrep(txt, '^', '\^{}');
        end

        function txt = latexText(txt)
            txt = GeometryChecker.escapeLatex(txt);
            txt = strrep(txt, '-', '{-}');
            txt = strrep(txt, ' ', '\ ');
        end

        function txt = latexHeading(txt)
            txt = ['$\mathrm{' GeometryChecker.latexText(txt) '}$'];
        end

        function vol = computeSignedVolume(vertices, faces)
            v1 = vertices(faces(:,1), :);
            v2 = vertices(faces(:,2), :);
            v3 = vertices(faces(:,3), :);
            vol = sum(dot(v1, cross(v2, v3, 2), 2)) / 6;
        end
    end
end
