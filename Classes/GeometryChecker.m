classdef GeometryChecker
    % GeometryChecker
    % Run mesh sanity checks before hydrostatic analysis.
    %
    % Main checks:
    %   - structural validity of vertices/faces arrays
    %   - face-index integrity
    %   - degenerate-face detection
    %   - boundary-edge detection (open surface)
    %   - non-manifold-edge detection
    %   - signed-volume sanity
    %   - bounding-box sanity
    %
    % Notes:
    %   - Boundary edges are warnings unless requireWatertight = true
    %   - Near-zero signed volume is only treated as fatal when watertightness
    %     is required; otherwise it is reported as a warning
    %   - Degenerate faces can optionally be removed automatically

    methods
        function [mesh, report] = check(~, mesh, options)
            if nargin < 3 || isempty(options)
                options = struct();
            end
            options = GeometryChecker.applyDefaults(options, struct( ...
                'areaTolerance', 1e-12, ...
                'indexTolerance', 1e-9, ...
                'removeDegenerateFaces', true, ...
                'requireWatertight', false, ...
                'volumeTolerance', 1e-10, ...
                'spanTolerance', 1e-9));

            report = GeometryChecker.initializeReport();

            if ~isstruct(mesh) || ~isfield(mesh, 'vertices') || ~isfield(mesh, 'faces')
                error('GeometryChecker:InvalidMesh', ...
                    'Mesh must be a struct with fields "vertices" and "faces".');
            end

            V = mesh.vertices;
            F = mesh.faces;

            validateattributes(V, {'numeric'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
                'GeometryChecker.check', 'mesh.vertices');
            validateattributes(F, {'numeric'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
                'GeometryChecker.check', 'mesh.faces');

            report.originalVertices = size(V, 1);
            report.originalFaces = size(F, 1);

            if isempty(V)
                error('GeometryChecker:EmptyVertices', ...
                    'Mesh contains no vertices.');
            end
            if isempty(F)
                error('GeometryChecker:EmptyFaces', ...
                    'Mesh contains no faces.');
            end

            nV = size(V, 1);

            % Validate face indices are integer-valued within tolerance.
            idxError = abs(F - round(F));
            badIntegerMask = idxError > options.indexTolerance;
            if any(badIntegerMask, 'all')
                [rowBad, colBad] = find(badIntegerMask, 1, 'first');
                error('GeometryChecker:NonIntegerFaceIndices', ...
                    ['Mesh faces contain non-integer vertex indices. ' ...
                     'Example: faces(%d,%d) = %.16g.'], ...
                    rowBad, colBad, F(rowBad, colBad));
            end

            % Safe to round now because we already validated closeness to integer.
            F = round(F);

            if any(F(:) < 1) || any(F(:) > nV)
                invalidFaces = find(any(F < 1, 2) | any(F > nV, 2));
                exampleFace = invalidFaces(1);
                error('GeometryChecker:FaceIndexOutOfRange', ...
                    ['Mesh face indices exceed vertex bounds [1, %d]. ' ...
                     'Found %d invalid faces; example face index: %d.'], ...
                    nV, numel(invalidFaces), exampleFace);
            end

            % Reject faces with repeated vertex indices before area computation.
            repeatedIndexMask = F(:,1) == F(:,2) | F(:,1) == F(:,3) | F(:,2) == F(:,3);
            if any(repeatedIndexMask)
                nRep = nnz(repeatedIndexMask);
                if options.removeDegenerateFaces
                    F = F(~repeatedIndexMask, :);
                    report.messages{end+1} = sprintf( ...
                        'Removed %d faces with repeated vertex indices.', nRep);
                else
                    report.isValid = false;
                    report.messages{end+1} = sprintf( ...
                        'Found %d faces with repeated vertex indices.', nRep);
                end
            end

            if isempty(F)
                report.isValid = false;
                report.messages{end+1} = 'No valid faces remain after repeated-index filtering.';
                mesh.faces = F;
                mesh.nVertices = size(V, 1);
                mesh.nFaces = 0;
                report = GeometryChecker.finalizeReport(report, mesh, V, F, options);
                return;
            end

            % Compute face areas.
            a = V(F(:,1), :);
            b = V(F(:,2), :);
            c = V(F(:,3), :);
            faceArea = 0.5 * vecnorm(cross(b - a, c - a, 2), 2, 2);
            degMask = faceArea <= options.areaTolerance;

            report.degenerateFacesFound = nnz(degMask);

            if any(degMask)
                nDeg = nnz(degMask);
                if options.removeDegenerateFaces
                    F = F(~degMask, :);
                    report.messages{end+1} = sprintf('Removed %d degenerate faces.', nDeg);
                else
                    report.isValid = false;
                    report.messages{end+1} = sprintf('Found %d degenerate faces.', nDeg);
                end
            end

            if isempty(F)
                report.isValid = false;
                report.messages{end+1} = 'No valid faces remain after degenerate-face removal.';
                mesh.faces = F;
                mesh.nVertices = size(V, 1);
                mesh.nFaces = 0;
                report = GeometryChecker.finalizeReport(report, mesh, V, F, options);
                return;
            end

            % Edge manifold analysis.
            edges = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
            edgesSorted = sort(edges, 2);
            [uniqueEdges, ~, edgeIds] = unique(edgesSorted, 'rows');
            edgeCount = accumarray(edgeIds, 1);

            boundaryMask = edgeCount == 1;
            nonManifoldMask = edgeCount > 2;

            nBoundaryEdges = nnz(boundaryMask);
            nNonManifoldEdges = nnz(nonManifoldMask);

            if nBoundaryEdges > 0
                report.messages{end+1} = sprintf( ...
                    'Mesh has %d boundary edges (open surface).', nBoundaryEdges);
                if options.requireWatertight
                    report.isValid = false;
                end
            end

            if nNonManifoldEdges > 0
                report.messages{end+1} = sprintf( ...
                    'Mesh has %d non-manifold edges.', nNonManifoldEdges);
                report.isValid = false;
            end

            % Signed volume.
            signedVolume = GeometryChecker.computeSignedVolume(V, F);

            if abs(signedVolume) < options.volumeTolerance
                if options.requireWatertight
                    report.messages{end+1} = ...
                        'Signed volume is near zero; watertight closed-body assumption is violated.';
                    report.isValid = false;
                else
                    report.messages{end+1} = ...
                        'Signed volume is near zero; this may be acceptable for an open surface but is unsuitable for closed-body hydrostatics.';
                end
            end

            % Bounding box.
            boundsMin = min(V, [], 1);
            boundsMax = max(V, [], 1);
            span = boundsMax - boundsMin;

            if any(span <= options.spanTolerance)
                report.messages{end+1} = ...
                    'Mesh bounding box is collapsed in one or more directions.';
                report.isValid = false;
            end

            % Update mesh.
            mesh.faces = F;
            mesh.nVertices = size(V, 1);
            mesh.nFaces = size(F, 1);

            % Populate report.
            report.nVertices = mesh.nVertices;
            report.nFaces = mesh.nFaces;
            report.boundaryEdges = nBoundaryEdges;
            report.nonManifoldEdges = nNonManifoldEdges;
            report.boundaryEdgeList = uniqueEdges(boundaryMask, :);
            report.nonManifoldEdgeList = uniqueEdges(nonManifoldMask, :);
            report.signedVolume = signedVolume;
            report.boundingBoxMin = boundsMin;
            report.boundingBoxMax = boundsMax;
            report.boundingBoxSpan = span;

            report.removedFaces = report.originalFaces - report.nFaces;
        end

        function writeReportToFile(~, outputTag, report)
            outDir = 'Output';
            if exist(outDir, 'dir') ~= 7
                mkdir(outDir);
            end

            txtFile = fullfile(outDir, ['MeshSanity_', outputTag, '.txt']);
            fid = fopen(txtFile, 'w');
            if fid == -1
                warning('GeometryChecker:FileOpen', ...
                    'Could not open "%s" for writing.', txtFile);
                return;
            end

            cleaner = onCleanup(@() fclose(fid));

            if isfield(report, 'boundingBoxSpan')
                span = report.boundingBoxSpan;
            else
                span = report.boundingBoxMax - report.boundingBoxMin;
            end

            fprintf(fid, 'NAUTILUS - Mesh Sanity Report\n');
            fprintf(fid, '%s\n', repmat('=', 1, 72));
            fprintf(fid, 'Mesh: %s\n', outputTag);
            fprintf(fid, 'Status: %s\n', GeometryChecker.statusLabel(report));

            fprintf(fid, '\nGeometry\n');
            fprintf(fid, '  Original vertices:  %d\n', GeometryChecker.getFieldOr(report, 'originalVertices', report.nVertices));
            fprintf(fid, '  Original faces:     %d\n', GeometryChecker.getFieldOr(report, 'originalFaces', report.nFaces));
            fprintf(fid, '  Final vertices:     %d\n', report.nVertices);
            fprintf(fid, '  Final faces:        %d\n', report.nFaces);
            fprintf(fid, '  Removed faces:      %d\n', GeometryChecker.getFieldOr(report, 'removedFaces', 0));
            fprintf(fid, '  |Signed volume|:    %.6f m^3\n', abs(report.signedVolume));
            fprintf(fid, '  Boundary edges:     %d\n', report.boundaryEdges);
            fprintf(fid, '  Non-manifold edges: %d\n', report.nonManifoldEdges);

            fprintf(fid, '\nBounding box\n');
            fprintf(fid, '  xmin/xmax:          %.6f / %.6f m\n', report.boundingBoxMin(1), report.boundingBoxMax(1));
            fprintf(fid, '  ymin/ymax:          %.6f / %.6f m\n', report.boundingBoxMin(2), report.boundingBoxMax(2));
            fprintf(fid, '  zmin/zmax:          %.6f / %.6f m\n', report.boundingBoxMin(3), report.boundingBoxMax(3));
            fprintf(fid, '  Span [L,B,D]:       %.6f  %.6f  %.6f m\n', span(1), span(2), span(3));

            if ~isempty(report.messages)
                fprintf(fid, '\nMessages\n');
                for k = 1:numel(report.messages)
                    fprintf(fid, '  - %s\n', report.messages{k});
                end
            end

            fprintf('  Written: %s\n', txtFile);
        end

        function visualizeMesh(~, outputTag, mesh, report, options)
            if nargin < 5 || isempty(options)
                options = struct();
            end
            options = GeometryChecker.applyDefaults(options, struct( ...
                'visualizeInput', true, ...
                'saveFigure', true, ...
                'showEdges', false, ...
                'highlightProblemEdges', true, ...
                'faceAlpha', 0.96));

            if ~options.visualizeInput
                return;
            end

            if ~isfield(mesh, 'vertices') || ~isfield(mesh, 'faces') || isempty(mesh.faces)
                warning('GeometryChecker:NoMeshToVisualize', ...
                    'Visualization skipped because the mesh is empty or invalid.');
                return;
            end

            outDir = 'Output';
            if exist(outDir, 'dir') ~= 7
                mkdir(outDir);
            end

            V = mesh.vertices;
            F = mesh.faces;
            span = report.boundingBoxMax - report.boundingBoxMin;

            fig = figure( ...
                'Color', 'w', ...
                'Position', [80 80 1360 860], ...
                'Renderer', 'painters', ...
                'InvertHardcopy', 'off');

            tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

            ax1 = nexttile(tl, 1, [2 1]);
            GeometryChecker.plotMesh(ax1, V, F, report, options);
            view(ax1, 34, 24);
            xlabel(ax1, '$x\,[\mathrm{m}]$', 'Interpreter', 'latex');
            ylabel(ax1, '$y\,[\mathrm{m}]$', 'Interpreter', 'latex');
            zlabel(ax1, '$z\,[\mathrm{m}]$', 'Interpreter', 'latex');
            title(ax1, GeometryChecker.latexHeading('Isometric Mesh View'), 'Interpreter', 'latex');

            if exist('camlight', 'file') == 2
                camlight(ax1, 'headlight');
                lighting(ax1, 'gouraud');
            end

            ax2 = nexttile(tl, 2);
            GeometryChecker.plotMesh(ax2, V, F, report, options);
            view(ax2, 0, 90);
            xlabel(ax2, '$x\,[\mathrm{m}]$', 'Interpreter', 'latex');
            ylabel(ax2, '$y\,[\mathrm{m}]$', 'Interpreter', 'latex');
            zlabel(ax2, '');
            title(ax2, GeometryChecker.latexHeading('Plan View'), 'Interpreter', 'latex');

            ax3 = nexttile(tl, 4);
            axis(ax3, 'off');

            lines = { ...
                ['Status: ' GeometryChecker.statusLabel(report)], ...
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

            text(ax3, 0.02, 0.98, strjoin(lines, newline), ...
                'Units', 'normalized', ...
                'Interpreter', 'none', ...
                'FontName', 'Times New Roman', ...
                'FontSize', 12, ...
                'VerticalAlignment', 'top', ...
                'Color', [0.15 0.15 0.15]);

            title(tl, GeometryChecker.latexHeading('Geometry Gate Before Hydrostatics'), ...
                'Interpreter', 'latex', 'FontSize', 18);

            GeometryChecker.hideAxesToolbars(fig);

            if options.saveFigure
                set(fig, 'PaperPositionMode', 'auto');
                outBase = fullfile(outDir, ['MeshOverview_', outputTag]);

                print(fig, [outBase, '.png'], '-dpng', '-r300');
                fprintf('  Saved:   %s.png\n', outBase);

                try
                    print(fig, [outBase, '.pdf'], '-dpdf', '-painters');
                    fprintf('  Saved:   %s.pdf\n', outBase);
                catch ME
                    warning('GeometryChecker:PDFExportFailed', ...
                        'PDF export failed for "%s": %s. PNG saved successfully.', outBase, ME.message);
                end
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
                fprintf('    Checks: %s\n', GeometryChecker.statusLabel(report));
            end
        end
    end

    methods (Static, Access = private)
        function report = initializeReport()
            report = struct();
            report.isValid = true;
            report.messages = {};
            report.originalVertices = 0;
            report.originalFaces = 0;
            report.nVertices = 0;
            report.nFaces = 0;
            report.removedFaces = 0;
            report.degenerateFacesFound = 0;
            report.boundaryEdges = 0;
            report.nonManifoldEdges = 0;
            report.boundaryEdgeList = zeros(0,2);
            report.nonManifoldEdgeList = zeros(0,2);
            report.signedVolume = NaN;
            report.boundingBoxMin = [NaN NaN NaN];
            report.boundingBoxMax = [NaN NaN NaN];
            report.boundingBoxSpan = [NaN NaN NaN];
        end

        function report = finalizeReport(report, mesh, V, F, ~)
            report.nVertices = size(V, 1);
            report.nFaces = size(F, 1);

            if ~isempty(V)
                report.boundingBoxMin = min(V, [], 1);
                report.boundingBoxMax = max(V, [], 1);
                report.boundingBoxSpan = report.boundingBoxMax - report.boundingBoxMin;
            end

            if ~isempty(F)
                report.signedVolume = GeometryChecker.computeSignedVolume(V, F);
            else
                report.signedVolume = 0;
            end

            mesh.nVertices = size(V, 1); %#ok<NASGU>
            mesh.nFaces = size(F, 1); %#ok<NASGU>
        end

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

            x = [xyz1(:,1), xyz2(:,1), nan(size(edgeList,1),1)]';
            y = [xyz1(:,2), xyz2(:,2), nan(size(edgeList,1),1)]';
            z = [xyz1(:,3), xyz2(:,3), nan(size(edgeList,1),1)]';

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
            if isempty(vertices) || isempty(faces)
                vol = 0;
                return;
            end
            v1 = vertices(faces(:,1), :);
            v2 = vertices(faces(:,2), :);
            v3 = vertices(faces(:,3), :);
            vol = sum(dot(v1, cross(v2, v3, 2), 2)) / 6;
        end

        function S = applyDefaults(S, defaults)
            fn = fieldnames(defaults);
            for k = 1:numel(fn)
                if ~isfield(S, fn{k})
                    S.(fn{k}) = defaults.(fn{k});
                end
            end
        end

        function value = getFieldOr(S, name, fallback)
            if isfield(S, name)
                value = S.(name);
            else
                value = fallback;
            end
        end
    end
end