classdef Pipeline
    %Pipeline Orchestrate mesh loading, checks, hydrostatics, and outputs.

    properties
        cfg struct
    end

    methods
        function obj = Pipeline(cfg)
            obj.cfg = Pipeline.withDefaults(cfg);
        end

        function runForHull(obj, hullName, T_vec)
            % --- Validate draft range ---
            if isempty(T_vec)
                error('Pipeline:EmptyDraftRange', 'T_vec must be a non-empty vector of draft values.');
            end
            if any(T_vec <= 0)
                error('Pipeline:NegativeDraft', 'All draft values must be positive (got min = %.4f m).', min(T_vec));
            end
            if numel(T_vec) > 1 && any(diff(T_vec) <= 0)
                warning('Pipeline:UnsortedDrafts', 'T_vec is not strictly increasing; results may be confusing.');
            end

            % --- Ensure Output directory exists ---
            if ~exist('Output', 'dir')
                mkdir('Output');
            end

            fprintf('\n');
            fprintf('  +-----------------------------------------------------------+\n');
            fprintf('  |  Hull : %-50s|\n', hullName);
            fprintf('  |  Drafts : %7.3f m -> %7.3f m  (%3d points)               |\n', T_vec(1), T_vec(end), numel(T_vec));
            fprintf('  +-----------------------------------------------------------+\n\n');

            fprintf('  [1/5] Loading geometry as mesh...\n');
            loader = GeometryLoader(obj.cfg);
            [mesh, inputData, baseName] = loader.loadHullMesh(hullName);
            fprintf('        Mesh loaded: %d vertices, %d faces\n\n', size(mesh.vertices,1), size(mesh.faces,1));

            outputTag = [baseName, '_mesh'];

            fprintf('  [2/5] Running geometry sanity checks...\n');
            checker = GeometryChecker();
            [mesh, report] = checker.check(mesh, obj.cfg.sanity);
            checker.printReport(baseName, report);
            checker.writeReportToFile(outputTag, report);
            fprintf('\n');

            fprintf('  [3/5] Rendering pre-analysis mesh overview...\n');
            visOptions = struct();
            visOptions.visualizeInput = obj.cfg.meshOptions.visualizeInput;
            visOptions.saveFigure = obj.cfg.meshOptions.saveInputFigure;
            visOptions.showEdges = obj.cfg.meshOptions.showInputEdges;
            checker.visualizeMesh(outputTag, mesh, report, visOptions);
            fprintf('\n');

            if ~report.isValid
                error('Pipeline:SanityCheckFailed', ...
                    'Geometry sanity checks failed for %s.', baseName);
            end

            fprintf('  [4/5] Running hydrostatics and GZ...\n');
            engine = HydrostaticsEngine(obj.cfg);
            results = engine.run(mesh, inputData, T_vec);

            fprintf('  [5/5] Post-processing plots and tables...\n');
            post = PostProcessor(obj.cfg);
            post.generate(outputTag, results, mesh, inputData);
            fprintf('  Post-processing done.\n\n');
        end
    end

    methods (Static, Access = private)
        function cfg = withDefaults(cfg)
            if nargin < 1 || isempty(cfg)
                cfg = struct();
            end

            if ~isfield(cfg, 'computeGZ'); cfg.computeGZ = true; end
            if ~isfield(cfg, 'tetaT_vec'); cfg.tetaT_vec = 0:5:60; end
            if ~isfield(cfg, 'refDraftIdx'); cfg.refDraftIdx = 1; end
            if ~isfield(cfg, 'tolVolume'); cfg.tolVolume = 1e-5; end
            if ~isfield(cfg, 'tolTrim'); cfg.tolTrim = 1e-4; end
            if ~isfield(cfg, 'rho'); cfg.rho = 1025; end
            if ~isfield(cfg, 'savePlots'); cfg.savePlots = true; end
            if ~isfield(cfg, 'Lwave'); cfg.Lwave = 0; end
            if ~isfield(cfg, 'Hwave'); cfg.Hwave = 0; end
            if ~isfield(cfg, 'Xcrest'); cfg.Xcrest = 0; end
            if ~isfield(cfg, 'Wavetype'); cfg.Wavetype = 2; end
            if ~isfield(cfg, 'autoSaveMesh'); cfg.autoSaveMesh = true; end
            if ~isfield(cfg, 'meshOptions'); cfg.meshOptions = struct(); end
            if ~isfield(cfg.meshOptions, 'plotSubmerged'); cfg.meshOptions.plotSubmerged = false; end
            if ~isfield(cfg.meshOptions, 'visualizeInput'); cfg.meshOptions.visualizeInput = true; end
            if ~isfield(cfg.meshOptions, 'saveInputFigure'); cfg.meshOptions.saveInputFigure = cfg.savePlots; end
            if ~isfield(cfg.meshOptions, 'showInputEdges'); cfg.meshOptions.showInputEdges = false; end
            if ~isfield(cfg, 'sanity'); cfg.sanity = struct(); end
            if ~isfield(cfg.sanity, 'areaTolerance'); cfg.sanity.areaTolerance = 1e-12; end
            if ~isfield(cfg.sanity, 'removeDegenerateFaces'); cfg.sanity.removeDegenerateFaces = true; end
            if ~isfield(cfg.sanity, 'requireWatertight'); cfg.sanity.requireWatertight = false; end
        end
    end
end
