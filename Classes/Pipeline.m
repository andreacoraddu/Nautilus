classdef Pipeline
    % Pipeline
    % Orchestrate mesh loading, geometry checks, hydrostatics, and outputs.
    %
    % Workflow:
    %   1) load hull input + mesh
    %   2) run geometry sanity checks
    %   3) visualize input mesh
    %   4) run hydrostatics / equilibrium / GZ
    %   5) generate plots and tables
    %
    % Design notes:
    %   - Geometry sanity failures stop the workflow before hydrostatics.
    %   - Geometry warnings are allowed and are carried into reporting.
    %   - Configuration defaults are centralized in withDefaults().

    properties
        cfg struct
    end

    methods
        function obj = Pipeline(cfg)
            if nargin < 1
                cfg = struct();
            end
            obj.cfg = Pipeline.withDefaults(cfg);
        end

        function runForHull(obj, hullName, T_vec)
            arguments
                obj
                hullName (1,:) char
                T_vec (:,1) double
            end

            T_vec = T_vec(:).';  % enforce row vector for consistent printing/indexing
            Pipeline.validateDraftVector(T_vec);
            outputDir = Pipeline.ensureOutputDirectory();

            fprintf('\n');
            fprintf('  +-----------------------------------------------------------+\n');
            fprintf('  |  Hull   : %-49s|\n', hullName);
            fprintf('  |  Drafts : %7.3f m -> %7.3f m  (%3d points)               |\n', ...
                T_vec(1), T_vec(end), numel(T_vec));
            fprintf('  +-----------------------------------------------------------+\n\n');

            %------------------------------------------------------------------
            % 1) Load geometry and input data
            %------------------------------------------------------------------
            fprintf('  [1/5] Loading geometry as mesh...\n');
            loader = GeometryLoader(obj.cfg);

            try
                [mesh, inputData, baseName] = loader.loadHullMesh(hullName);
            catch ME
                error('Pipeline:LoadFailed', ...
                    'Failed while loading hull "%s".\nOriginal error: %s', ...
                    hullName, ME.message);
            end

            fprintf('        Mesh loaded: %d vertices, %d faces\n\n', ...
                size(mesh.vertices,1), size(mesh.faces,1));

            outputTag = [baseName, '_mesh'];

            %------------------------------------------------------------------
            % 2) Geometry sanity checks
            %------------------------------------------------------------------
            fprintf('  [2/5] Running geometry sanity checks...\n');
            checker = GeometryChecker();

            try
                [mesh, report] = checker.check(mesh, obj.cfg.sanity);
            catch ME
                error('Pipeline:GeometryCheckFailed', ...
                    'Geometry checks crashed for hull "%s".\nOriginal error: %s', ...
                    baseName, ME.message);
            end

            checker.printReport(baseName, report);
            checker.writeReportToFile(outputTag, report);
            fprintf('\n');

            %------------------------------------------------------------------
            % 3) Visualize input geometry
            %------------------------------------------------------------------
            fprintf('  [3/5] Rendering pre-analysis mesh overview...\n');

            visOptions = struct();
            visOptions.visualizeInput      = obj.cfg.meshOptions.visualizeInput;
            visOptions.saveFigure          = obj.cfg.meshOptions.saveInputFigure;
            visOptions.showEdges           = obj.cfg.meshOptions.showInputEdges;
            visOptions.highlightProblemEdges = true;
            visOptions.faceAlpha           = 0.96;

            try
                checker.visualizeMesh(outputTag, mesh, report, visOptions);
            catch ME
                warning('Pipeline:VisualizationFailed', ...
                    'Input mesh visualization failed for "%s": %s', ...
                    baseName, ME.message);
            end
            fprintf('\n');

            % Stop here if geometry is not acceptable
            if ~report.isValid
                error('Pipeline:SanityCheckFailed', ...
                    'Geometry sanity checks failed for "%s". Hydrostatics was not started.', ...
                    baseName);
            end

            %------------------------------------------------------------------
            % 4) Hydrostatics / GZ
            %------------------------------------------------------------------
            fprintf('  [4/5] Running hydrostatics and GZ...\n');
            engine = HydrostaticsEngine(obj.cfg);

            try
                results = engine.run(mesh, inputData, T_vec);
            catch ME
                error('Pipeline:HydrostaticsFailed', ...
                    'Hydrostatics stage failed for "%s".\nOriginal error: %s', ...
                    baseName, ME.message);
            end

            %------------------------------------------------------------------
            % 5) Post-processing
            %------------------------------------------------------------------
            fprintf('  [5/5] Post-processing plots and tables...\n');
            post = PostProcessor(obj.cfg);

            try
                post.generate(outputTag, results, mesh, inputData);
            catch ME
                error('Pipeline:PostProcessingFailed', ...
                    'Post-processing failed for "%s".\nOriginal error: %s', ...
                    baseName, ME.message);
            end

            fprintf('  Post-processing done.\n');
            fprintf('  Output folder: %s\n\n', outputDir);
        end
    end

    methods (Static, Access = private)
        function cfg = withDefaults(cfg)
            if nargin < 1 || isempty(cfg)
                cfg = struct();
            end

            %-----------------------------
            % Global analysis options
            %-----------------------------
            if ~isfield(cfg, 'computeGZ');   cfg.computeGZ = true;    end
            if ~isfield(cfg, 'tetaT_vec');   cfg.tetaT_vec = 0:5:60;  end
            if ~isfield(cfg, 'refDraftIdx'); cfg.refDraftIdx = 1;     end
            if ~isfield(cfg, 'tolVolume');   cfg.tolVolume = 1e-5;    end
            if ~isfield(cfg, 'tolTrim');     cfg.tolTrim = 1e-4;      end
            if ~isfield(cfg, 'rho');         cfg.rho = 1025;          end
            if ~isfield(cfg, 'savePlots');   cfg.savePlots = true;    end

            %-----------------------------
            % Wave / operating condition options
            %-----------------------------
            if ~isfield(cfg, 'Lwave');    cfg.Lwave = 0; end
            if ~isfield(cfg, 'Hwave');    cfg.Hwave = 0; end
            if ~isfield(cfg, 'Xcrest');   cfg.Xcrest = 0; end
            if ~isfield(cfg, 'Wavetype'); cfg.Wavetype = 2; end

            %-----------------------------
            % Mesh loading / output options
            %-----------------------------
            if ~isfield(cfg, 'autoSaveMesh'); cfg.autoSaveMesh = true; end
            if ~isfield(cfg, 'cleanMesh'); cfg.cleanMesh = true; end
            if ~isfield(cfg, 'allowZeroBasedObj'); cfg.allowZeroBasedObj = false; end
            if ~isfield(cfg, 'geometricTolerance'); cfg.geometricTolerance = 1e-12; end

            if ~isfield(cfg, 'meshOptions') || ~isstruct(cfg.meshOptions)
                cfg.meshOptions = struct();
            end

            if ~isfield(cfg.meshOptions, 'plotSubmerged');   cfg.meshOptions.plotSubmerged = false; end
            if ~isfield(cfg.meshOptions, 'visualizeInput');  cfg.meshOptions.visualizeInput = true; end
            if ~isfield(cfg.meshOptions, 'saveInputFigure'); cfg.meshOptions.saveInputFigure = cfg.savePlots; end
            if ~isfield(cfg.meshOptions, 'showInputEdges');  cfg.meshOptions.showInputEdges = false; end

            %-----------------------------
            % Geometry sanity options
            %-----------------------------
            if ~isfield(cfg, 'sanity') || ~isstruct(cfg.sanity)
                cfg.sanity = struct();
            end

            if ~isfield(cfg.sanity, 'areaTolerance');         cfg.sanity.areaTolerance = 1e-12; end
            if ~isfield(cfg.sanity, 'indexTolerance');        cfg.sanity.indexTolerance = 1e-9; end
            if ~isfield(cfg.sanity, 'removeDegenerateFaces'); cfg.sanity.removeDegenerateFaces = true; end
            if ~isfield(cfg.sanity, 'requireWatertight');     cfg.sanity.requireWatertight = false; end
            if ~isfield(cfg.sanity, 'volumeTolerance');       cfg.sanity.volumeTolerance = 1e-10; end
            if ~isfield(cfg.sanity, 'spanTolerance');         cfg.sanity.spanTolerance = 1e-9; end

            %-----------------------------
            % Final defensive validation
            %-----------------------------
            validateattributes(cfg.tolVolume, {'numeric'}, {'scalar','positive','finite'});
            validateattributes(cfg.tolTrim, {'numeric'}, {'scalar','positive','finite'});
            validateattributes(cfg.rho, {'numeric'}, {'scalar','positive','finite'});
            validateattributes(cfg.refDraftIdx, {'numeric'}, {'scalar','integer','positive','finite'});
        end

        function validateDraftVector(T_vec)
            if isempty(T_vec)
                error('Pipeline:EmptyDraftRange', ...
                    'T_vec must be a non-empty vector of draft values.');
            end

            if any(~isfinite(T_vec))
                error('Pipeline:InvalidDraftRange', ...
                    'T_vec must contain only finite numeric values.');
            end

            if any(T_vec <= 0)
                error('Pipeline:NegativeDraft', ...
                    'All draft values must be strictly positive (minimum found: %.6f m).', ...
                    min(T_vec));
            end

            if numel(T_vec) > 1 && any(diff(T_vec) <= 0)
                warning('Pipeline:UnsortedDrafts', ...
                    'T_vec is not strictly increasing; downstream summaries may be harder to interpret.');
            end
        end

        function outputDir = ensureOutputDirectory()
            outputDir = 'Output';
            if exist(outputDir, 'dir') ~= 7
                mkdir(outputDir);
            end
        end
    end
end