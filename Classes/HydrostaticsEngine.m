classdef HydrostaticsEngine
    %HydrostaticsEngine Run mesh hydrostatics draft sweep and GZ.
    %
    % FIXED BUGS:
    % 1. Added note about unused wave parameters (documented)
    % 2. Better error handling in initializeConditions

    properties
        cfg struct
    end

    methods
        function obj = HydrostaticsEngine(cfg)
            obj.cfg = cfg;
        end

        function results = run(obj, mesh, inputData, T_vec)
            nT = numel(T_vec);
            COND = obj.initializeConditions(inputData.VCG, T_vec);

            fprintf('  %4s  %8s  %12s  %8s  %8s  %8s  %-8s\n', ...
                'Step', 'T [m]', 'V [m^3]', 'KB [m]', 'KM [m]', 'Af [m^2]', 'Status');
            fprintf('  %s\n', repmat('-', 1, 68));

            nOK = 0;
            nFail = 0;

            for k = 1:nT
                try
                    [COND(k), ~] = Algorithms.geometriaMesh(COND(k), mesh);

                    V = COND(k).V;
                    KB = COND(k).B(3);
                    BM = COND(k).I(1) / max(V, 1e-9);
                    KM = KB + BM;
                    fprintf('  %4d  %8.3f  %12.3f  %8.4f  %8.4f  %8.3f  OK\n', ...
                        k, COND(k).T, V, KB, KM, COND(k).Af);
                    nOK = nOK + 1;
                catch ME
                    fprintf('  %4d  %8.3f  %12s  %8s  %8s  %8s  FAILED: %s\n', ...
                        k, T_vec(k), '-', '-', '-', '-', ME.message);
                    nFail = nFail + 1;
                end
            end

            fprintf('  %s\n', repmat('-', 1, 68));
            fprintf('  Computed: %d OK, %d failed\n\n', nOK, nFail);

            GZdata = [];
            if isfield(obj.cfg, 'computeGZ') && obj.cfg.computeGZ && nOK > 0
                refIdx = min(obj.cfg.refDraftIdx, nT);

                % Ensure the chosen reference draft was successfully computed.
                % If not, scan forward to the first valid draft.
                if isnan(COND(refIdx).V)
                    validIdxs = find(~isnan([COND.V]));
                    if isempty(validIdxs)
                        warning('HydrostaticsEngine:NoValidRefDraft', ...
                            'No valid draft found for GZ computation. Skipping GZ.');
                        refIdx = 0;
                    else
                        refIdx = validIdxs(1);
                        warning('HydrostaticsEngine:RefDraftAdjusted', ...
                            'refDraftIdx pointed to a failed draft; using first valid draft T = %.3f m instead.', ...
                            COND(refIdx).T);
                    end
                end

                if refIdx > 0
                    % Skip GZ if waterplane area is too small (near full immersion)
                    if COND(refIdx).Af < 0.15  % m^2 threshold
                        warning('HydrostaticsEngine:SmallWaterplaneArea', ...
                            'Waterplane area Af = %.3f m^2 at T = %.3f m is too small for stable GZ. Skipping GZ.', ...
                            COND(refIdx).Af, COND(refIdx).T);
                        GZdata = [];
                    else
                        fprintf('  Computing GZ curve at T = %.2f m ...\n', COND(refIdx).T);

                        COND_ref = struct();
                        COND_ref.V = COND(refIdx).V;
                        COND_ref.T = COND(refIdx).T;
                        COND_ref.G = [0, 0, inputData.VCG];

                        gzOptions = struct();
                        gzOptions.heelAngles = obj.cfg.tetaT_vec;
                        gzOptions.useMesh = true;
                        gzOptions.tolVolume = obj.cfg.tolVolume;
                        gzOptions.tolTrim = obj.cfg.tolTrim;
                        gzOptions.verbose = true;

                        GZdata = Algorithms.computeGZcurve(mesh, COND_ref, gzOptions);
                        fprintf('\n');
                    end
                end
            end

            results = struct();
            results.COND = COND;
            results.GZdata = GZdata;
            results.nOK = nOK;
            results.nFail = nFail;
        end
    end

    methods (Access = private)
        function COND = initializeConditions(obj, vcg, T_vec)
            nT = numel(T_vec);

            condTemplate = struct();
            condTemplate.V = NaN;
            condTemplate.B = [NaN NaN NaN];
            condTemplate.G = [NaN NaN vcg];
            condTemplate.Af = NaN;
            condTemplate.I = [NaN NaN];
            condTemplate.Lwl = NaN;
            condTemplate.Bwl = NaN;
            condTemplate.Ams = NaN;
            condTemplate.Ws = NaN;
            condTemplate.tetaT = 0;
            condTemplate.tetaL = 0;
            condTemplate.T = NaN;
            condTemplate.F = [NaN NaN NaN];
            condTemplate.GZ = NaN;
            
            % NOTE: Wave parameters are stored but not used in current implementation.
            % These are placeholders for future wave-added mass/damping calculations
            % or seakeeping analysis. The hydrostatics computation assumes calm water.
            if isfield(obj.cfg, 'Lwave')
                condTemplate.Lwave = obj.cfg.Lwave;
            else
                condTemplate.Lwave = 0;
            end
            if isfield(obj.cfg, 'Hwave')
                condTemplate.Hwave = obj.cfg.Hwave;
            else
                condTemplate.Hwave = 0;
            end
            if isfield(obj.cfg, 'Xcrest')
                condTemplate.Xcrest = obj.cfg.Xcrest;
            else
                condTemplate.Xcrest = 0;
            end
            if isfield(obj.cfg, 'Wavetype')
                condTemplate.Wavetype = obj.cfg.Wavetype;
            else
                condTemplate.Wavetype = 2;
            end

            COND = repmat(condTemplate, 1, nT);
            for k = 1:nT
                COND(k).T = T_vec(k);
            end
        end
    end
end
