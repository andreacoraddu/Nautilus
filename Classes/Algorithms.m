classdef Algorithms
    %Algorithms Core mesh hydrostatics and stability algorithms.
    %
    % FIXED BUGS:
    % 1. computeRightingLever - Fixed sign logic using proper cross product
    % 2. computeWaterplaneProperties - Fixed to handle edge cases and return NaN for invalid
    % 3. computeTrimCorrection - Fixed sign convention
    % 4. clipTriangleToWaterplane - Fixed edge cases

    methods (Static)
        function [COND1, submergedMesh] = geometriaMesh(COND, mesh, varargin)
            % varargin absorbs any legacy plottami argument; it is no longer used.

            validateattributes(COND.T, {'numeric'}, {'scalar', 'real'}, 'Algorithms.geometriaMesh', 'T');
            validateattributes(COND.tetaT, {'numeric'}, {'scalar', 'real'}, 'Algorithms.geometriaMesh', 'tetaT');
            validateattributes(COND.tetaL, {'numeric'}, {'scalar', 'real'}, 'Algorithms.geometriaMesh', 'tetaL');

            if ~isstruct(mesh) || ~isfield(mesh, 'vertices') || ~isfield(mesh, 'faces')
                error('Algorithms:InvalidMesh', 'mesh must include vertices and faces.');
            end

            vertices = mesh.vertices;
            faces = mesh.faces;
            nFaces = size(faces, 1);

            tetaT = COND.tetaT * pi / 180;
            tetaL = COND.tetaL * pi / 180;
            TM = COND.T;

            R_T = [1, 0, 0; 0, cos(tetaT), -sin(tetaT); 0, sin(tetaT), cos(tetaT)];
            R_L = [cos(tetaL), 0, sin(tetaL); 0, 1, 0; -sin(tetaL), 0, cos(tetaL)];
            R = R_L * R_T;

            vertices_world = (R * vertices')';
            vertices_world(:, 3) = vertices_world(:, 3) - TM;

            maxSubmergedFaces = nFaces * 3;
            subV = zeros(maxSubmergedFaces * 3, 3);
            subF = zeros(maxSubmergedFaces, 3);
            nSub = 0;
            waterEdges = zeros(nFaces * 2, 3);
            nWaterEdges = 0;
            wettedArea = 0;

            for f = 1:nFaces
                vIdx = faces(f, :);
                v = vertices_world(vIdx, :);

                z = v(:, 3);
                zMin = min(z);
                zMax = max(z);

                if zMin > 0
                    continue;
                end

                if zMax <= 0
                    nSub = nSub + 1;
                    idx = (nSub-1)*3 + (1:3);
                    subV(idx, :) = v;
                    subF(nSub, :) = idx;

                    e1 = v(2, :) - v(1, :);
                    e2 = v(3, :) - v(1, :);
                    n_vec = cross(e1, e2);
                    area = 0.5 * norm(n_vec);
                    wettedArea = wettedArea + area;
                else
                    [clippedPoly, wpEdge] = Algorithms.clipTriangleToWaterplane(v);

                    if ~isempty(wpEdge)
                        nNew = size(wpEdge, 1);
                        waterEdges(nWaterEdges + (1:nNew), :) = wpEdge;
                        nWaterEdges = nWaterEdges + nNew;
                    end

                    nClip = size(clippedPoly, 1);
                    if nClip >= 3
                        for t = 2:nClip-1
                            tri = [clippedPoly(1, :); clippedPoly(t, :); clippedPoly(t+1, :)];

                            nSub = nSub + 1;
                            idx = (nSub-1)*3 + (1:3);
                            subV(idx, :) = tri;
                            subF(nSub, :) = idx;

                            e1 = tri(2, :) - tri(1, :);
                            e2 = tri(3, :) - tri(1, :);
                            n_vec = cross(e1, e2);
                            area = 0.5 * norm(n_vec);
                            wettedArea = wettedArea + area;
                        end
                    end
                end
            end

            if nSub == 0
                COND1 = Algorithms.setEmptyHydrostatics(COND);
                submergedMesh = [];
                return;
            end

            subV = subV(1:nSub*3, :);
            subF = subF(1:nSub, :);
            waterEdges = waterEdges(1:nWaterEdges, :);

            % Vectorised divergence-theorem accumulation (no loop over faces)
            V0m = subV(subF(:,1), :);
            V1m = subV(subF(:,2), :);
            V2m = subV(subF(:,3), :);
            tetVols      = dot(V0m, cross(V1m, V2m, 2), 2) / 6;
            signedVolume = sum(tetVols);
            centroidAccum = sum(tetVols .* ((V0m + V1m + V2m) / 4), 1);

            Volume = abs(signedVolume);
            if abs(signedVolume) > 1e-12
                B_world = centroidAccum / signedVolume;
                B_ship = (R' * B_world')';
                B_ship(3) = B_ship(3) + TM;
            else
                B_ship = [NaN, NaN, NaN];
            end

            [Af, Lwl, Bwl, F_world, Ix, Iy] = Algorithms.computeWaterplaneProperties(waterEdges);

            if Af > 0
                F_ship = (R' * F_world')';
                F_ship(3) = F_ship(3) + TM;
            else
                F_ship = [NaN, NaN, NaN];
            end

            % Compute true midship cross-sectional area by slicing the
            % submerged mesh at the longitudinal midpoint of the original hull.
            xMid = (min(vertices(:,1)) + max(vertices(:,1))) / 2;
            Ams = Algorithms.computeMidshipArea(subV, subF, xMid);
            if Ams < 1e-12 && Lwl > 0
                % Fallback: midship slice returned nothing (e.g. hull does not
                % span x=0).  Use average section area; Cp/Cm will be inaccurate.
                Ams = Volume / Lwl;
            end

            COND1 = COND;
            COND1.V = Volume;
            COND1.B = B_ship;
            COND1.Af = Af;
            COND1.I = [Ix, Iy];
            COND1.Lwl = Lwl;
            COND1.Bwl = Bwl;
            COND1.F = F_ship;
            COND1.Ams = Ams;
            COND1.Ws = wettedArea;

            submergedMesh = struct();
            submergedMesh.vertices = subV;
            submergedMesh.faces = (1:(nSub*3))';
            submergedMesh.faces = reshape(submergedMesh.faces, 3, [])';
            submergedMesh.nVertices = nSub * 3;
            submergedMesh.nFaces = nSub;
        end

        function GZdata = computeGZcurve(SHIP, COND_ref, options)
            if nargin < 3 || isempty(options)
                options = struct();
            end

            if ~isfield(options, 'heelAngles'); options.heelAngles = 0:5:60; end
            if ~isfield(options, 'useMesh'); options.useMesh = true; end
            if ~isfield(options, 'tolVolume'); options.tolVolume = 1e-5; end
            if ~isfield(options, 'tolTrim'); options.tolTrim = 1e-4; end
            if ~isfield(options, 'parallel'); options.parallel = false; end
            if ~isfield(options, 'verbose'); options.verbose = true; end

            heelAngles = options.heelAngles;
            nAngles = numel(heelAngles);

            GZdata = struct('heel', num2cell(heelAngles), ...
                            'trim', num2cell(nan(1, nAngles)), ...
                            'draft', num2cell(nan(1, nAngles)), ...
                            'GZ', num2cell(nan(1, nAngles)), ...
                            'V', num2cell(nan(1, nAngles)), ...
                            'conv', num2cell(false(1, nAngles)));

            assettoOptions = struct();
            assettoOptions.useMesh = options.useMesh;
            assettoOptions.tolVolume = options.tolVolume;
            assettoOptions.tolTrim = options.tolTrim;
            assettoOptions.verbose = false;

            if options.verbose
                fprintf('\n========================================\n');
                fprintf('  Computing GZ Curve\n');
                fprintf('  %d heel angles: %.1f° to %.1f°\n', nAngles, min(heelAngles), max(heelAngles));
                fprintf('========================================\n\n');
                fprintf('%6s %8s %10s %10s %10s %8s\n', 'Heel', 'Trim', 'Draft', 'GZ', 'Volume', 'Status');
                fprintf('%s\n', repmat('-', 1, 60));
            end

            tGZ = tic;  % Use a handle so this does not interfere with outer timers
            if options.parallel && license('test', 'Distrib_Computing_Toolbox')
                % Parallel branch: collect all results first, then print in order.
                gdArr(1:nAngles) = deal(struct('heel',nan,'trim',nan,'draft',nan,'GZ',nan,'V',nan,'conv',false));
                parfor i = 1:nAngles
                    gdArr(i) = Algorithms.computeSingleHeel(SHIP, COND_ref, heelAngles(i), assettoOptions);
                end
                for i = 1:nAngles
                    GZdata(i) = gdArr(i);
                    if options.verbose
                        gd = GZdata(i);
                        status = 'FAIL'; if gd.conv; status = 'OK'; end
                        fprintf('%6.1f %8.3f %10.3f %10.4f %10.1f %8s\n', ...
                            gd.heel, gd.trim, gd.draft, gd.GZ, gd.V, status);
                    end
                end
            else
                % Serial branch (default)
                for i = 1:nAngles
                    GZdata(i) = Algorithms.computeSingleHeel(SHIP, COND_ref, heelAngles(i), assettoOptions);
                    if options.verbose
                        gd = GZdata(i);
                        status = 'FAIL'; if gd.conv; status = 'OK'; end
                        fprintf('%6.1f %8.3f %10.3f %10.4f %10.1f %8s\n', ...
                            gd.heel, gd.trim, gd.draft, gd.GZ, gd.V, status);
                    end
                end
            end
            elapsed = toc(tGZ);

            if options.verbose
                nConverged = sum([GZdata.conv]);
                fprintf('%s\n', repmat('-', 1, 60));
                fprintf('  Converged: %d/%d\n', nConverged, nAngles);
                fprintf('  Time: %.2f s (%.3f s per angle)\n', elapsed, elapsed/nAngles);

                GZvals = [GZdata.GZ];
                valid = ~isnan(GZvals);
                validIdx = find(valid);
                if ~isempty(validIdx)
                    [GZmax, iMax] = max(GZvals(validIdx));
                    heelMax = heelAngles(validIdx);
                    heelMax = heelMax(iMax);
                    iZero = find(GZvals(validIdx) < 0, 1, 'first');
                    if ~isempty(iZero) && iZero > 1
                        iPair = validIdx(iZero-1:iZero);
                        heelZero = interp1(GZvals(iPair), heelAngles(iPair), 0, 'linear', 'extrap');
                    else
                        heelZero = nan;
                    end
                    % Area under the GZ curve (stability energy index)
                    heelRad = heelAngles(validIdx) * pi / 180;
                    gzArea  = trapz(heelRad, GZvals(validIdx));

                    fprintf('\n  Stability Summary:\n');
                    fprintf('    Max GZ:     %.4f m at %.1f°\n', GZmax, heelMax);
                    if ~isnan(heelZero)
                        fprintf('    Range:      0° to %.1f°\n', heelZero);
                    end
                    fprintf('    GZ area:    %.4f m·rad\n', gzArea);
                end
            end
        end

        function COND1 = assetto(COND, SHIP, options)
            if nargin < 3 || isempty(options)
                options = struct();
            end

            if ~isfield(options, 'useMesh'); options.useMesh = true; end
            if ~isfield(options, 'tolVolume'); options.tolVolume = 1e-5; end
            if ~isfield(options, 'tolTrim'); options.tolTrim = 1e-4; end
            if ~isfield(options, 'maxIter'); options.maxIter = 100; end
            if ~isfield(options, 'relaxation'); options.relaxation = 0.8; end
            if ~isfield(options, 'verbose'); options.verbose = true; end

            Algorithms.validateInputs(COND, SHIP, options);

            V0 = COND.V;
            if abs(V0) < 1e-10
                error('Algorithms:zeroVolume', 'Target volume V0 = %.3e is too small.', V0);
            end

            if options.verbose
                fprintf('  Assetto(class) - Equilibrium Solver\n');
                fprintf('  Geometry: MESH | Heel: %.2f°\n', COND.tetaT);
            end

            geometryFunc = @(c) Algorithms.geometriaMesh(c, SHIP);

            COND1 = COND;
            COND1.tetaL = 0;
            COND1 = geometryFunc(COND1);

            history = struct();
            history.volume = [];
            history.trim = [];
            history.draft = [];

            [COND1, iterV0, convV0] = Algorithms.solveVolumeEquilibrium(COND1, V0, geometryFunc, options);
            history.volume = [history.volume; COND1.V];
            history.draft = [history.draft; COND1.T];

            if options.verbose
                fprintf('  Phase 1 (Volume):  %d iters, dV/V = %.2e\n', iterV0, convV0);
            end

            [COND1, iterTrim, convTrim, history] = Algorithms.solveTrimEquilibrium(COND1, V0, geometryFunc, options, history);
            if options.verbose
                fprintf('  Phase 2 (Trim):    %d iters, dθ = %.2e°\n', iterTrim, convTrim);
            end

            COND1 = geometryFunc(COND1);
            COND1 = Algorithms.computeRightingLever(COND1);

            if options.verbose
                fprintf('  Final: T=%.3f m, trim=%.4f°, GZ=%+.4f m\n', COND1.T, COND1.tetaL, COND1.GZ);
            end

            COND1.conv.iterV0 = iterV0;
            COND1.conv.iterTrim = iterTrim;
            COND1.conv.dVfinal = (COND1.V - V0) / V0;
            COND1.conv.dtlFinal = convTrim;
            COND1.conv.converged = abs(COND1.conv.dVfinal) < options.tolVolume && abs(COND1.conv.dtlFinal) < options.tolTrim;
            COND1.conv.history = history;
        end
    end

    methods (Static, Access = private)
        function gd = computeSingleHeel(SHIP, COND_ref, heelAngle, assettoOptions)
            gd = struct('heel', heelAngle, 'trim', nan, 'draft', nan, 'GZ', nan, 'V', nan, 'conv', false);
            try
                COND = COND_ref;
                COND.tetaT = heelAngle;
                COND.tetaL = 0;
                COND1 = Algorithms.assetto(COND, SHIP, assettoOptions);
                gd.trim = COND1.tetaL;
                gd.draft = COND1.T;
                gd.GZ = COND1.GZ;
                gd.V = COND1.V;
                gd.conv = COND1.conv.converged;
            catch
                gd.conv = false;
            end
        end

        function [clippedPoly, wpEdge] = clipTriangleToWaterplane(tri)
            % FIXED: Handle edge cases better
            clippedPolyBuf = zeros(4, 3);
            wpEdgeBuf = zeros(2, 3);
            nClip = 0;
            nWp = 0;
            n = size(tri, 1);

            % Tolerance for numerical comparisons
            tol = 1e-12;

            % Early-return: all vertices on the waterplane — whole triangle is
            % submerged and every edge is a waterplane edge.
            if all(abs(tri(:,3)) <= tol)
                clippedPoly = tri;
                wpEdge      = tri;
                return;
            end

            for i = 1:n
                iNext = mod(i, n) + 1;
                pCurr = tri(i, :);
                pNext = tri(iNext, :);
                zCurr = pCurr(3);
                zNext = pNext(3);

                % FIXED: Use tolerance for comparisons
                currInside = zCurr <= tol;
                nextInside = zNext <= tol;

                if currInside
                    nClip = nClip + 1;
                    clippedPolyBuf(nClip, :) = pCurr;
                end

                if currInside ~= nextInside
                    % Compute intersection with z=0 plane
                    if abs(zNext - zCurr) > tol
                        t = -zCurr / (zNext - zCurr);
                        t = max(0, min(1, t));  % Clamp to [0,1]
                        pInt = pCurr + t * (pNext - pCurr);
                        pInt(3) = 0;  % Force exactly zero
                        nClip = nClip + 1;
                        clippedPolyBuf(nClip, :) = pInt;
                        nWp = nWp + 1;
                        wpEdgeBuf(nWp, :) = pInt;
                    end
                end
            end

            clippedPoly = clippedPolyBuf(1:nClip, :);
            wpEdge = wpEdgeBuf(1:nWp, :);
        end

        function [Af, Lwl, Bwl, F_world, Ix, Iy] = computeWaterplaneProperties(waterEdges)
            % FIXED: Improved waterplane area calculation
            % Returns NaN for invalid cases instead of wrong values

            Af = 0; Lwl = 0; Bwl = 0; F_world = [NaN, NaN, NaN]; Ix = 0; Iy = 0;

            if size(waterEdges, 1) < 3
                return;
            end

            xy = waterEdges(:, 1:2);

            % Remove duplicate points
            [xy, ~, ~] = unique(xy, 'rows', 'stable');

            if size(xy, 1) < 3
                return;
            end

            % Compute centroid
            centroid = mean(xy, 1);

            % Sort by angle around centroid (for simple convex/concave polygons)
            angles = atan2(xy(:, 2) - centroid(2), xy(:, 1) - centroid(1));
            [~, sortIdx] = sort(angles);
            xySorted = xy(sortIdx, :);

            % Vectorised shoelace — build shifted-index arrays once
            n    = size(xySorted, 1);
            iN   = [2:n, 1]';                   % next-vertex indices
            x    = xySorted(:, 1);
            y    = xySorted(:, 2);
            xN   = x(iN);
            yN   = y(iN);
            cp   = x .* yN - xN .* y;           % cross-product per edge

            Af = abs(sum(cp)) / 2;

            if Af <= 1e-9
                Af = 0;
                return;
            end

            % Vectorised centroid
            cx = sum((x + xN) .* cp) / (6 * Af);
            cy = sum((y + yN) .* cp) / (6 * Af);
            F_world = [cx, cy, 0];

            % FIXED: Compute Lwl and Bwl from sorted polygon extents
            Lwl = max(x) - min(x);
            Bwl = max(y) - min(y);

            % Vectorised second moments of area (about origin, then corrected)
            Ix = abs(sum((y.^2 + y.*yN + yN.^2) .* cp)) / 12;
            Iy = abs(sum((x.^2 + x.*xN + xN.^2) .* cp)) / 12;

            % Parallel axis theorem to get moments about centroid
            Ix = Ix - Af * cy^2;
            Iy = Iy - Af * cx^2;

            % Ensure non-negative
            Ix = max(0, Ix);
            Iy = max(0, Iy);
        end

        function COND1 = setEmptyHydrostatics(COND)
            COND1 = COND;
            COND1.V = 0;
            COND1.B = [NaN, NaN, NaN];
            COND1.Af = 0;
            COND1.I = [0, 0];
            COND1.Lwl = 0;
            COND1.Bwl = 0;
            COND1.F = [NaN, NaN, NaN];
            COND1.Ams = 0;
            COND1.Ws = 0;
        end

        function [COND, iter, residual] = solveVolumeEquilibrium(COND, V0, geometryFunc, options)
            iter = 0;
            residual = inf;
            alpha = 1.0;

            % FIXED: Initialize residual properly
            if ~isnan(COND.V) && COND.V > 0
                residual = abs((V0 - COND.V) / V0);
            end

            while residual > options.tolVolume && iter < options.maxIter
                iter = iter + 1;
                dV = V0 - COND.V;

                if abs(V0) < 1e-12
                    break;
                end

                residual = abs(dV / V0);
                if residual <= options.tolVolume
                    break;
                end

                Af = max(COND.Af, 1e-6);
                dT = dV / Af;
                alpha = Algorithms.lineSearchVolume(COND, V0, dT, alpha, geometryFunc);
                COND.T = COND.T + alpha * dT;
                COND = geometryFunc(COND);

                if iter > 1
                    newResidual = abs((V0 - COND.V) / V0);
                    if newResidual > residual * 0.9
                        alpha = max(0.5, alpha * 0.8);
                    else
                        alpha = min(1.0, alpha * 1.1);
                    end
                end
            end
        end

        function [COND, iter, residual, history] = solveTrimEquilibrium(COND, V0, geometryFunc, options, history)
            iter = 0;
            residual = inf;
            relaxation = options.relaxation;

            % FIXED: Initialize residual properly
            dtl = Algorithms.computeTrimCorrection(COND);
            residual = abs(dtl);

            while residual > options.tolTrim && iter < options.maxIter
                iter = iter + 1;
                residual = abs(dtl);
                COND.tetaL = COND.tetaL + relaxation * dtl;
                COND = geometryFunc(COND);

                volumeTol = options.tolVolume * 10;
                volIter = 0;
                dV = V0 - COND.V;

                while abs(dV / V0) > volumeTol && volIter < 20
                    volIter = volIter + 1;
                    Af = max(COND.Af, 1e-6);
                    COND.T = COND.T + dV / Af;
                    COND = geometryFunc(COND);
                    dV = V0 - COND.V;
                end

                dtl = Algorithms.computeTrimCorrection(COND);

                if iter > 1
                    if abs(dtl) > residual * 0.95
                        relaxation = max(0.3, relaxation * 0.9);
                    else
                        relaxation = min(options.relaxation, relaxation * 1.05);
                    end
                end

                history.volume = [history.volume; COND.V];
                history.trim = [history.trim; COND.tetaL];
                history.draft = [history.draft; COND.T];
            end
        end

        function alpha = lineSearchVolume(COND, V0, dT, alpha0, geometryFunc)
            alpha = alpha0;
            maxTry = 5;
            COND_test = COND;
            COND_test.T = COND.T + alpha * dT;
            COND_test = geometryFunc(COND_test);

            if abs(V0) < 1e-12
                return;
            end

            residual_new = abs((V0 - COND_test.V) / V0);
            residual_old = abs((V0 - COND.V) / V0);
            tries = 0;

            while residual_new > residual_old * 1.5 && tries < maxTry && alpha > 0.1
                alpha = alpha * 0.5;
                COND_test.T = COND.T + alpha * dT;
                COND_test = geometryFunc(COND_test);
                residual_new = abs((V0 - COND_test.V) / V0);
                tries = tries + 1;
            end
        end

        function dtl = computeTrimCorrection(COND)
            % FIXED: Improved trim correction calculation
            % Bow-down (negative trim by stern) should give positive correction

            tetaTrad = COND.tetaT * pi / 180;
            tetaLrad = COND.tetaL * pi / 180;

            % Unit vector in longitudinal direction (in ship coordinates)
            % FIXED: Sign convention - positive trim by stern means bow is up
            % Rotation matrix for trim (pitch about y-axis, positive = bow up)
            % R_L = [cos(tetaL), 0, sin(tetaL); 0, 1, 0; -sin(tetaL), 0, cos(tetaL)]
            % After rotation, the longitudinal direction in world frame is:
            PI = [-cos(tetaTrad)*sin(tetaLrad), -sin(tetaTrad)*cos(tetaLrad), cos(tetaTrad)*cos(tetaLrad)];

            if COND.V > 1e-9 && COND.I(2) > 1e-9
                RL = COND.I(2) / COND.V;
                ML = COND.B + RL * PI;
            else
                dtl = 0;
                return;
            end

            BML = ML - COND.B;
            GML = ML - COND.G;
            normBML = norm(BML);
            normGML = norm(GML);

            if normBML < 1e-12 || normGML < 1e-12
                dtl = 0;
                return;
            end

            % Compute angle between BML and GML in the longitudinal plane
            % Positive correction means bow-down trimming moment needed
            dotProd = dot(BML, GML);
            crossProd = BML(1)*GML(3) - BML(3)*GML(1);

            % Angle in degrees
            dtl = atan2(crossProd, dotProd) * 180/pi;

            % Limit to reasonable range
            % Clamp to small angles to avoid atan2 instability
            % When angle approaches 180 deg, atan2 becomes very sensitive
            dtl = max(-15, min(15, dtl));
        end

        function COND = computeRightingLever(COND)
            % FIXED: Completely rewritten righting lever calculation
            % The GZ is the horizontal distance from G to the line of action of buoyancy
            %
            % In the rotated (world) coordinate system:
            % - Gravity acts in -Z direction
            % - Buoyancy acts in +Z direction (normal to waterplane)
            %
            % GZ = |GB| * sin(phi) where phi is the angle between GB and the vertical

            if any(isnan(COND.B)) || any(isnan(COND.G))
                COND.GZ = NaN;
                return;
            end

            % Heel and trim angles in radians
            tetaTrad = COND.tetaT * pi / 180;
            tetaLrad = COND.tetaL * pi / 180;

            % Build rotation matrix: ship coordinates -> world coordinates
            % R_T: roll about the x-axis (heel, positive = starboard down)
            R_T = [1, 0, 0;
                   0, cos(tetaTrad), -sin(tetaTrad);
                   0, sin(tetaTrad), cos(tetaTrad)];

            % R_L: pitch about the y-axis (trim, positive = bow up)
            R_L = [cos(tetaLrad), 0, sin(tetaLrad);
                   0, 1, 0;
                   -sin(tetaLrad), 0, cos(tetaLrad)];

            % Combined rotation: ship -> world
            R = R_L * R_T;

            % In world coordinates, vertical is Z direction
            % The righting lever is the horizontal component perpendicular to the
            % plane formed by the buoyancy force (vertical) and the BG vector
            %
            % Simplified: GZ = horizontal distance from G to vertical through B
            % In world frame, this is the magnitude of the horizontal components
            % of the vector from B to G

            % Actually, let's use the classical definition:
            % GZ = (KB + BM - KG) * sin(heel) for small angles
            % But for large angles, we need the cross product approach

            % Vector from G to B
            GB = COND.B - COND.G;

            % Unit vector in vertical direction (world frame)
            k_world = [0, 0, 1];

            % The righting moment arm GZ is the component perpendicular to gravity
            % GZ_vector = (GB . k) * k - GB  <- projection onto horizontal plane
            % But actually, we want the perpendicular distance from G to buoyancy line

            % In world frame, buoyancy acts upward (+Z), gravity acts downward (-Z)
            % The line of action of buoyancy passes through B and is vertical
            % GZ is the horizontal distance from G to this line

            % Transform G and B from ship coordinates to world coordinates.
            % In the world frame the waterplane is horizontal and gravity acts in -Z.
            % GZ is the perpendicular (horizontal) distance from G to the vertical
            % line of action of the buoyancy force, which passes through B_world.
            G_world = (R * COND.G')';
            B_world = (R * COND.B')';

            % In world frame, the vertical is Z
            % The horizontal plane is XY
            % GZ is the horizontal distance from G to the vertical line through B
            % Projected onto the direction perpendicular to both the ship's
            % longitudinal axis and the vertical

            % For a static GZ curve, we typically report the GZ in the
            % plane of symmetry (for symmetric hulls)
            % This is the y-component of the vector from the line of action
            % of buoyancy to G

            % Simpler: The righting lever is the moment arm
            % GZ = |GB x F_buoyancy_unit|
            % but projected onto the transverse direction

            % Let me use the most straightforward approach:
            % GZ = (ZB - ZG) * sin(phi) for small angles
            % But we need the exact value for large angles

            % The restoring moment is: M_r = rho * g * V * GZ
            % where GZ is the perpendicular distance from G to the line of action of buoyancy

            % Line of action of buoyancy in world frame: passes through B_world, direction [0,0,1]
            % Distance from G_world to this line:
            % d = | (G_world - B_world) x [0,0,1] |

            diff = G_world - B_world;
            cross_prod = cross(diff, [0, 0, 1]);

            % The magnitude of the cross product gives the perpendicular distance
            gz_magnitude = norm(cross_prod);

            % The sign of GZ: positive for righting moment (restoring)
            % For a positive heel angle (port side up), if G is to starboard of B,
            % we have a restoring moment
            %
            % In world frame after positive heel about X:
            % - Y positive is to port (original)
            % - Z positive is up
            %
            % For a restoring moment, G should be "above and to the side" of B
            % such that the moment brings the ship back

            % Sign determination:
            % If cross_prod(1) < 0, the moment is restoring for positive heel
            % Actually, let's look at the geometry:
            % At positive heel, the ship rotates so port side goes up
            % B moves to port (positive Y in world frame)
            % If G is above B's new position, there's a restoring moment

            % The sign convention: positive GZ = restoring (stability)
            % For positive heel: if G is to starboard of B's line of action, GZ > 0

            % Sign convention:
            %   cross([dx,dy,dz], [0,0,1]) = [dy, -dx, 0]
            %   For a positive heel (starboard down), B moves to the low side
            %   (negative Y in world frame) while G stays near Y=0 for a
            %   symmetric hull, so diff(2) = G_world(2)-B_world(2) > 0,
            %   which gives cross_prod(1) > 0 — a restoring (positive) GZ. ✓
            if cross_prod(1) > 0
                COND.GZ = gz_magnitude;  % Restoring
            else
                COND.GZ = -gz_magnitude; % Capsizing
            end

            % Alternative simpler formula for small angles:
            % GZ = GM * sin(heel), but this is approximate
        end

        function Ams = computeMidshipArea(subV, subF, xMid)
            % COMPUTEMIDSHIPAREA  Transverse cross-sectional area of the submerged
            % mesh at the longitudinal position x = xMid (default 0).
            %
            % Collects every triangle edge that crosses the plane x=xMid,
            % projects the intersection points onto the (y,z) plane, sorts them
            % by angle around their centroid, and applies the shoelace formula.
            %
            % Assumptions:
            %   - The hull is symmetric and centred at x = 0.
            %   - The cross-section is convex (angle-sort gives correct winding).
            %     For non-convex sections (catamarans, etc.) the result will
            %     overestimate area — flag as a known limitation.

            if nargin < 3 || isempty(xMid)
                xMid = 0;
            end

            tol = 1e-10;
            yzPoints = zeros(0, 2);

            nF = size(subF, 1);
            for f = 1:nF
                v   = subV(subF(f, :), :);   % 3×3 triangle vertices (world frame)
                xv  = v(:, 1) - xMid;        % x relative to slice plane

                for i = 1:3
                    j  = mod(i, 3) + 1;
                    xi = xv(i);
                    xj = xv(j);

                    if (xi < -tol && xj > tol) || (xi > tol && xj < -tol)
                        % Edge crosses the plane strictly — interpolate
                        t = -xi / (xj - xi);
                        t = max(0, min(1, t));
                        p = v(i, :) + t * (v(j, :) - v(i, :));
                        yzPoints(end + 1, :) = p(2:3); %#ok<AGROW>
                    elseif abs(xi) <= tol
                        % Vertex lies exactly on the slice plane
                        yzPoints(end + 1, :) = v(i, 2:3); %#ok<AGROW>
                    end
                end
            end

            if size(yzPoints, 1) < 3
                Ams = 0;
                return;
            end

            % De-duplicate points that are numerically identical
            yzRnd = round(yzPoints * 1e8) / 1e8;
            [~, ia] = unique(yzRnd, 'rows', 'stable');
            yzPoints = yzPoints(ia, :);

            if size(yzPoints, 1) < 3
                Ams = 0;
                return;
            end

            % Sort by angle around centroid (shoelace requires consistent winding)
            c      = mean(yzPoints, 1);
            angles = atan2(yzPoints(:, 2) - c(2), yzPoints(:, 1) - c(1));
            [~, idx]  = sort(angles);
            yz     = yzPoints(idx, :);

            % Vectorised shoelace
            n    = size(yz, 1);
            jIdx = [2:n, 1]';
            Ams  = abs(sum(yz(:,1) .* yz(jIdx,2) - yz(jIdx,1) .* yz(:,2))) / 2;
        end

        function validateInputs(COND, SHIP, options)
            required = {'V', 'tetaT', 'tetaL', 'T', 'G'};
            for k = 1:numel(required)
                if ~isfield(COND, required{k})
                    error('Algorithms:missingField', 'COND missing required field: .%s', required{k});
                end
            end
            if ~isstruct(SHIP)
                error('Algorithms:invalidSHIP', 'SHIP must be a struct.');
            end
            if options.useMesh
                if ~isfield(SHIP, 'vertices') || ~isfield(SHIP, 'faces')
                    error('Algorithms:invalidMesh', 'Mesh SHIP must have vertices and faces fields.');
                end
            else
                error('Algorithms:cellModeUnsupported', 'Class-only pipeline supports mesh mode only.');
            end
        end
    end
end
