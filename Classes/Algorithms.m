classdef Algorithms
    % Algorithms
    % Core hydrostatics and stability algorithms for triangulated hull meshes.
    %
    % Main improvements:
    %   - centralised tolerances
    %   - safer clipping of partially submerged triangles
    %   - better defensive checks for invalid intermediate states
    %   - more robust equilibrium iteration logic
    %   - clearer righting lever computation and sign convention
    %   - more robust handling of degenerate waterplane and section cases
    %
    % Important limitation:
    %   Waterplane and midship-section reconstruction still rely on point-cloud
    %   ordering, which is acceptable for simple single-hull geometries but not
    %   fully reliable for disconnected or strongly non-convex sections.

    properties (Constant, Access = private)
        TOL_GEOM   = 1e-12;
        TOL_AREA   = 1e-10;
        TOL_VOLUME = 1e-12;
        TOL_ANGLE  = 1e-12;
        MAX_TRIM_STEP_DEG = 10;
        MIN_WP_AREA = 1e-8;
        MIN_VOL = 1e-10;
    end

    methods (Static)
        function [COND1, submergedMesh] = geometriaMesh(COND, mesh, varargin) %#ok<INUSD>
            validateattributes(COND.T,     {'numeric'}, {'scalar','real','finite'}, 'Algorithms.geometriaMesh', 'COND.T');
            validateattributes(COND.tetaT, {'numeric'}, {'scalar','real','finite'}, 'Algorithms.geometriaMesh', 'COND.tetaT');
            validateattributes(COND.tetaL, {'numeric'}, {'scalar','real','finite'}, 'Algorithms.geometriaMesh', 'COND.tetaL');

            if ~isstruct(mesh) || ~isfield(mesh, 'vertices') || ~isfield(mesh, 'faces')
                error('Algorithms:InvalidMesh', 'mesh must contain fields "vertices" and "faces".');
            end

            V = double(mesh.vertices);
            F = double(mesh.faces);

            if isempty(V) || isempty(F) || size(V,2) ~= 3 || size(F,2) ~= 3
                error('Algorithms:InvalidMeshShape', 'Mesh vertices/faces must be Nx3 / Mx3 arrays.');
            end

            nFaces = size(F, 1);

            phi   = deg2rad(COND.tetaT); % heel / roll
            theta = deg2rad(COND.tetaL); % trim / pitch
            Tm    = COND.T;

            R = Algorithms.buildRotation(theta, phi);

            % Ship frame -> world frame, then shift waterplane to z=0
            Vw = (R * V')';
            Vw(:,3) = Vw(:,3) - Tm;

            maxSubFaces = 3 * nFaces;
            subV = zeros(3 * maxSubFaces, 3);
            subF = zeros(maxSubFaces, 3);

            % Waterplane intersection points only; not edges-as-polygons
            wpPts = zeros(2 * nFaces, 3);

            nSub = 0;
            nWp  = 0;
            wettedArea = 0;

            for f = 1:nFaces
                tri = Vw(F(f,:), :);
                z = tri(:,3);

                if min(z) > Algorithms.TOL_GEOM
                    continue;
                end

                if max(z) <= Algorithms.TOL_GEOM
                    [subV, subF, nSub, wettedArea] = Algorithms.appendTriangle(subV, subF, nSub, wettedArea, tri);
                    continue;
                end

                [clippedPoly, wpLocal] = Algorithms.clipTriangleToWaterplane(tri);

                if ~isempty(wpLocal)
                    nNew = size(wpLocal,1);
                    wpPts(nWp + (1:nNew), :) = wpLocal;
                    nWp = nWp + nNew;
                end

                nClip = size(clippedPoly,1);
                if nClip >= 3
                    for k = 2:(nClip-1)
                        triSub = [clippedPoly(1,:); clippedPoly(k,:); clippedPoly(k+1,:)];
                        [subV, subF, nSub, wettedArea] = Algorithms.appendTriangle(subV, subF, nSub, wettedArea, triSub);
                    end
                end
            end

            if nSub == 0
                COND1 = Algorithms.setEmptyHydrostatics(COND);
                submergedMesh = [];
                return;
            end

            subV = subV(1:3*nSub, :);
            subF = subF(1:nSub, :);
            wpPts = wpPts(1:nWp, :);

            % Divergence-theorem accumulation on submerged triangulation
            V0m = subV(subF(:,1), :);
            V1m = subV(subF(:,2), :);
            V2m = subV(subF(:,3), :);

            tetVol = dot(V0m, cross(V1m, V2m, 2), 2) / 6;
            signedVolume = sum(tetVol);

            centroidAccum = sum(tetVol .* ((V0m + V1m + V2m) / 4), 1);

            Volume = abs(signedVolume);

            if abs(signedVolume) > Algorithms.TOL_VOLUME
                B_world = centroidAccum / signedVolume;
                B_ship = (R' * B_world')';
                B_ship(3) = B_ship(3) + Tm;
            else
                B_ship = [NaN, NaN, NaN];
            end

            [Af, Lwl, Bwl, F_world, Ix, Iy] = Algorithms.computeWaterplaneProperties(wpPts);

            if Af > 0 && all(isfinite(F_world))
                F_ship = (R' * F_world')';
                F_ship(3) = F_ship(3) + Tm;
            else
                F_ship = [NaN, NaN, NaN];
            end

            xMid = 0.5 * (min(V(:,1)) + max(V(:,1)));
            Ams = Algorithms.computeMidshipArea(subV, subF, xMid);

            if Ams < Algorithms.TOL_AREA && Lwl > Algorithms.TOL_GEOM
                Ams = Volume / Lwl;
            end

            COND1 = COND;
            COND1.V   = Volume;
            COND1.B   = B_ship;
            COND1.Af  = Af;
            COND1.I   = [Ix, Iy];
            COND1.Lwl = Lwl;
            COND1.Bwl = Bwl;
            COND1.F   = F_ship;
            COND1.Ams = Ams;
            COND1.Ws  = wettedArea;

            submergedMesh = struct();
            submergedMesh.vertices  = subV;
            submergedMesh.faces     = reshape(1:(3*nSub), 3, [])';
            submergedMesh.nVertices = size(subV,1);
            submergedMesh.nFaces    = size(subF,1);
        end

        function GZdata = computeGZcurve(SHIP, COND_ref, options)
            if nargin < 3 || isempty(options)
                options = struct();
            end

            options = Algorithms.applyDefaults(options, struct( ...
                'heelAngles', 0:5:60, ...
                'useMesh', true, ...
                'tolVolume', 1e-5, ...
                'tolTrim', 1e-4, ...
                'parallel', false, ...
                'verbose', true));

            heelAngles = options.heelAngles(:).';
            nAngles = numel(heelAngles);

            GZdata = repmat(struct( ...
                'heel', nan, 'trim', nan, 'draft', nan, ...
                'GZ', nan, 'V', nan, 'conv', false), 1, nAngles);

            assettoOptions = struct( ...
                'useMesh', options.useMesh, ...
                'tolVolume', options.tolVolume, ...
                'tolTrim', options.tolTrim, ...
                'verbose', false);

            if options.verbose
                fprintf('\n========================================\n');
                fprintf('  Computing GZ Curve\n');
                fprintf('  %d heel angles: %.1f° to %.1f°\n', nAngles, min(heelAngles), max(heelAngles));
                fprintf('========================================\n\n');
                fprintf('%6s %8s %10s %10s %10s %8s\n', 'Heel', 'Trim', 'Draft', 'GZ', 'Volume', 'Status');
                fprintf('%s\n', repmat('-', 1, 60));
            end

            t0 = tic;

            if options.parallel && license('test', 'Distrib_Computing_Toolbox')
                tmp(1:nAngles) = GZdata(1);
                parfor i = 1:nAngles
                    tmp(i) = Algorithms.computeSingleHeel(SHIP, COND_ref, heelAngles(i), assettoOptions);
                end
                GZdata = tmp;
            else
                for i = 1:nAngles
                    GZdata(i) = Algorithms.computeSingleHeel(SHIP, COND_ref, heelAngles(i), assettoOptions);
                end
            end

            elapsed = toc(t0);

            if options.verbose
                for i = 1:nAngles
                    gd = GZdata(i);
                    status = 'FAIL';
                    if gd.conv, status = 'OK'; end
                    fprintf('%6.1f %8.3f %10.3f %10.4f %10.1f %8s\n', ...
                        gd.heel, gd.trim, gd.draft, gd.GZ, gd.V, status);
                end

                fprintf('%s\n', repmat('-', 1, 60));
                fprintf('  Converged: %d/%d\n', sum([GZdata.conv]), nAngles);
                fprintf('  Time: %.2f s (%.3f s per angle)\n', elapsed, elapsed / max(nAngles,1));

                GZvals = [GZdata.GZ];
                valid = isfinite(GZvals);
                if any(valid)
                    hv = heelAngles(valid);
                    gv = GZvals(valid);

                    [GZmax, idxMax] = max(gv);
                    heelMax = hv(idxMax);

                    heelZero = nan;
                    idxSign = find(gv < 0, 1, 'first');
                    if ~isempty(idxSign) && idxSign > 1
                        pair = idxSign-1:idxSign;
                        heelZero = interp1(gv(pair), hv(pair), 0, 'linear', 'extrap');
                    end

                    gzArea = trapz(deg2rad(hv), gv);

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

            options = Algorithms.applyDefaults(options, struct( ...
                'useMesh', true, ...
                'tolVolume', 1e-5, ...
                'tolTrim', 1e-4, ...
                'maxIter', 100, ...
                'relaxation', 0.8, ...
                'verbose', true));

            Algorithms.validateInputs(COND, SHIP, options);

            V0 = COND.V;
            if ~isfinite(V0) || abs(V0) < Algorithms.MIN_VOL
                error('Algorithms:ZeroTargetVolume', 'Target displacement volume V0 is too small or invalid.');
            end

            if options.verbose
                fprintf('  Assetto(class) - Equilibrium Solver\n');
                fprintf('  Geometry: MESH | Heel: %.2f°\n', COND.tetaT);
            end

            geometryFunc = @(c) Algorithms.geometriaMesh(c, SHIP);

            COND1 = COND;
            COND1.tetaL = 0;
            COND1 = geometryFunc(COND1);

            history = struct('volume', [], 'trim', [], 'draft', []);

            [COND1, iterV0, convV0] = Algorithms.solveVolumeEquilibrium(COND1, V0, geometryFunc, options);
            history.volume(end+1,1) = COND1.V;
            history.draft(end+1,1)  = COND1.T;

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

            COND1.conv = struct();
            COND1.conv.iterV0     = iterV0;
            COND1.conv.iterTrim   = iterTrim;
            COND1.conv.dVfinal    = (COND1.V - V0) / V0;
            COND1.conv.dtlFinal   = convTrim;
            COND1.conv.converged  = abs(COND1.conv.dVfinal) < options.tolVolume && abs(COND1.conv.dtlFinal) < options.tolTrim;
            COND1.conv.history    = history;
        end
    end

    methods (Static, Access = private)
        function R = buildRotation(theta, phi)
            % theta = trim about y
            % phi   = heel about x
            R_T = [1, 0, 0;
                   0, cos(phi), -sin(phi);
                   0, sin(phi),  cos(phi)];

            R_L = [ cos(theta), 0, sin(theta);
                    0,          1, 0;
                   -sin(theta), 0, cos(theta)];

            R = R_L * R_T;
        end

        function [subV, subF, nSub, wettedArea] = appendTriangle(subV, subF, nSub, wettedArea, tri)
            e1 = tri(2,:) - tri(1,:);
            e2 = tri(3,:) - tri(1,:);
            area = 0.5 * norm(cross(e1, e2));

            if area <= Algorithms.TOL_AREA
                return;
            end

            nSub = nSub + 1;
            idx = (nSub-1)*3 + (1:3);
            subV(idx,:) = tri;
            subF(nSub,:) = idx;
            wettedArea = wettedArea + area;
        end

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
                gd.conv = isfield(COND1, 'conv') && isfield(COND1.conv, 'converged') && COND1.conv.converged;
            catch
                gd.conv = false;
            end
        end

        function [clippedPoly, wpPoints] = clipTriangleToWaterplane(tri)
            tol = Algorithms.TOL_GEOM;

            clippedBuf = zeros(4,3);
            wpBuf = zeros(2,3);
            nClip = 0;
            nWp = 0;

            if all(abs(tri(:,3)) <= tol)
                clippedPoly = tri;
                wpPoints = tri(:, :);
                return;
            end

            for i = 1:3
                j = mod(i,3) + 1;

                p1 = tri(i,:);
                p2 = tri(j,:);
                z1 = p1(3);
                z2 = p2(3);

                in1 = z1 <= tol;
                in2 = z2 <= tol;

                if in1
                    nClip = nClip + 1;
                    clippedBuf(nClip,:) = p1;
                end

                if xor(in1, in2)
                    dz = z2 - z1;
                    if abs(dz) > tol
                        t = -z1 / dz;
                        t = max(0, min(1, t));
                        pint = p1 + t * (p2 - p1);
                        pint(3) = 0;

                        nClip = nClip + 1;
                        clippedBuf(nClip,:) = pint;

                        nWp = nWp + 1;
                        wpBuf(nWp,:) = pint;
                    end
                elseif abs(z1) <= tol && abs(z2) <= tol
                    % whole edge on waterplane; keep endpoints as candidate points
                    if nWp + 2 <= size(wpBuf,1)
                        nWp = nWp + 1;
                        wpBuf(nWp,:) = p1;
                        if norm(p2 - p1) > tol
                            nWp = nWp + 1;
                            wpBuf(nWp,:) = p2;
                        end
                    end
                end
            end

            clippedPoly = clippedBuf(1:nClip,:);
            wpPoints = wpBuf(1:nWp,:);
        end

        function [Af, Lwl, Bwl, F_world, Ix, Iy] = computeWaterplaneProperties(wpPts)
            Af = 0;
            Lwl = 0;
            Bwl = 0;
            F_world = [NaN, NaN, NaN];
            Ix = 0;
            Iy = 0;

            if size(wpPts,1) < 3
                return;
            end

            xy = wpPts(:,1:2);
            xy = Algorithms.uniqueTolRows(xy, 1e-9);

            if size(xy,1) < 3
                return;
            end

            % Convex hull is safer than centroid-angle sorting on raw noisy points.
            % It may under-represent concave geometries, but it avoids self-crossing
            % polygons and catastrophic area errors.
            try
                k = convhull(xy(:,1), xy(:,2));
            catch
                return;
            end

            if numel(k) < 4
                return;
            end

            poly = xy(k(1:end-1), :);
            x = poly(:,1);
            y = poly(:,2);

            n = numel(x);
            j = [2:n, 1]';
            xj = x(j);
            yj = y(j);

            cp = x .* yj - xj .* y;
            A_signed = 0.5 * sum(cp);
            Af = abs(A_signed);

            if Af < Algorithms.MIN_WP_AREA
                Af = 0;
                return;
            end

            denom = 6 * A_signed;
            if abs(denom) < Algorithms.TOL_GEOM
                Af = 0;
                return;
            end

            cx = sum((x + xj) .* cp) / denom;
            cy = sum((y + yj) .* cp) / denom;

            F_world = [cx, cy, 0];
            Lwl = max(x) - min(x);
            Bwl = max(y) - min(y);

            Ix0 = sum((y.^2 + y.*yj + yj.^2) .* cp) / 12;
            Iy0 = sum((x.^2 + x.*xj + xj.^2) .* cp) / 12;

            Ix = abs(Ix0 - A_signed * cy^2);
            Iy = abs(Iy0 - A_signed * cx^2);

            Ix = max(Ix, 0);
            Iy = max(Iy, 0);
        end

        function COND1 = setEmptyHydrostatics(COND)
            COND1 = COND;
            COND1.V   = 0;
            COND1.B   = [NaN, NaN, NaN];
            COND1.Af  = 0;
            COND1.I   = [0, 0];
            COND1.Lwl = 0;
            COND1.Bwl = 0;
            COND1.F   = [NaN, NaN, NaN];
            COND1.Ams = 0;
            COND1.Ws  = 0;
        end

        function [COND, iter, residual] = solveVolumeEquilibrium(COND, V0, geometryFunc, options)
            iter = 0;
            residual = inf;
            alpha = 1.0;

            if isfinite(COND.V) && COND.V > Algorithms.MIN_VOL
                residual = abs((V0 - COND.V) / V0);
            end

            while residual > options.tolVolume && iter < options.maxIter
                iter = iter + 1;

                if ~isfinite(COND.V) || COND.Af < Algorithms.MIN_WP_AREA
                    break;
                end

                dV = V0 - COND.V;
                residual = abs(dV / V0);

                if residual <= options.tolVolume
                    break;
                end

                dT = dV / max(COND.Af, Algorithms.MIN_WP_AREA);

                alpha = Algorithms.lineSearchVolume(COND, V0, dT, alpha, geometryFunc);

                Tnew = COND.T + alpha * dT;
                if ~isfinite(Tnew)
                    break;
                end

                COND.T = Tnew;
                COND = geometryFunc(COND);

                if isfinite(COND.V)
                    newResidual = abs((V0 - COND.V) / V0);
                    if newResidual > residual * 0.95
                        alpha = max(0.2, 0.8 * alpha);
                    else
                        alpha = min(1.0, 1.05 * alpha);
                    end
                    residual = newResidual;
                else
                    break;
                end
            end
        end

        function [COND, iter, residual, history] = solveTrimEquilibrium(COND, V0, geometryFunc, options, history)
            iter = 0;
            dtl = Algorithms.computeTrimCorrection(COND);
            residual = abs(dtl);
            relaxation = options.relaxation;

            while residual > options.tolTrim && iter < options.maxIter
                iter = iter + 1;

                if ~isfinite(dtl)
                    break;
                end

                oldResidual = residual;

                step = relaxation * dtl;
                step = max(-Algorithms.MAX_TRIM_STEP_DEG, min(Algorithms.MAX_TRIM_STEP_DEG, step));

                COND.tetaL = COND.tetaL + step;
                COND = geometryFunc(COND);

                if ~isfinite(COND.V) || COND.Af < Algorithms.MIN_WP_AREA
                    break;
                end

                volIter = 0;
                volumeTol = 10 * options.tolVolume;
                dV = V0 - COND.V;

                while abs(dV / V0) > volumeTol && volIter < 20
                    volIter = volIter + 1;

                    dT = dV / max(COND.Af, Algorithms.MIN_WP_AREA);
                    COND.T = COND.T + dT;
                    COND = geometryFunc(COND);

                    if ~isfinite(COND.V) || COND.Af < Algorithms.MIN_WP_AREA
                        break;
                    end

                    dV = V0 - COND.V;
                end

                dtl = Algorithms.computeTrimCorrection(COND);
                residual = abs(dtl);

                if iter > 1
                    if residual > 0.95 * oldResidual
                        relaxation = max(0.25, 0.85 * relaxation);
                    else
                        relaxation = min(options.relaxation, 1.02 * relaxation);
                    end
                end

                history.volume(end+1,1) = COND.V;
                history.trim(end+1,1)   = COND.tetaL;
                history.draft(end+1,1)  = COND.T;
            end
        end

        function alpha = lineSearchVolume(COND, V0, dT, alpha0, geometryFunc)
            alpha = alpha0;
            maxTry = 5;

            if ~isfinite(dT) || abs(V0) < Algorithms.MIN_VOL
                return;
            end

            residualOld = abs((V0 - COND.V) / V0);

            for k = 1:maxTry
                test = COND;
                test.T = COND.T + alpha * dT;
                test = geometryFunc(test);

                if ~isfinite(test.V)
                    alpha = 0.5 * alpha;
                    continue;
                end

                residualNew = abs((V0 - test.V) / V0);
                if residualNew <= 1.5 * residualOld || alpha <= 0.1
                    return;
                end

                alpha = 0.5 * alpha;
            end
        end

        function dtl = computeTrimCorrection(COND)
            if ~isfield(COND, 'B') || ~isfield(COND, 'G') || ~isfield(COND, 'I') || ~isfield(COND, 'V')
                dtl = 0;
                return;
            end

            if any(~isfinite(COND.B)) || any(~isfinite(COND.G)) || ~isfinite(COND.V) || COND.V <= Algorithms.MIN_VOL
                dtl = 0;
                return;
            end

            if numel(COND.I) < 2 || ~isfinite(COND.I(2)) || COND.I(2) <= Algorithms.TOL_GEOM
                dtl = 0;
                return;
            end

            phi   = deg2rad(COND.tetaT);
            theta = deg2rad(COND.tetaL);

            % Approximate longitudinal metacentre direction in ship frame/world-aligned logic
            PI = [-cos(phi) * sin(theta), ...
                  -sin(phi) * cos(theta), ...
                   cos(phi) * cos(theta)];

            RL = COND.I(2) / COND.V;
            ML = COND.B + RL * PI;

            BML = ML - COND.B;
            GML = ML - COND.G;

            n1 = norm(BML);
            n2 = norm(GML);

            if n1 < Algorithms.TOL_GEOM || n2 < Algorithms.TOL_GEOM
                dtl = 0;
                return;
            end

            dotProd = dot(BML, GML);
            crossXZ = BML(1) * GML(3) - BML(3) * GML(1);

            dtl = rad2deg(atan2(crossXZ, dotProd));
            dtl = max(-Algorithms.MAX_TRIM_STEP_DEG, min(Algorithms.MAX_TRIM_STEP_DEG, dtl));
        end

        function COND = computeRightingLever(COND)
            if any(~isfinite(COND.B)) || any(~isfinite(COND.G))
                COND.GZ = NaN;
                return;
            end

            phi   = deg2rad(COND.tetaT);
            theta = deg2rad(COND.tetaL);
            R = Algorithms.buildRotation(theta, phi);

            G_world = (R * COND.G(:)).';
            B_world = (R * COND.B(:)).';

            % Distance from G to vertical buoyancy line through B
            d = G_world - B_world;
            gzMag = hypot(d(1), d(2));

            if gzMag < Algorithms.TOL_GEOM
                COND.GZ = 0;
                return;
            end

            % Sign convention:
            % positive GZ should correspond to restoring moment for positive heel.
            % For roll about x with gravity in -z, the restoring moment sign is
            % governed by the x-component of r x F, where r is from G to B and F is +z.
            r = B_world - G_world;
            m = cross(r, [0 0 1]);
            mx = m(1);

            if COND.tetaT >= 0
                sgn = sign(mx);
            else
                sgn = -sign(mx);
            end

            if sgn == 0
                sgn = 1;
            end

            COND.GZ = sgn * gzMag;
        end

        function Ams = computeMidshipArea(subV, subF, xMid)
            if nargin < 3 || isempty(xMid)
                xMid = 0;
            end

            tol = 1e-10;
            yzPoints = zeros(0,2);

            nF = size(subF,1);
            for f = 1:nF
                tri = subV(subF(f,:), :);
                xv = tri(:,1) - xMid;

                for i = 1:3
                    j = mod(i,3) + 1;
                    xi = xv(i);
                    xj = xv(j);

                    if (xi < -tol && xj > tol) || (xi > tol && xj < -tol)
                        t = -xi / (xj - xi);
                        t = max(0, min(1, t));
                        p = tri(i,:) + t * (tri(j,:) - tri(i,:));
                        yzPoints(end+1,:) = p(2:3); %#ok<AGROW>
                    elseif abs(xi) <= tol
                        yzPoints(end+1,:) = tri(i,2:3); %#ok<AGROW>
                    end
                end
            end

            if size(yzPoints,1) < 3
                Ams = 0;
                return;
            end

            yzPoints = Algorithms.uniqueTolRows(yzPoints, 1e-8);
            if size(yzPoints,1) < 3
                Ams = 0;
                return;
            end

            % Use convex hull for robustness against self-crossing sort orders
            try
                k = convhull(yzPoints(:,1), yzPoints(:,2));
            catch
                Ams = 0;
                return;
            end

            if numel(k) < 4
                Ams = 0;
                return;
            end

            yz = yzPoints(k(1:end-1), :);
            x = yz(:,1);
            y = yz(:,2);

            n = numel(x);
            j = [2:n, 1]';
            Ams = 0.5 * abs(sum(x .* y(j) - x(j) .* y));
        end

        function validateInputs(COND, SHIP, options)
            required = {'V','tetaT','tetaL','T','G'};
            for k = 1:numel(required)
                if ~isfield(COND, required{k})
                    error('Algorithms:MissingField', 'COND is missing required field ".%s".', required{k});
                end
            end

            if ~isstruct(SHIP)
                error('Algorithms:InvalidSHIP', 'SHIP must be a struct.');
            end

            if ~isfield(options, 'useMesh') || ~options.useMesh
                error('Algorithms:CellModeUnsupported', 'This pipeline currently supports mesh mode only.');
            end

            if ~isfield(SHIP, 'vertices') || ~isfield(SHIP, 'faces')
                error('Algorithms:InvalidMesh', 'SHIP must contain "vertices" and "faces".');
            end

            if ~isnumeric(COND.G) || numel(COND.G) ~= 3 || any(~isfinite(COND.G))
                error('Algorithms:InvalidG', 'COND.G must be a finite 1x3 or 3x1 numeric vector.');
            end
        end

        function S = applyDefaults(S, defaults)
            names = fieldnames(defaults);
            for i = 1:numel(names)
                if ~isfield(S, names{i})
                    S.(names{i}) = defaults.(names{i});
                end
            end
        end

        function A = uniqueTolRows(A, tol)
            if isempty(A)
                return;
            end
            Ar = round(A / tol) * tol;
            [~, ia] = unique(Ar, 'rows', 'stable');
            A = A(ia, :);
        end
    end
end