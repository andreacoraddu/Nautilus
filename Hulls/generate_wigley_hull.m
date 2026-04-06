function generate_wigley_hull()
%GENERATE_WIGLEY_HULL Create a watertight Wigley-hull benchmark mesh.
%
% Coordinate system:
%   x : stern (-L/2) to bow (+L/2)
%   y : port  (-)     to starboard (+)
%   z : keel   0      to waterline (+D)
%
% Hull definition:
%   y = +/- (B/2) * (1 - (2x/L)^2) * (1 - (z/D)^2)
%
% where:
%   x in [-L/2, L/2]
%   z in [0, D]
%
% Theoretical displaced volume at draft T:
%   V(T) = (2/3) * L * B * T * (1 - (1/3)*(T/D)^2)
%
% Output:
%   Wigley_hull_mesh.obj saved in the same folder as this script.

    hullDir = fileparts(mfilename('fullpath'));
    outputFile = fullfile(hullDir, 'Wigley_hull_mesh.obj');

    % ---------------------------------------------------------------------
    % Hull parameters (must match Input_Wigley.m)
    % ---------------------------------------------------------------------
    L = 4.0;      % Length [m]
    B = 0.4;      % Beam   [m]
    D = 0.35;     % Depth  [m]

    % Mesh resolution
    nx = 40;      % Longitudinal divisions
    nz = 20;      % Vertical divisions

    fprintf('Generating Wigley hull...\n');
    fprintf('  L = %.2f m, B = %.2f m, D = %.2f m\n', L, B, D);
    fprintf('  Coordinate system: Z=0 at KEEL, Z=D at waterline\n');

    % ---------------------------------------------------------------------
    % Parametric grid
    % ---------------------------------------------------------------------
    % u : longitudinal parameter in [-1, 1]
    % w : vertical parameter in [0, 1]
    u = linspace(-1, 1, nx + 1);
    w = linspace(0, 1, nz + 1);

    % ---------------------------------------------------------------------
    % Vertex generation
    % ---------------------------------------------------------------------
    % Shared vertices are used whenever the half-breadth collapses to the
    % centerplane, so that the final surface remains manifold.
    %
    % vertexMap(j,i,1) = port-side vertex index
    % vertexMap(j,i,2) = starboard-side vertex index
    vertices = zeros(0, 3);
    vertexMap = zeros(nz + 1, nx + 1, 2);

    vidx = 0;
    tolCenter = 1e-12;

    for j = 1:(nz + 1)
        for i = 1:(nx + 1)
            u_i = u(i);
            w_j = w(j);

            x_i = 0.5 * L * u_i;
            z_j = D * w_j;

            % Half-breadth magnitude from the prescribed Wigley definition
            yMag = 0.5 * B * (1 - u_i^2) * (1 - w_j^2);

            if yMag < tolCenter
                % Shared centerplane vertex
                vidx = vidx + 1;
                vertices(vidx, :) = [x_i, 0.0, z_j];
                vertexMap(j, i, 1) = vidx;
                vertexMap(j, i, 2) = vidx;
            else
                % Port side: y < 0
                vidx = vidx + 1;
                vertices(vidx, :) = [x_i, -yMag, z_j];
                vertexMap(j, i, 1) = vidx;

                % Starboard side: y > 0
                vidx = vidx + 1;
                vertices(vidx, :) = [x_i, +yMag, z_j];
                vertexMap(j, i, 2) = vidx;
            end
        end
    end

    % ---------------------------------------------------------------------
    % Face generation
    % ---------------------------------------------------------------------
    faces = zeros(0, 3);

    % ---- Port side ------------------------------------------------------
    % Port is y < 0. We first create a consistent triangulation, then we
    % globally correct orientation later using the signed volume check.
    for j = 1:nz
        for i = 1:nx
            v1 = vertexMap(j,   i,   1);
            v2 = vertexMap(j,   i+1, 1);
            v3 = vertexMap(j+1, i+1, 1);
            v4 = vertexMap(j+1, i,   1);

            % Skip fully degenerate cells
            if numel(unique([v1 v2 v3 v4])) < 3
                continue;
            end

            % Avoid pathological diagonal choice where the patch collapses
            % near the centerplane.
            onCenter = abs(vertices([v1, v3, v4], 2)) < tolCenter;
            if all(onCenter) && abs(vertices(v2, 2)) >= tolCenter
                faces = [faces; v1, v4, v2]; %#ok<AGROW>
                faces = [faces; v2, v4, v3]; %#ok<AGROW>
            else
                faces = [faces; v1, v2, v3]; %#ok<AGROW>
                faces = [faces; v1, v3, v4]; %#ok<AGROW>
            end
        end
    end

    % ---- Starboard side -------------------------------------------------
    for j = 1:nz
        for i = 1:nx
            v1 = vertexMap(j,   i,   2);
            v2 = vertexMap(j,   i+1, 2);
            v3 = vertexMap(j+1, i+1, 2);
            v4 = vertexMap(j+1, i,   2);

            if numel(unique([v1 v2 v3 v4])) < 3
                continue;
            end

            onCenter = abs(vertices([v1, v3, v4], 2)) < tolCenter;
            if all(onCenter) && abs(vertices(v2, 2)) >= tolCenter
                faces = [faces; v1, v2, v4]; %#ok<AGROW>
                faces = [faces; v2, v3, v4]; %#ok<AGROW>
            else
                faces = [faces; v1, v3, v2]; %#ok<AGROW>
                faces = [faces; v1, v4, v3]; %#ok<AGROW>
            end
        end
    end

    % ---- Bottom closure at z = 0 ---------------------------------------
    keelRow = 1;
    for i = 1:nx
        portCurrent = vertexMap(keelRow, i,   1);
        portNext    = vertexMap(keelRow, i+1, 1);
        stbdCurrent = vertexMap(keelRow, i,   2);
        stbdNext    = vertexMap(keelRow, i+1, 2);

        if portCurrent == stbdCurrent && portNext ~= stbdNext
            faces = [faces; portCurrent, portNext, stbdNext]; %#ok<AGROW>
        elseif portCurrent ~= stbdCurrent && portNext == stbdNext
            faces = [faces; portCurrent, portNext, stbdCurrent]; %#ok<AGROW>
        elseif portCurrent ~= stbdCurrent && portNext ~= stbdNext
            faces = [faces; portCurrent, portNext, stbdNext]; %#ok<AGROW>
            faces = [faces; portCurrent, stbdNext, stbdCurrent]; %#ok<AGROW>
        end
    end

    fprintf('  Generated %d vertices, %d faces\n', size(vertices, 1), size(faces, 1));

    % ---------------------------------------------------------------------
    % Remove any accidental degenerate faces
    % ---------------------------------------------------------------------
    faces = remove_degenerate_faces(vertices, faces, 1e-14);
    fprintf('  After cleanup: %d faces\n', size(faces, 1));

    % ---------------------------------------------------------------------
    % Mesh checks before orientation correction
    % ---------------------------------------------------------------------
    report = check_mesh(vertices, faces);

    % If the surface is closed but inward-oriented, flip all triangles.
    if report.isWatertight && report.isManifold && report.signedVolume < 0
        fprintf('  Flipping global face orientation to enforce outward normals...\n');
        faces = faces(:, [1 3 2]);
        report = check_mesh(vertices, faces);
    end

    % ---------------------------------------------------------------------
    % Theoretical reference volume
    % ---------------------------------------------------------------------
    T_design = 0.25;
    V_theo = (2/3) * L * B * T_design * (1 - (1/3) * (T_design / D)^2);
    fprintf('  Theoretical volume at T = %.2f m: %.6f m^3\n', T_design, V_theo);

    % ---------------------------------------------------------------------
    % Write OBJ
    % ---------------------------------------------------------------------
    write_obj_file(outputFile, vertices, faces, L, B, D);
    fprintf('  Saved: %s\n', outputFile);
    fprintf('Done.\n');
end


function faces = remove_degenerate_faces(vertices, faces, areaTol)
%REMOVE_DEGENERATE_FACES Remove repeated-index and near-zero-area triangles.

    if isempty(faces)
        return;
    end

    % Repeated vertex indices
    badIdx = faces(:,1) == faces(:,2) | ...
             faces(:,1) == faces(:,3) | ...
             faces(:,2) == faces(:,3);

    faces = faces(~badIdx, :);

    if isempty(faces)
        return;
    end

    v1 = vertices(faces(:,1), :);
    v2 = vertices(faces(:,2), :);
    v3 = vertices(faces(:,3), :);

    area2 = vecnorm(cross(v2 - v1, v3 - v1, 2), 2, 2);  % twice area
    good = area2 > areaTol;

    nRemoved = nnz(~good);
    if nRemoved > 0
        fprintf('  Removed %d degenerate faces.\n', nRemoved);
    end

    faces = faces(good, :);
end


function report = check_mesh(vertices, faces)
%CHECK_MESH Run topological and signed-volume checks.

    report = struct();

    edges = [
        faces(:, [1 2]);
        faces(:, [2 3]);
        faces(:, [3 1])
    ];

    edgesSorted = sort(edges, 2);
    [uniqueEdges, ~, edgeIds] = unique(edgesSorted, 'rows');
    edgeCount = accumarray(edgeIds, 1);

    nBoundary = nnz(edgeCount == 1);
    nNonManifold = nnz(edgeCount > 2);

    report.nEdges = size(uniqueEdges, 1);
    report.nBoundary = nBoundary;
    report.nNonManifold = nNonManifold;
    report.isWatertight = (nBoundary == 0);
    report.isManifold = (nNonManifold == 0);

    signedVol = compute_signed_volume(vertices, faces);
    report.signedVolume = signedVol;
    report.volume = abs(signedVol);

    fprintf('  Mesh check:\n');
    fprintf('    Unique edges       : %d\n', report.nEdges);
    fprintf('    Boundary edges     : %d\n', report.nBoundary);
    fprintf('    Non-manifold edges : %d\n', report.nNonManifold);

    if report.isWatertight
        fprintf('    Watertight         : YES\n');
    else
        fprintf('    Watertight         : NO\n');
    end

    if report.isManifold
        fprintf('    Manifold           : YES\n');
    else
        fprintf('    Manifold           : NO\n');
    end

    fprintf('    Signed volume      : %.6f m^3\n', report.signedVolume);
    fprintf('    Absolute volume    : %.6f m^3\n', report.volume);

    if report.signedVolume > 0
        fprintf('    Face orientation   : outward\n');
    elseif report.signedVolume < 0
        fprintf('    Face orientation   : inward\n');
    else
        fprintf('    Face orientation   : degenerate / zero-volume\n');
    end
end


function vol = compute_signed_volume(vertices, faces)
%COMPUTE_SIGNED_VOLUME Compute signed enclosed volume.

    v1 = vertices(faces(:,1), :);
    v2 = vertices(faces(:,2), :);
    v3 = vertices(faces(:,3), :);

    vol = sum(dot(v1, cross(v2, v3, 2), 2)) / 6;
end


function write_obj_file(filename, vertices, faces, L, B, D)
%WRITE_OBJ_FILE Export the triangulated mesh to OBJ format.

    fid = fopen(filename, 'w');
    if fid < 0
        error('generate_wigley_hull:FileOpenFailed', 'Cannot open file: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid, '# Wigley Hull Mesh\n');
    fprintf(fid, '# Generated by NAUTILUS\n');
    fprintf(fid, '# L = %.6f m\n', L);
    fprintf(fid, '# B = %.6f m\n', B);
    fprintf(fid, '# D = %.6f m\n', D);
    fprintf(fid, '# Coordinate system: Z = 0 at keel, Z = D at waterline\n');
    fprintf(fid, '# Outward-oriented triangular surface mesh\n');
    fprintf(fid, '#\n');

    for i = 1:size(vertices, 1)
        fprintf(fid, 'v %.8f %.8f %.8f\n', vertices(i,1), vertices(i,2), vertices(i,3));
    end

    fprintf(fid, '\n');

    for i = 1:size(faces, 1)
        fprintf(fid, 'f %d %d %d\n', faces(i,1), faces(i,2), faces(i,3));
    end
end