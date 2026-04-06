function generate_box_barge()
%GENERATE_BOX_BARGE Create a watertight outward-oriented box barge mesh.
%
% Coordinate system:
%   x : stern (-L/2) to bow (+L/2)
%   y : port  (-B/2) to starboard (+B/2)
%   z : keel   0     to deck (+D)
%
% This convention matches the NAUTILUS hydrostatics workflow:
%   - z = 0 at keel
%   - positive z upward
%   - draft T measured upward from keel to waterline
%
% Output:
%   BoxBarge_mesh.obj saved in the same folder as this script.

    L = 10.0;   % Length [m]
    B = 4.0;    % Beam   [m]
    D = 2.0;    % Depth  [m]

    fprintf('Generating Box Barge hull...\n');
    fprintf('  L = %.2f m, B = %.2f m, D = %.2f m\n', L, B, D);
    fprintf('  Coordinate system: Z=0 at KEEL, Z=D at deck\n');

    % ---------------------------------------------------------------------
    % Vertex definition
    % ---------------------------------------------------------------------
    % 8 box corners
    %
    % Bottom (keel, z = 0)
    %   1: stern, port
    %   2: bow,   port
    %   3: bow,   starboard
    %   4: stern, starboard
    %
    % Top (deck, z = D)
    %   5: stern, port
    %   6: bow,   port
    %   7: bow,   starboard
    %   8: stern, starboard
    vertices = [
        -L/2, -B/2, 0;   % 1
         L/2, -B/2, 0;   % 2
         L/2,  B/2, 0;   % 3
        -L/2,  B/2, 0;   % 4
        -L/2, -B/2, D;   % 5
         L/2, -B/2, D;   % 6
         L/2,  B/2, D;   % 7
        -L/2,  B/2, D;   % 8
    ];

    % ---------------------------------------------------------------------
    % Face definition
    % ---------------------------------------------------------------------
    % 12 triangles, outward-oriented.
    %
    % Face normals:
    %   bottom    -> -z
    %   top       -> +z
    %   port      -> -y
    %   starboard -> +y
    %   stern     -> -x
    %   bow       -> +x
    faces = [
        % Bottom face (keel, z = 0), outward normal = -z
        1, 3, 2;
        1, 4, 3;

        % Top face (deck, z = D), outward normal = +z
        5, 6, 7;
        5, 7, 8;

        % Port side (y = -B/2), outward normal = -y
        1, 2, 6;
        1, 6, 5;

        % Starboard side (y = +B/2), outward normal = +y
        4, 7, 3;
        4, 8, 7;

        % Stern (x = -L/2), outward normal = -x
        1, 8, 4;
        1, 5, 8;

        % Bow (x = +L/2), outward normal = +x
        2, 7, 6;
        2, 3, 7;
    ];

    fprintf('  Generated %d vertices, %d faces\n', size(vertices, 1), size(faces, 1));

    % ---------------------------------------------------------------------
    % Mesh checks
    % ---------------------------------------------------------------------
    report = check_mesh(vertices, faces, L, B, D);

    if ~report.isWatertight
        warning('Mesh is not watertight.');
    end
    if report.nNonManifold > 0
        warning('Mesh has non-manifold edges.');
    end
    if report.signedVolume < 0
        warning('Mesh has negative signed volume; face winding is inward.');
    end

    % ---------------------------------------------------------------------
    % Write OBJ file
    % ---------------------------------------------------------------------
    output_file = fullfile(fileparts(mfilename('fullpath')), 'BoxBarge_mesh.obj');
    write_obj_file(output_file, vertices, faces, L, B, D);
    fprintf('  Saved: %s\n', output_file);

    % ---------------------------------------------------------------------
    % Print theoretical hydrostatic properties
    % ---------------------------------------------------------------------
    fprintf('\n  Theoretical hydrostatic properties:\n');
    for T = [0.5, 1.0, 1.5]
        V = L * B * T;
        KB = T / 2;
        fprintf('    T = %.1f m : V = %.2f m^3, KB = %.2f m\n', T, V, KB);
    end

    fprintf('\nDone.\n');
end


function report = check_mesh(vertices, faces, L, B, D)
%CHECK_MESH Run basic topological and geometric checks.

    report = struct();

    % ---------------------------------------------------------------------
    % Edge manifold / watertightness check
    % ---------------------------------------------------------------------
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

    fprintf('  Mesh check:\n');
    fprintf('    Unique edges       : %d\n', report.nEdges);
    fprintf('    Boundary edges     : %d\n', nBoundary);
    fprintf('    Non-manifold edges : %d\n', nNonManifold);

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

    % ---------------------------------------------------------------------
    % Signed volume check
    % ---------------------------------------------------------------------
    signedVol = compute_signed_volume(vertices, faces);
    absVol = abs(signedVol);
    expectedVol = L * B * D;
    relErrPct = 100 * abs(absVol - expectedVol) / max(expectedVol, eps);

    report.signedVolume = signedVol;
    report.volume = absVol;
    report.expectedVolume = expectedVol;
    report.volumeErrorPct = relErrPct;

    fprintf('    Signed volume      : %.6f m^3\n', signedVol);
    fprintf('    Absolute volume    : %.6f m^3\n', absVol);
    fprintf('    Expected volume    : %.6f m^3\n', expectedVol);
    fprintf('    Volume error       : %.6f %%\n', relErrPct);

    if signedVol > 0
        fprintf('    Face orientation   : outward\n');
    elseif signedVol < 0
        fprintf('    Face orientation   : inward\n');
    else
        fprintf('    Face orientation   : degenerate / zero-volume\n');
    end

    % ---------------------------------------------------------------------
    % Bounding box check
    % ---------------------------------------------------------------------
    xyzMin = min(vertices, [], 1);
    xyzMax = max(vertices, [], 1);
    span = xyzMax - xyzMin;

    fprintf('    Bounding box       : [%.3f, %.3f] x [%.3f, %.3f] x [%.3f, %.3f]\n', ...
        xyzMin(1), xyzMax(1), xyzMin(2), xyzMax(2), xyzMin(3), xyzMax(3));
    fprintf('    Span [L,B,D]       : [%.3f, %.3f, %.3f] m\n', span(1), span(2), span(3));

    if any(abs(span - [L, B, D]) > 1e-12)
        warning('Bounding box span does not match the prescribed dimensions exactly.');
    end
end


function vol = compute_signed_volume(vertices, faces)
%COMPUTE_SIGNED_VOLUME Compute signed volume using tetrahedra to the origin.

    v1 = vertices(faces(:,1), :);
    v2 = vertices(faces(:,2), :);
    v3 = vertices(faces(:,3), :);

    vol = sum(dot(v1, cross(v2, v3, 2), 2)) / 6;
end


function write_obj_file(filename, vertices, faces, L, B, D)
%WRITE_OBJ_FILE Export vertices and triangular faces to OBJ format.

    fid = fopen(filename, 'w');
    if fid < 0
        error('generate_box_barge:FileOpenFailed', 'Cannot open file: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid, '# Box Barge Mesh\n');
    fprintf(fid, '# Generated by NAUTILUS\n');
    fprintf(fid, '# L = %.6f m\n', L);
    fprintf(fid, '# B = %.6f m\n', B);
    fprintf(fid, '# D = %.6f m\n', D);
    fprintf(fid, '# Coordinate system: Z = 0 at keel, Z = D at deck\n');
    fprintf(fid, '# Outward-oriented triangular surface mesh\n');
    fprintf(fid, '#\n');

    for i = 1:size(vertices, 1)
        fprintf(fid, 'v %.6f %.6f %.6f\n', vertices(i, 1), vertices(i, 2), vertices(i, 3));
    end

    fprintf(fid, '\n');

    for i = 1:size(faces, 1)
        fprintf(fid, 'f %d %d %d\n', faces(i, 1), faces(i, 2), faces(i, 3));
    end
end