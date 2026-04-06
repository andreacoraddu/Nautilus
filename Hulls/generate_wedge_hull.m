function generate_wedge_hull()
%GENERATE_WEDGE_HULL Create a watertight outward-oriented wedge-prism mesh.
%
% Coordinate system:
%   x : stern (-L/2) to bow (+L/2)
%   y : port  (-B/2) to starboard (+B/2)
%   z : keel   0     to deck (+D)
%
% Cross-section in the y-z plane is triangular:
%   keel point        : (y = 0,      z = 0)
%   deck port edge    : (y = -B/2,   z = D)
%   deck stbd edge    : (y = +B/2,   z = D)
%
% The section is extruded uniformly along x over length L.
%
% Output:
%   Wedge_hull_mesh.obj saved in the same folder as this script.

    L = 8.0;
    B = 3.0;
    D = 2.0;

    fprintf('Generating Wedge hull...\n');
    fprintf('  L = %.2f m, B = %.2f m, D = %.2f m\n', L, B, D);
    fprintf('  Coordinate system: Z=0 at KEEL, Z=D at deck\n');

    x0 = -L/2;
    x1 =  L/2;

    % ---------------------------------------------------------------------
    % Vertex definition
    % ---------------------------------------------------------------------
    % Stern section (x = x0)
    ks = [x0,  0.0,   0.0];   % 1: stern keel
    ps = [x0, -B/2,   D  ];   % 2: stern port deck edge
    ss = [x0,  B/2,   D  ];   % 3: stern starboard deck edge

    % Bow section (x = x1)
    kb = [x1,  0.0,   0.0];   % 4: bow keel
    pb = [x1, -B/2,   D  ];   % 5: bow port deck edge
    sb = [x1,  B/2,   D  ];   % 6: bow starboard deck edge

    vertices = [ks; ps; ss; kb; pb; sb];

    % ---------------------------------------------------------------------
    % Face definition
    % ---------------------------------------------------------------------
    % 8 triangles, outward-oriented.
    %
    % Face normals:
    %   stern cap   -> -x
    %   bow cap     -> +x
    %   port face   -> outward to port / down-left
    %   stbd face   -> outward to starboard / down-right
    %   deck face   -> +z
    faces = [
        % Stern triangular cap (outward normal = -x)
        1, 3, 2;

        % Bow triangular cap (outward normal = +x)
        4, 5, 6;

        % Port sloped face
        1, 2, 5;
        1, 5, 4;

        % Starboard sloped face
        1, 4, 6;
        1, 6, 3;

        % Deck face (z = D, outward normal = +z)
        2, 3, 6;
        2, 6, 5;
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
    output_file = fullfile(fileparts(mfilename('fullpath')), 'Wedge_hull_mesh.obj');
    write_obj_file(output_file, vertices, faces, L, B, D);
    fprintf('  Saved: %s\n', output_file);

    % ---------------------------------------------------------------------
    % Print theoretical hydrostatic properties
    % ---------------------------------------------------------------------
    fprintf('\n  Theoretical hydrostatic properties:\n');
    for T = [0.5, 1.0, 1.5]
        % Width at the waterline z = T
        Bwl = B * (T / D);

        % Submerged cross-sectional area of similar triangle
        Asec = B * T^2 / (2 * D);

        % Displacement volume
        V = L * Asec;

        % Vertical centre of buoyancy above keel
        KB = (2/3) * T;

        % Waterplane area
        Af = L * Bwl;

        fprintf('    T = %.1f m : V = %.4f m^3, KB = %.4f m, Af = %.4f m^2, Bwl = %.4f m\n', ...
            T, V, KB, Af, Bwl);
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
    expectedVol = L * B * D / 2;
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

    if abs(span(1) - L) > 1e-12 || abs(span(2) - B) > 1e-12 || abs(span(3) - D) > 1e-12
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
        error('generate_wedge_hull:FileOpenFailed', 'Cannot open file: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid, '# Wedge Hull Mesh\n');
    fprintf(fid, '# Generated by NAUTILUS\n');
    fprintf(fid, '# L = %.6f m\n', L);
    fprintf(fid, '# B = %.6f m\n', B);
    fprintf(fid, '# D = %.6f m\n', D);
    fprintf(fid, '# Coordinate system: Z = 0 at keel, Z = D at deck\n');
    fprintf(fid, '# Outward-oriented triangular surface mesh\n');
    fprintf(fid, '#\n');

    for i = 1:size(vertices, 1)
        fprintf(fid, 'v %.10f %.10f %.10f\n', vertices(i,1), vertices(i,2), vertices(i,3));
    end

    fprintf(fid, '\n');

    for i = 1:size(faces, 1)
        fprintf(fid, 'f %d %d %d\n', faces(i,1), faces(i,2), faces(i,3));
    end
end