function generate_wigley_hull()
%GENERATE_WIGLEY_HULL Creates a mathematically accurate Wigley hull mesh
% Coordinate system: Z=0 at KEEL, Z=D at waterline (matching original code)
%
% The Wigley hull is defined by:
%   y = +/- (B/2) * (1 - (2x/L)^2) * (1 - (z/D)^2)
% where x is longitudinal (-L/2 to L/2), z is vertical (0 to D)
%
% Theoretical volume at draft T: V = (2/3) * L * B * T * (1 - (1/3)*(T/D)^2)

    hullDir = fileparts(mfilename('fullpath'));
    output_file = fullfile(hullDir, 'Wigley_hull_mesh.obj');

    % Hull parameters (must match Input_Wigley.m)
    L = 4.0;        % Length [m]
    B = 0.4;        % Beam [m]
    D = 0.35;       % Depth [m]

    % Mesh resolution
    nx = 40;        % Longitudinal divisions
    nz = 20;        % Vertical divisions

    fprintf('Generating Wigley hull...\n');
    fprintf('  L = %.2f m, B = %.2f m, D = %.2f m\n', L, B, D);
    fprintf('  Coordinate system: Z=0 at KEEL, Z=D at waterline\n');

    % Create vertices on the hull surface
    vertices = [];

    % Parametric coordinates
    u = linspace(-1, 1, nx+1);  % -1 to 1 (stern to bow)
    w = linspace(0, 1, nz+1);   % 0 to 1 (keel to waterline)

    % Create vertices for port and starboard sides.
    % When the half-breadth collapses to the centerline, both sides share
    % the same vertex so the stem/stern/waterline seams are manifold.
    vertex_map = zeros(nz+1, nx+1, 2);  % (z, x, side) -> vertex index
    % side: 1 = port (y>0), 2 = starboard (y<0)

    vidx = 0;
    for j = 1:nz+1  % z direction (keel to waterline)
        for i = 1:nx+1  % x direction
            u_i = u(i);
            w_j = w(j);

            x_i = (L/2) * u_i;
            z_j = D * w_j;

            % Wigley hull equation: y = (B/2) * (1 - u^2) * (1 - w^2)
            y_mag = (B/2) * (1 - u_i^2) * (1 - w_j^2);

            % Share exact centerline vertices so the generated mesh stays
            % watertight after triangulation.
            if y_mag < 1e-10
                vidx = vidx + 1;
                vertices(vidx, :) = [x_i, 0, z_j];
                vertex_map(j, i, 1) = vidx;
                vertex_map(j, i, 2) = vidx;
            else
                % Port side (y >= 0)
                vidx = vidx + 1;
                vertices(vidx, :) = [x_i, y_mag, z_j];
                vertex_map(j, i, 1) = vidx;

                % Starboard side (y <= 0)
                vidx = vidx + 1;
                vertices(vidx, :) = [x_i, -y_mag, z_j];
                vertex_map(j, i, 2) = vidx;
            end
        end
    end

    % Create faces
    faces = [];

    % Port side faces
    for j = 1:nz
        for i = 1:nx
            v1 = vertex_map(j, i, 1);
            v2 = vertex_map(j, i+1, 1);
            v3 = vertex_map(j+1, i+1, 1);
            v4 = vertex_map(j+1, i, 1);

            % Skip degenerate faces
            if v1 == v2 || v2 == v3 || v3 == v4 || v4 == v1
                continue;
            end

            onCenterline = abs(vertices([v1, v3, v4], 2)) < 1e-12;
            if all(onCenterline) && abs(vertices(v2, 2)) >= 1e-12
                % The stern-top corner collapses onto the centerplane.
                % Use the opposite diagonal so port/stbd patches do not
                % overlap on the same interior edge.
                faces = [faces; v1, v4, v2];
                faces = [faces; v2, v4, v3];
            else
                % Two triangles, outward normal (pointing to port/y+)
                faces = [faces; v1, v3, v2];
                faces = [faces; v1, v4, v3];
            end
        end
    end

    % Starboard side faces (reverse winding for outward normal pointing to starboard/y-)
    for j = 1:nz
        for i = 1:nx
            v1 = vertex_map(j, i, 2);
            v2 = vertex_map(j, i+1, 2);
            v3 = vertex_map(j+1, i+1, 2);
            v4 = vertex_map(j+1, i, 2);

            % Skip degenerate faces
            if v1 == v2 || v2 == v3 || v3 == v4 || v4 == v1
                continue;
            end

            onCenterline = abs(vertices([v1, v3, v4], 2)) < 1e-12;
            if all(onCenterline) && abs(vertices(v2, 2)) >= 1e-12
                faces = [faces; v1, v2, v4];
                faces = [faces; v2, v3, v4];
            else
                % Two triangles, reversed winding for outward normal
                faces = [faces; v1, v2, v3];
                faces = [faces; v1, v3, v4];
            end
        end
    end

    % Bottom face at z=0
    keel_z_idx = 1;
    for i = 1:nx
        port_current = vertex_map(keel_z_idx, i, 1);
        port_next = vertex_map(keel_z_idx, i+1, 1);
        stbd_current = vertex_map(keel_z_idx, i, 2);
        stbd_next = vertex_map(keel_z_idx, i+1, 2);

        % The bottom strip collapses to a triangle at the pointed bow/stern.
        % Outward normal points in -z direction (down).
        if port_current == stbd_current && port_next ~= stbd_next
            faces = [faces; port_current, port_next, stbd_next];
        elseif port_current ~= stbd_current && port_next == stbd_next
            faces = [faces; port_current, port_next, stbd_current];
        elseif port_current ~= stbd_current && port_next ~= stbd_next
            faces = [faces; port_current, port_next, stbd_next];
            faces = [faces; port_current, stbd_next, stbd_current];
        end
    end

    fprintf('  Generated %d vertices, %d faces\n', size(vertices, 1), size(faces, 1));

    % Check mesh
    check_mesh(vertices, faces);

    % Calculate theoretical volume at design draft
    T_design = 0.25;
    V_theo = (2/3) * L * B * T_design * (1 - (1/3)*(T_design/D)^2);
    fprintf('  Theoretical volume at T=%.2f m: %.4f m³\n', T_design, V_theo);

    % Write OBJ file
    write_obj_file(output_file, vertices, faces, L, B, D);
    fprintf('  Saved: %s\n', output_file);
end

function check_mesh(vertices, faces)
    % Check for boundary edges
    edges = [faces(:,[1 2]); faces(:,[2 3]); faces(:,[3 1])];
    edgesSorted = sort(edges, 2);
    [~, ~, edgeIds] = unique(edgesSorted, 'rows');
    edgeCount = accumarray(edgeIds, 1);
    nBoundary = nnz(edgeCount == 1);
    nNonManifold = nnz(edgeCount > 2);

    fprintf('  Mesh check: %d boundary edges, %d non-manifold edges\n', nBoundary, nNonManifold);
    if nBoundary == 0
        fprintf('  Mesh is watertight.\n');
    else
        fprintf('  Warning: Mesh is not watertight!\n');
    end

    % Calculate volume
    vol = compute_volume(vertices, faces);
    fprintf('  Mesh volume: %.4f m³\n', vol);
end

function vol = compute_volume(vertices, faces)
    v1 = vertices(faces(:,1), :);
    v2 = vertices(faces(:,2), :);
    v3 = vertices(faces(:,3), :);
    vol = sum(dot(v1, cross(v2, v3, 2), 2)) / 6;
    vol = abs(vol);
end

function write_obj_file(filename, vertices, faces, L, B, D)
    fid = fopen(filename, 'w');
    if fid < 0
        error('Cannot open file: %s', filename);
    end

    fprintf(fid, '# Wigley Hull Mesh\n');
    fprintf(fid, '# Generated by NAUTILUS\n');
    fprintf(fid, '# L=%.2f, B=%.2f, D=%.2f\n', L, B, D);
    fprintf(fid, '# Coordinate system: Z=0 at KEEL, Z=D at waterline\n');
    fprintf(fid, '#\n');

    for i = 1:size(vertices, 1)
        fprintf(fid, 'v %.6f %.6f %.6f\n', vertices(i, 1), vertices(i, 2), vertices(i, 3));
    end

    for i = 1:size(faces, 1)
        fprintf(fid, 'f %d %d %d\n', faces(i, 1), faces(i, 2), faces(i, 3));
    end

    fclose(fid);
end
