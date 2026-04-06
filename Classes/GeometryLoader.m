classdef GeometryLoader
    %GeometryLoader Load hull input and mesh geometry from multiple formats.
    %
    % FIXED BUGS:
    % 1. readObjMesh - Now handles 0-based and negative indices properly
    % 2. readStlMesh - Better error handling for binary STL files
    % 3. loadInput - Better error messages for missing VCG

    properties
        autoSaveMesh (1,1) logical = true
        meshOptions struct = struct()
    end

    methods
        function obj = GeometryLoader(cfg)
            if nargin < 1
                cfg = struct();
            end
            if isfield(cfg, 'autoSaveMesh')
                obj.autoSaveMesh = logical(cfg.autoSaveMesh);
            end
            if isfield(cfg, 'meshOptions')
                obj.meshOptions = cfg.meshOptions;
            end
        end

        function [mesh, inputData, baseName] = loadHullMesh(obj, hullName)
            baseName = regexprep(hullName, '_mesh$', '');
            inputScript = GeometryLoader.resolveInputScript(baseName);

            inputData = obj.loadInput(inputScript, baseName);

            [meshFile, meshFormat] = GeometryLoader.resolveMeshFile(baseName);

            switch meshFormat
                case 'mat'
                    meshData = load(meshFile, 'mesh');
                    if ~isfield(meshData, 'mesh')
                        error('GeometryLoader:InvalidMeshMat', ...
                            'MAT mesh file %s must contain a variable named mesh.', meshFile);
                    end
                    mesh = meshData.mesh;

                case 'obj'
                    mesh = GeometryLoader.readObjMesh(meshFile);

                case 'stl'
                    mesh = GeometryLoader.readStlMesh(meshFile);

                otherwise
                    error('GeometryLoader:UnsupportedFormat', ...
                        'Unsupported mesh format %s for file %s.', meshFormat, meshFile);
            end

            mesh = GeometryLoader.normalizeMesh(mesh);
        end
    end

    methods (Static, Access = private)
        function inputData = loadInput(inputScript, baseName)
            inputData = struct();
            
            % Run the input script
            try
                run(fullfile('Hulls', [inputScript, '.m']));
            catch ME
                error('GeometryLoader:InputScriptError', ...
                    'Failed to run input script %s: %s', inputScript, ME.message);
            end

            if ~isfield(inputData, 'VCG')
                if exist('INPUT', 'var') && isfield(INPUT, 'VCG')
                    inputData = INPUT;
                else
                    error('GeometryLoader:MissingVCG', ...
                        'Input file for %s must define INPUT.VCG or inputData.VCG.', baseName);
                end
            end
        end

        function inputScript = resolveInputScript(baseName)
            candidates = { ...
                ['Input_', baseName], ...
                ['Input_', regexprep(baseName, '_hull$', '')], ...
                ['Input_', regexprep(baseName, '_mesh$', '')], ...
                ['Input_', regexprep(regexprep(baseName, '_mesh$', ''), '_hull$', '')] ...
            };

            found = '';
            for k = 1:numel(candidates)
                c = candidates{k};
                if exist(fullfile('Hulls', [c, '.m']), 'file')
                    found = c;
                    break;
                end
            end

            if isempty(found)
                error('GeometryLoader:MissingInput', ...
                    'No input script found for hull %s.', baseName);
            end

            inputScript = found;
        end

        function [meshFile, meshFormat] = resolveMeshFile(baseName)
            candidates = { ...
                struct('path', fullfile('Hulls', [baseName, '_mesh.obj']), 'fmt', 'obj'), ...
                struct('path', fullfile('Hulls', [baseName, '_mesh.stl']), 'fmt', 'stl'), ...
                struct('path', fullfile('Hulls', [baseName, '.obj']), 'fmt', 'obj'), ...
                struct('path', fullfile('Hulls', [baseName, '.stl']), 'fmt', 'stl'), ...
                struct('path', fullfile('Hulls', [baseName, '_mesh.mat']), 'fmt', 'mat') ...
            };

            for k = 1:numel(candidates)
                c = candidates{k};
                if exist(c.path, 'file')
                    meshFile = c.path;
                    meshFormat = c.fmt;
                    return;
                end
            end

            error('GeometryLoader:MissingGeometry', ...
                ['Mesh file not found for %s. Expected one of: %s, %s, %s, %s, %s. ' ...
                'Section .mat hull files are not used in this workflow.'], ...
                baseName, ...
                fullfile('Hulls', [baseName, '_mesh.mat']), ...
                fullfile('Hulls', [baseName, '_mesh.obj']), ...
                fullfile('Hulls', [baseName, '_mesh.stl']), ...
                fullfile('Hulls', [baseName, '.obj']), ...
                fullfile('Hulls', [baseName, '.stl']));
        end

        function mesh = normalizeMesh(mesh)
            if ~isstruct(mesh) || ~isfield(mesh, 'vertices') || ~isfield(mesh, 'faces')
                error('GeometryLoader:InvalidMesh', ...
                    'Loaded mesh must contain vertices and faces.');
            end

            mesh.vertices = double(mesh.vertices);
            mesh.faces = double(mesh.faces);

            if size(mesh.vertices, 2) ~= 3 || size(mesh.faces, 2) ~= 3
                error('GeometryLoader:InvalidMeshShape', ...
                    'Mesh vertices and faces must be Nx3 and Mx3 arrays.');
            end

            mesh.nVertices = size(mesh.vertices, 1);
            mesh.nFaces = size(mesh.faces, 1);
        end

        function mesh = readObjMesh(objPath)
            fid = fopen(objPath, 'r');
            if fid < 0
                error('GeometryLoader:ObjOpenFailed', ...
                    'Cannot open OBJ file: %s', objPath);
            end

            C = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
            fclose(fid);

            lines = C{1};
            vertices = zeros(0, 3);
            faces = zeros(0, 3);
            rawFaces = {};
            sawZeroBasedIndex = false;
            
            % Track normals and textures (for error reporting)
            nNormals = 0;
            nTexCoords = 0;

            for i = 1:numel(lines)
                line = strtrim(lines{i});
                if isempty(line) || startsWith(line, '#')
                    continue;
                end

                if startsWith(line, 'v ')
                    vals = sscanf(line(2:end), '%f %f %f');
                    if numel(vals) >= 3
                        vertices(end+1, :) = vals(1:3)'; %#ok<AGROW>
                    end
                elseif startsWith(line, 'vn ')
                    nNormals = nNormals + 1;
                elseif startsWith(line, 'vt ')
                    nTexCoords = nTexCoords + 1;
                elseif startsWith(line, 'f ')
                    parts = strsplit(strtrim(line(2:end)));
                    if numel(parts) < 3
                        continue;
                    end

                    idx = zeros(numel(parts), 1);
                    for p = 1:numel(parts)
                        token = parts{p};
                        slashPos = strfind(token, '/');
                        if ~isempty(slashPos)
                            token = token(1:slashPos(1)-1);
                        end
                        idx(p) = str2double(token);
                    end
                    sawZeroBasedIndex = sawZeroBasedIndex || any(idx == 0);
                    rawFaces{end+1} = idx; %#ok<AGROW>
                end
            end

            nV = size(vertices, 1);
            for i = 1:numel(rawFaces)
                idx = rawFaces{i};

                % FIXED: Handle OBJ index conventions
                % - 1-based: standard positive index
                % - 0-based: convert the full non-negative file to 1-based
                % - Negative: relative to end of vertex list
                if sawZeroBasedIndex
                    nonNegativeMask = idx >= 0;
                    idx(nonNegativeMask) = idx(nonNegativeMask) + 1;
                end
                for p = 1:numel(idx)
                    if idx(p) < 0
                        % Negative index - relative to end
                        idx(p) = nV + idx(p) + 1;
                    end
                end

                % Remove invalid indices
                validIdx = idx(~isnan(idx) & idx >= 1 & idx <= nV);

                if numel(validIdx) < 3
                    continue;
                end

                % Triangulate polygonal face as fan around first vertex.
                for p = 2:(numel(validIdx)-1)
                    faces(end+1, :) = [validIdx(1), validIdx(p), validIdx(p+1)]; %#ok<AGROW>
                end
            end
            
            % Report parsing summary
            if nNormals > 0 || nTexCoords > 0
                fprintf('  [OBJ Parser] Vertices: %d, Faces: %d, Normals: %d, TexCoords: %d\n', ...
                    size(vertices, 1), size(faces, 1), nNormals, nTexCoords);
            end

            mesh = struct('vertices', vertices, 'faces', faces);
        end

        function mesh = readStlMesh(stlPath)
            % FIXED: Handle both modern MATLAB built-in stlread (R2019b+, returns
            % a triangulation object) and legacy FEX stlread (returns [faces, vertices]).

            stlreadAvailable = exist('stlread', 'file') == 2 || ...
                               exist('stlread', 'builtin') == 5;

            if stlreadAvailable
                try
                    out = stlread(stlPath);
                    if isa(out, 'triangulation')
                        % R2019b+ built-in: returns a triangulation object
                        mesh = struct('vertices', double(out.Points), ...
                                      'faces',    double(out.ConnectivityList));
                    else
                        % Legacy FEX stlread: first output is faces matrix.
                        % Call again with two outputs to also get vertices.
                        [F, V] = stlread(stlPath);
                        mesh = struct('vertices', double(V), 'faces', double(F));
                    end
                    return;
                catch ME
                    warning('GeometryLoader:StlReadFailed', ...
                        'stlread failed: %s. Falling back to ASCII parser.', ME.message);
                end
            end

            % Try to detect binary vs ASCII
            fid = fopen(stlPath, 'r');
            if fid < 0
                error('GeometryLoader:StlOpenFailed', ...
                    'Cannot open STL file: %s', stlPath);
            end
            
            % Read first 80 bytes (header)
            header = fread(fid, 80, '*char')';
            fseek(fid, 0, 'eof');
            fileSize = ftell(fid);
            fclose(fid);
            
            % Check if it starts with "solid" (likely ASCII)
            isAscii = ~isempty(header) && startsWith(strtrim(header), 'solid');
            
            if ~isAscii
                % Try to verify binary format
                % Binary STL: 80-byte header + 4-byte triangle count + triangles
                if fileSize >= 84
                    fid = fopen(stlPath, 'r');
                    fseek(fid, 80, 'bof');
                    nTri = fread(fid, 1, 'uint32');
                    fclose(fid);
                    
                    expectedSize = 80 + 4 + nTri * 50;  % 50 bytes per triangle
                    if abs(expectedSize - fileSize) <= 2  % Allow small padding difference
                        error('GeometryLoader:BinarySTL', ...
                            'File %s appears to be binary STL. Please install stlread or convert to ASCII.', stlPath);
                    end
                end
            end

            % Parse as ASCII
            fid = fopen(stlPath, 'r');
            C = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
            fclose(fid);

            lines = C{1};
            verts = zeros(0, 3);

            for i = 1:numel(lines)
                line = strtrim(lines{i});
                if startsWith(line, 'vertex ')
                    vals = sscanf(line(7:end), '%f %f %f');
                    if numel(vals) == 3
                        verts(end+1, :) = vals'; %#ok<AGROW>
                    end
                end
            end

            if isempty(verts) || mod(size(verts,1), 3) ~= 0
                error('GeometryLoader:StlParseFailed', ...
                    ['ASCII STL parse failed for %s. ' ...
                    'File may be binary or corrupted.'], stlPath);
            end

            nTri = size(verts,1) / 3;
            [V, ~, ic] = unique(verts, 'rows', 'stable');
            F = reshape(ic, 3, nTri)';
            mesh = struct('vertices', V, 'faces', F);
        end
    end
end
