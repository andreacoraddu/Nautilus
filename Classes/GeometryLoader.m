classdef GeometryLoader
    % GeometryLoader
    % Load hull input data and triangulated mesh geometry from multiple formats.
    %
    % Supported mesh formats:
    %   - OBJ
    %   - STL
    %   - MAT (must contain variable 'mesh' with fields vertices and faces)
    %
    % Main improvements over the previous version:
    %   1) Safer OBJ parsing:
    %      - standard 1-based OBJ indices are handled correctly
    %      - negative relative indices are handled correctly
    %      - index 0 is treated as invalid by default
    %      - optional nonstandard 0-based OBJ support is explicit, not guessed
    %
    %   2) More robust input-script loading:
    %      - input script is run via absolute path
    %      - both INPUT and inputData conventions are supported
    %      - clearer validation of required fields
    %
    %   3) Better mesh normalization and validation:
    %      - finite numeric checks
    %      - integer face validation
    %      - duplicate/degenerate face removal
    %
    %   4) Better STL loading:
    %      - supports built-in triangulation-returning stlread
    %      - supports legacy [F,V] stlread
    %      - binary STL is detected more reliably
    %
    % Notes:
    %   - This loader assumes triangular faces in the normalized mesh.
    %   - Polygonal OBJ faces are triangulated by fan decomposition.
    %   - This class does not attempt to repair severely broken topology.
    
    properties
        autoSaveMesh (1,1) logical = true
        meshOptions struct = struct()

        % If true, allow nonstandard OBJ files with 0-based positive indices.
        % Default is false because standard OBJ uses 1-based indexing.
        allowZeroBasedObj (1,1) logical = false

        % If true, remove repeated and degenerate triangles during normalization.
        cleanMesh (1,1) logical = true

        % Tolerance used in mesh sanity checks.
        geometricTolerance (1,1) double = 1e-12
    end

    methods
        function obj = GeometryLoader(cfg)
            if nargin < 1 || isempty(cfg)
                cfg = struct();
            end

            if isfield(cfg, 'autoSaveMesh')
                obj.autoSaveMesh = logical(cfg.autoSaveMesh);
            end
            if isfield(cfg, 'meshOptions')
                obj.meshOptions = cfg.meshOptions;
            end
            if isfield(cfg, 'allowZeroBasedObj')
                obj.allowZeroBasedObj = logical(cfg.allowZeroBasedObj);
            end
            if isfield(cfg, 'cleanMesh')
                obj.cleanMesh = logical(cfg.cleanMesh);
            end
            if isfield(cfg, 'geometricTolerance')
                obj.geometricTolerance = double(cfg.geometricTolerance);
            end
        end

        function [mesh, inputData, baseName] = loadHullMesh(obj, hullName)
            arguments
                obj
                hullName (1,:) char
            end

            baseName = regexprep(strtrim(hullName), '_mesh$', '');

            repoRoot = GeometryLoader.getRepoRoot();
            hullsDir  = fullfile(repoRoot, 'Hulls');

            inputScript = GeometryLoader.resolveInputScript(baseName, hullsDir);
            inputData   = GeometryLoader.loadInput(inputScript, baseName, hullsDir);

            [meshFile, meshFormat] = GeometryLoader.resolveMeshFile(baseName, hullsDir);

            switch lower(meshFormat)
                case 'mat'
                    S = load(meshFile, 'mesh');
                    if ~isfield(S, 'mesh')
                        error('GeometryLoader:InvalidMeshMat', ...
                            'MAT mesh file "%s" must contain a variable named "mesh".', meshFile);
                    end
                    mesh = S.mesh;

                case 'obj'
                    mesh = GeometryLoader.readObjMesh(meshFile, obj.allowZeroBasedObj);

                case 'stl'
                    mesh = GeometryLoader.readStlMesh(meshFile);

                otherwise
                    error('GeometryLoader:UnsupportedFormat', ...
                        'Unsupported mesh format "%s" for file "%s".', meshFormat, meshFile);
            end

            mesh = GeometryLoader.normalizeMesh(mesh, obj.cleanMesh, obj.geometricTolerance);

            % Store source metadata for traceability
            mesh.sourceFile   = meshFile;
            mesh.sourceFormat = meshFormat;
            mesh.baseName     = baseName;
        end
    end

    methods (Static, Access = private)
        function repoRoot = getRepoRoot()
            % Try to locate the repository root from this class location.
            thisFile = mfilename('fullpath');
            classDir = fileparts(thisFile);
            repoRoot = fileparts(classDir);
        end

        function inputData = loadInput(inputScript, baseName, hullsDir)
            inputPath = fullfile(hullsDir, [inputScript, '.m']);
            if exist(inputPath, 'file') ~= 2
                error('GeometryLoader:MissingInput', ...
                    'Input script not found: "%s".', inputPath);
            end

            % Run in local function workspace so that variables defined in the
            % script are visible here but do not leak elsewhere.
            try
                run(inputPath);
            catch ME
                error('GeometryLoader:InputScriptError', ...
                    'Failed to run input script "%s" for hull "%s".\nOriginal error: %s', ...
                    inputScript, baseName, ME.message);
            end

            % Accept either:
            %   - inputData struct
            %   - INPUT struct
            if exist('inputData', 'var') == 1 && isstruct(inputData)
                candidate = inputData;
            elseif exist('INPUT', 'var') == 1 && isstruct(INPUT)
                candidate = INPUT;
            else
                error('GeometryLoader:MissingInputStruct', ...
                    ['Input file for "%s" must define either a struct named INPUT ' ...
                     'or a struct named inputData.'], baseName);
            end

            % Validate required fields
            requiredFields = {'VCG'};
            missing = requiredFields(~isfield(candidate, requiredFields));

            if ~isempty(missing)
                error('GeometryLoader:MissingRequiredInputField', ...
                    'Input file for "%s" is missing required field(s): %s.', ...
                    baseName, strjoin(missing, ', '));
            end

            % Basic type validation
            if ~isnumeric(candidate.VCG) || ~isscalar(candidate.VCG) || ~isfinite(candidate.VCG)
                error('GeometryLoader:InvalidVCG', ...
                    'INPUT.VCG or inputData.VCG for "%s" must be a finite numeric scalar.', baseName);
            end

            inputData = candidate;
        end

        function inputScript = resolveInputScript(baseName, hullsDir)
            stripped1 = regexprep(baseName, '_hull$', '');
            stripped2 = regexprep(baseName, '_mesh$', '');
            stripped3 = regexprep(stripped2, '_hull$', '');

            candidates = { ...
                ['Input_', baseName], ...
                ['Input_', stripped1], ...
                ['Input_', stripped2], ...
                ['Input_', stripped3] ...
            };

            candidates = unique(candidates, 'stable');

            for k = 1:numel(candidates)
                candidatePath = fullfile(hullsDir, [candidates{k}, '.m']);
                if exist(candidatePath, 'file') == 2
                    inputScript = candidates{k};
                    return;
                end
            end

            error('GeometryLoader:MissingInput', ...
                'No input script found for hull "%s" in folder "%s".', baseName, hullsDir);
        end

        function [meshFile, meshFormat] = resolveMeshFile(baseName, hullsDir)
            candidates = { ...
                struct('path', fullfile(hullsDir, [baseName, '_mesh.obj']), 'fmt', 'obj'), ...
                struct('path', fullfile(hullsDir, [baseName, '_mesh.stl']), 'fmt', 'stl'), ...
                struct('path', fullfile(hullsDir, [baseName, '_mesh.mat']), 'fmt', 'mat'), ...
                struct('path', fullfile(hullsDir, [baseName, '.obj']),      'fmt', 'obj'), ...
                struct('path', fullfile(hullsDir, [baseName, '.stl']),      'fmt', 'stl'), ...
                struct('path', fullfile(hullsDir, [baseName, '.mat']),      'fmt', 'mat') ...
            };

            for k = 1:numel(candidates)
                if exist(candidates{k}.path, 'file') == 2
                    meshFile   = candidates{k}.path;
                    meshFormat = candidates{k}.fmt;
                    return;
                end
            end

            expectedFiles = cellfun(@(c) c.path, candidates, 'UniformOutput', false);

            error('GeometryLoader:MissingGeometry', ...
                'Mesh file not found for "%s". Expected one of:%s\n  - %s', ...
                baseName, newline, strjoin(expectedFiles, [newline, '  - ']));
        end

        function mesh = normalizeMesh(mesh, cleanMesh, tol)
            if ~isstruct(mesh)
                error('GeometryLoader:InvalidMesh', ...
                    'Loaded mesh must be a struct.');
            end
            if ~isfield(mesh, 'vertices') || ~isfield(mesh, 'faces')
                error('GeometryLoader:InvalidMesh', ...
                    'Loaded mesh must contain fields "vertices" and "faces".');
            end

            V = double(mesh.vertices);
            F = double(mesh.faces);

            if isempty(V) || isempty(F)
                error('GeometryLoader:EmptyMesh', ...
                    'Loaded mesh contains empty vertices or faces.');
            end

            if ~ismatrix(V) || size(V, 2) ~= 3
                error('GeometryLoader:InvalidVertexArray', ...
                    'Mesh vertices must be an Nx3 numeric array.');
            end
            if ~ismatrix(F) || size(F, 2) ~= 3
                error('GeometryLoader:InvalidFaceArray', ...
                    'Mesh faces must be an Mx3 numeric array of triangle indices.');
            end

            if any(~isfinite(V), 'all')
                error('GeometryLoader:NonFiniteVertices', ...
                    'Mesh vertices contain NaN or Inf values.');
            end
            if any(~isfinite(F), 'all')
                error('GeometryLoader:NonFiniteFaces', ...
                    'Mesh faces contain NaN or Inf values.');
            end

            % Face indices must be integer-valued
            if any(abs(F - round(F)) > 0, 'all')
                error('GeometryLoader:NonIntegerFaceIndices', ...
                    'Mesh faces contain non-integer indices.');
            end
            F = round(F);

            nV = size(V, 1);
            if any(F(:) < 1) || any(F(:) > nV)
                badMin = min(F(:));
                badMax = max(F(:));
                error('GeometryLoader:FaceIndexOutOfRange', ...
                    ['Mesh faces reference vertices outside valid range [1, %d]. ' ...
                     'Observed min/max indices: [%d, %d].'], nV, badMin, badMax);
            end

            if cleanMesh
                % Remove repeated vertex references within the same face
                same12 = F(:,1) == F(:,2);
                same13 = F(:,1) == F(:,3);
                same23 = F(:,2) == F(:,3);
                nonDegenerateByIndex = ~(same12 | same13 | same23);

                if any(~nonDegenerateByIndex)
                    warning('GeometryLoader:DegenerateFacesRemoved', ...
                        'Removed %d faces with repeated vertex indices.', ...
                        nnz(~nonDegenerateByIndex));
                    F = F(nonDegenerateByIndex, :);
                end

                % Remove geometrically degenerate triangles
                A = V(F(:,1), :);
                B = V(F(:,2), :);
                C = V(F(:,3), :);
                triArea2 = vecnorm(cross(B - A, C - A, 2), 2, 2); % twice area
                nonDegenerateByArea = triArea2 > tol;

                if any(~nonDegenerateByArea)
                    warning('GeometryLoader:ZeroAreaFacesRemoved', ...
                        'Removed %d near-zero-area faces.', nnz(~nonDegenerateByArea));
                    F = F(nonDegenerateByArea, :);
                end

                % Remove exact duplicate faces irrespective of node order
                Fs = sort(F, 2);
                [~, ia] = unique(Fs, 'rows', 'stable');
                if numel(ia) < size(F, 1)
                    warning('GeometryLoader:DuplicateFacesRemoved', ...
                        'Removed %d duplicate faces.', size(F,1) - numel(ia));
                    F = F(ia, :);
                end
            end

            if isempty(F)
                error('GeometryLoader:NoValidFaces', ...
                    'Mesh does not contain any valid triangular faces after normalization.');
            end

            mesh.vertices  = V;
            mesh.faces     = F;
            mesh.nVertices = size(V, 1);
            mesh.nFaces    = size(F, 1);
        end

        function mesh = readObjMesh(objPath, allowZeroBasedObj)
            fid = fopen(objPath, 'r');
            if fid < 0
                error('GeometryLoader:ObjOpenFailed', ...
                    'Cannot open OBJ file: "%s".', objPath);
            end

            cleaner = onCleanup(@() fclose(fid));
            C = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
            lines = C{1};

            vertices = zeros(0, 3);
            faces    = zeros(0, 3);

            nNormals   = 0;
            nTexCoords = 0;
            nObjects   = 0;
            nGroups    = 0;

            for i = 1:numel(lines)
                line = strtrim(lines{i});
                if isempty(line) || startsWith(line, '#')
                    continue;
                end

                if startsWith(line, 'v ')
                    vals = sscanf(line(2:end), '%f %f %f');
                    if numel(vals) < 3
                        warning('GeometryLoader:MalformedVertexLine', ...
                            'Ignoring malformed vertex line %d in "%s".', i, objPath);
                        continue;
                    end
                    vertices(end+1, :) = vals(1:3).'; %#ok<AGROW>

                elseif startsWith(line, 'vn ')
                    nNormals = nNormals + 1;

                elseif startsWith(line, 'vt ')
                    nTexCoords = nTexCoords + 1;

                elseif startsWith(line, 'o ')
                    nObjects = nObjects + 1;

                elseif startsWith(line, 'g ')
                    nGroups = nGroups + 1;

                elseif startsWith(line, 'f ')
                    tokens = strsplit(strtrim(line(2:end)));
                    if numel(tokens) < 3
                        warning('GeometryLoader:MalformedFaceLine', ...
                            'Ignoring malformed face line %d in "%s".', i, objPath);
                        continue;
                    end

                    idx = zeros(numel(tokens), 1);
                    for p = 1:numel(tokens)
                        idx(p) = GeometryLoader.parseObjVertexIndex(tokens{p}, size(vertices,1), allowZeroBasedObj, objPath, i);
                    end

                    if any(isnan(idx))
                        warning('GeometryLoader:InvalidFaceSkipped', ...
                            'Skipping invalid face at line %d in "%s".', i, objPath);
                        continue;
                    end

                    % Fan triangulation
                    for p = 2:(numel(idx)-1)
                        tri = [idx(1), idx(p), idx(p+1)];
                        if numel(unique(tri)) == 3
                            faces(end+1, :) = tri; %#ok<AGROW>
                        end
                    end
                end
            end

            if isempty(vertices)
                error('GeometryLoader:NoVerticesInOBJ', ...
                    'OBJ file "%s" contains no vertices.', objPath);
            end
            if isempty(faces)
                error('GeometryLoader:NoFacesInOBJ', ...
                    'OBJ file "%s" contains no valid faces.', objPath);
            end

            if nNormals > 0 || nTexCoords > 0 || nObjects > 0 || nGroups > 0
                fprintf(['  [OBJ] V=%d, F=%d, VN=%d, VT=%d, Objects=%d, Groups=%d\n'], ...
                    size(vertices, 1), size(faces, 1), nNormals, nTexCoords, nObjects, nGroups);
            end

            mesh = struct('vertices', vertices, 'faces', faces);
        end

        function idx = parseObjVertexIndex(token, nVerticesSoFar, allowZeroBasedObj, objPath, lineNumber)
            % Parse OBJ face token forms:
            %   f v1 v2 v3
            %   f v1/vt1 v2/vt2 v3/vt3
            %   f v1//vn1 ...
            %   f v1/vt1/vn1 ...
            %
            % Only the vertex index is used here.

            parts = strsplit(token, '/');
            head  = strtrim(parts{1});

            if isempty(head)
                idx = NaN;
                return;
            end

            raw = str2double(head);
            if ~isfinite(raw) || abs(raw - round(raw)) > 0
                idx = NaN;
                return;
            end
            raw = round(raw);

            if raw > 0
                % Standard OBJ 1-based index
                idx = raw;

            elseif raw < 0
                % Relative indexing: -1 means most recent vertex
                idx = nVerticesSoFar + raw + 1;

            else
                % raw == 0 is invalid in standard OBJ
                if allowZeroBasedObj
                    idx = 1; % placeholder, corrected below
                    % In a 0-based exporter, token 0 refers to the first vertex
                    idx = raw + 1;
                else
                    error('GeometryLoader:ZeroBasedOBJIndexNotAllowed', ...
                        ['OBJ file "%s" line %d contains vertex index 0. ' ...
                         'Standard OBJ is 1-based. Set cfg.allowZeroBasedObj = true ' ...
                         'only if you know this file uses nonstandard 0-based indexing.'], ...
                         objPath, lineNumber);
                end
            end

            if idx < 1 || idx > nVerticesSoFar
                error('GeometryLoader:OBJIndexOutOfRange', ...
                    ['OBJ file "%s" line %d references vertex %d, but only %d vertices ' ...
                     'have been defined at that point.'], ...
                    objPath, lineNumber, idx, nVerticesSoFar);
            end
        end

        function mesh = readStlMesh(stlPath)
            stlreadAvailable = exist('stlread', 'file') == 2 || ...
                               exist('stlread', 'builtin') == 5;

            if stlreadAvailable
                try
                    out = stlread(stlPath);

                    if isa(out, 'triangulation')
                        mesh = struct( ...
                            'vertices', double(out.Points), ...
                            'faces',    double(out.ConnectivityList));
                        return;
                    else
                        [F, V] = stlread(stlPath);
                        mesh = struct( ...
                            'vertices', double(V), ...
                            'faces',    double(F));
                        return;
                    end
                catch ME
                    warning('GeometryLoader:StlReadFailed', ...
                        'stlread failed for "%s": %s. Falling back to internal parser.', ...
                        stlPath, ME.message);
                end
            end

            [isBinaryLikely, fileSize, nTriHeader] = GeometryLoader.detectBinaryStl(stlPath);

            if isBinaryLikely
                error('GeometryLoader:BinarySTLRequiresReader', ...
                    ['File "%s" appears to be a binary STL (size=%d bytes, header count=%d). ' ...
                     'Install a working stlread or convert the file to ASCII STL.'], ...
                    stlPath, fileSize, nTriHeader);
            end

            mesh = GeometryLoader.readAsciiStl(stlPath);
        end

        function [isBinaryLikely, fileSize, nTriHeader] = detectBinaryStl(stlPath)
            fid = fopen(stlPath, 'r');
            if fid < 0
                error('GeometryLoader:StlOpenFailed', ...
                    'Cannot open STL file: "%s".', stlPath);
            end
            cleaner = onCleanup(@() fclose(fid));

            header = fread(fid, 80, '*uint8');
            fseek(fid, 0, 'eof');
            fileSize = ftell(fid);

            isBinaryLikely = false;
            nTriHeader = NaN;

            if fileSize < 84
                return;
            end

            fseek(fid, 0, 'bof');
            firstBytes = fread(fid, min(fileSize, 256), '*char')';
            startsWithSolid = startsWith(strtrim(firstBytes), 'solid');

            fseek(fid, 80, 'bof');
            nTriHeader = fread(fid, 1, 'uint32');

            expectedSize = 84 + 50 * double(nTriHeader);

            % If the file size exactly matches binary STL layout, treat it as binary.
            if expectedSize == fileSize
                isBinaryLikely = true;
                return;
            end

            % If it does not start with "solid", it is very likely binary.
            if ~startsWithSolid
                isBinaryLikely = true;
            end
        end

        function mesh = readAsciiStl(stlPath)
            fid = fopen(stlPath, 'r');
            if fid < 0
                error('GeometryLoader:StlOpenFailed', ...
                    'Cannot open STL file: "%s".', stlPath);
            end
            cleaner = onCleanup(@() fclose(fid));

            C = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
            lines = C{1};

            verts = zeros(0, 3);

            for i = 1:numel(lines)
                line = strtrim(lines{i});
                if startsWith(line, 'vertex ')
                    vals = sscanf(line(7:end), '%f %f %f');
                    if numel(vals) ~= 3
                        warning('GeometryLoader:MalformedStlVertex', ...
                            'Ignoring malformed STL vertex at line %d in "%s".', i, stlPath);
                        continue;
                    end
                    verts(end+1, :) = vals.'; %#ok<AGROW>
                end
            end

            if isempty(verts)
                error('GeometryLoader:AsciiStlParseFailed', ...
                    'ASCII STL parse failed for "%s": no vertex records were found.', stlPath);
            end
            if mod(size(verts,1), 3) ~= 0
                error('GeometryLoader:AsciiStlParseFailed', ...
                    ['ASCII STL parse failed for "%s": number of vertex records (%d) ' ...
                     'is not a multiple of 3.'], stlPath, size(verts,1));
            end

            nTri = size(verts,1) / 3;

            % Deduplicate vertices
            [V, ~, ic] = unique(verts, 'rows', 'stable');
            F = reshape(ic, 3, nTri).';

            if isempty(V) || isempty(F)
                error('GeometryLoader:AsciiStlParseFailed', ...
                    'ASCII STL parse failed for "%s": could not build a valid mesh.', stlPath);
            end

            mesh = struct('vertices', double(V), 'faces', double(F));
        end
    end
end