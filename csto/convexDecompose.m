function new_faces2 = convexDecompose(new_faces1, intersection_points, obstacle_params)

    new_faces2 = new_faces1;

    for i = 1:numel(new_faces1)
        voronoi_cell = new_faces1{i};

        % No (or insufficient) intersection points -> no decomposition
        if isempty(intersection_points{i}) || size(intersection_points{i},1) < 3
            continue;
        end

        nI = size(intersection_points{i}, 1);

        if nI <= 6
            % ===== Branch A: keep your original ASV idea (wrapped in a function) =====
            piecesCells = asvDecomposeOneCell(voronoi_cell, intersection_points{i}, obstacle_params);

            % If no multiple pieces are produced, keep the original behavior
            if isempty(piecesCells)
                continue;
            end

            % Append new pieces and remove the original cell
            for k = 1:numel(piecesCells)
                new_faces2{end+1} = piecesCells{k}; %#ok<AGROW>
            end
            new_faces2{i} = [];

        else
            % ===== Branch B: too many intersections -> fast box-difference decomposition =====
            [cellMin, cellMax] = cellToAABB(voronoi_cell);

            % Find intersecting obstacles
            hitObs = [];
            for k = 1:numel(obstacle_params)
                obsMin = obstacle_params(k).start;
                obsMax = obsMin + obstacle_params(k).size;
                if aabbIntersect(cellMin, cellMax, obsMin, obsMax)
                    hitObs(end+1) = k; %#ok<AGROW>
                end
            end

            if isempty(hitObs)
                continue;
            end

            % Sequentially subtract obstacles
            pieces = {struct('min', cellMin, 'max', cellMax)};

            for kk = 1:numel(hitObs)
                k = hitObs(kk);
                obsMin = obstacle_params(k).start;
                obsMax = obsMin + obstacle_params(k).size;

                newPieces = {};
                for p = 1:numel(pieces)
                    P = pieces{p};
                    subPieces = subtractBoxByBox(P.min, P.max, obsMin, obsMax);
                    newPieces = [newPieces, subPieces]; %#ok<AGROW>
                end
                pieces = newPieces;

                if isempty(pieces)
                    break;
                end
            end

            if isempty(pieces)
                new_faces2{i} = [];
                continue;
            end

            % Convert each box piece to the "cell" format
            for p = 1:numel(pieces)
                P = pieces{p};
                if any((P.max - P.min) < 1e-10)
                    continue;
                end
                new_faces2{end+1} = boxToCell(P.min, P.max); %#ok<AGROW>
            end

            new_faces2{i} = [];
        end
    end

    new_faces2 = new_faces2(~cellfun('isempty', new_faces2));
end


function piecesCells = asvDecomposeOneCell(voronoi_cell, intersection_points_i, obstacle_params)

    piecesCells = {};

    % Merge points (face vertices + intersection points)
    faces_matrix = vertcat(voronoi_cell{2:end});
    combinedPoints = [faces_matrix; intersection_points_i];

    tol = 1e-6;
    combinedPoints = unique(round(combinedPoints/tol)*tol, 'rows');

    % Add obstacle surface/inside points for stabilization (your original logic)
    obstacle_points1 = classify_points_in_obstacle(combinedPoints, obstacle_params);

    faces_cells = {};
    [~, faces_cells] = recursive_asv_full(combinedPoints, 0, faces_cells, obstacle_points1);

    % Need at least positive/negative layers, otherwise concave handling is not possible
    if numel(faces_cells) < 2
        return;
    end

    % Handle concave faces and generate convex hull face sets
    [merged_faces, ~] = processConcaveFaces(faces_cells);
    if isempty(merged_faces)
        return;
    end

    % Coplanar merge (key fix: do NOT force 1:4 truncation)
    merged_faces = mergeCoplanarFacesInBodies(merged_faces);

    % merged_faces may contain multiple bodies, each is {face1, face2, ...}
    if iscell(merged_faces) && numel(merged_faces) > 1
        for j = 1:numel(merged_faces)
            current_faces = merged_faces{j};
            all_points = unique(vertcat(current_faces{:}), 'rows');
            center = mean(all_points, 1);
            piecesCells{end+1} = [{center}; current_faces(:)]; %#ok<AGROW>
        end
    else
        % If only one body is produced, keep original behavior: do not replace
        return;
    end
end

function merged_faces = mergeCoplanarFacesInBodies(merged_faces)

    for l = 1:numel(merged_faces)
        polyhedron = merged_faces{l};
        combined_faces = {};
        processed_faces = false(length(polyhedron), 1);

        for m = 1:length(polyhedron)
            if processed_faces(m)
                continue;
            end

            face1 = polyhedron{m};
            normal1 = computeNormal(face1);
            coplanar_vertices = face1;

            for n = m+1:length(polyhedron)
                if processed_faces(n)
                    continue;
                end
                face2 = polyhedron{n};
                normal2 = computeNormal(face2);

                if isCoplanar(normal1, normal2, face1, face2)
                    coplanar_vertices = [coplanar_vertices; face2]; %#ok<AGROW>
                    processed_faces(n) = true;
                end
            end

            coplanar_vertices = unique(coplanar_vertices, 'rows');

            % Output by 2D convex-hull ordering (may be 5/6-gon; do not truncate)
            sorted_vertices = sortCoplanarVertices3D(coplanar_vertices);

            % convhull returns a closed polygon; drop repeated endpoint
            if size(sorted_vertices,1) >= 2 && all(abs(sorted_vertices(1,:) - sorted_vertices(end,:)) < 1e-12)
                sorted_vertices(end,:) = [];
            end

            combined_faces{end+1} = sorted_vertices; %#ok<AGROW>
        end

        merged_faces{l} = combined_faces(:);
    end
end

%% Your recursive ASV (fix: pass obstacle_points1 in recursive call)
function [node, faces_cells] = recursive_asv_full(combinedPoints, depth, faces_cells, obstacle_points1)

    if isempty(combinedPoints)
        node = [];
        return;
    end

    combinedPoints = unique(combinedPoints, 'rows');

    node = struct('Vertices', combinedPoints, 'Depth', depth, ...
                  'IsPositive', mod(depth, 2) == 0, ...
                  'Children', [], 'Volume', []);

    try
        [K, volume] = convhulln(combinedPoints);
    catch
        node = [];
        return;
    end

    hullVertices = unique(K(:));
    convexHullPoints = combinedPoints(hullVertices, :);

    remainingPoints = setdiff(combinedPoints, convexHullPoints, 'rows');

    % If too few points remain, add obstacle points (your original logic)
    if size(remainingPoints, 1) < 7 && size(remainingPoints, 1) > 0
        remainingPoints = unique([remainingPoints; obstacle_points1], 'rows');
    end

    node.Volume = volume;

    % Save faces for current layer
    faces_cells{end+1} = calculate_faces_vertices_as_cells(combinedPoints, K);

    if isempty(remainingPoints)
        node.Children = [];
    else
        [child_node, faces_cells] = recursive_asv_full(remainingPoints, depth + 1, faces_cells, obstacle_points1);
        node.Children = {child_node};
    end
end

function faces_cells = calculate_faces_vertices_as_cells(vertices, K)
% Merge coplanar triangles (keep your original idea)
    numFaces = size(K, 1);
    faces_cells = {};
    faceNormals = zeros(numFaces, 3);
    facePoints = vertices(K(:, 1), :);

    for i = 1:numFaces
        faceVertices = vertices(K(i, :), :);
        v1 = faceVertices(2, :) - faceVertices(1, :);
        v2 = faceVertices(3, :) - faceVertices(1, :);
        n = cross(v1, v2);
        if norm(n) < 1e-12
            n = [0 0 0];
        else
            n = n / norm(n);
        end
        faceNormals(i, :) = n;
    end

    merged = false(numFaces, 1);

    for i = 1:numFaces
        if merged(i), continue; end

        currentFaceVertices = vertices(K(i, :), :);

        for j = i+1:numFaces
            if merged(j), continue; end
            if abs(dot(faceNormals(i, :), faceNormals(j, :))) > 0.999
                d = dot(faceNormals(i, :), (facePoints(j, :) - facePoints(i, :)));
                if abs(d) < 1e-5
                    currentFaceVertices = unique([currentFaceVertices; vertices(K(j, :), :)], 'rows', 'stable');
                    merged(j) = true;
                end
            end
        end

        faces_cells{end+1} = num2cell(currentFaceVertices, 2);
    end
end

%% Concave-face processing (kept as-is; plotting remains commented out)
function [merged_faces, removed_faces] = processConcaveFaces(faces_cells)

    concaveFaces = faces_cells{2};

    convexVertices = [];
    for i = 1:length(faces_cells{1})
        convexVertices = [convexVertices; cell2mat(faces_cells{1}{i})]; %#ok<AGROW>
    end
    convexVertices = unique(convexVertices, 'rows');

    concaveVertices = [];
    for i = 1:length(faces_cells{2})
        concaveVertices = [concaveVertices; cell2mat(faces_cells{2}{i})]; %#ok<AGROW>
    end
    concaveVertices = unique(concaveVertices, 'rows');

    merged_faces = {};
    removed_faces = {};

    for i = 1:length(concaveFaces)
        faceVertices = cell2mat(concaveFaces{i});

        if size(faceVertices,1) < 3
            removed_faces{end+1} = faceVertices; %#ok<AGROW>
            continue;
        end

        normalVector = cross(faceVertices(2,:) - faceVertices(1,:), faceVertices(3,:) - faceVertices(1,:));
        if norm(normalVector) < 1e-12
            removed_faces{end+1} = faceVertices; %#ok<AGROW>
            continue;
        end
        normalVector = normalVector / norm(normalVector);

        faceCenter = mean(faceVertices, 1);

        [posConvexSide, negConvexSide, posConcaveSide, negConcaveSide] = ...
            classifyVerticesBySide(faceCenter, convexVertices, concaveVertices, normalVector);

        if isempty(posConcaveSide) && ~isempty(posConvexSide) && ~isempty(negConcaveSide)
            combinedVertices = [faceVertices; posConvexSide];
        elseif isempty(negConcaveSide) && ~isempty(negConvexSide) && ~isempty(posConcaveSide)
            combinedVertices = [faceVertices; negConvexSide];
        else
            removed_faces{end+1} = faceVertices; %#ok<AGROW>
            continue;
        end

        combinedVertices = unique(combinedVertices, 'rows');

        try
            K = convhulln(combinedVertices);
        catch
            removed_faces{end+1} = faceVertices; %#ok<AGROW>
            continue;
        end

        current_convex_hull = {};
        for j = 1:size(K, 1)
            face = combinedVertices(K(j, :), :);
            current_convex_hull{end+1} = face; %#ok<AGROW>
            % plotFace(face); % Enable if visualization is needed
        end

        merged_faces{end+1} = current_convex_hull; %#ok<AGROW>
    end
end

function [posConvexSide, negConvexSide, posConcaveSide, negConcaveSide] = ...
    classifyVerticesBySide(faceCenter, convexVertices, concaveVertices, normalVector)

    posConvexSide = []; negConvexSide = [];
    posConcaveSide = []; negConcaveSide = [];

    for j = 1:size(convexVertices, 1)
        vertex = convexVertices(j, :);
        projection = dot(vertex - faceCenter, normalVector);
        if projection > 0
            posConvexSide = [posConvexSide; vertex]; %#ok<AGROW>
        elseif projection < 0
            negConvexSide = [negConvexSide; vertex]; %#ok<AGROW>
        end
    end

    for j = 1:size(concaveVertices, 1)
        vertex = concaveVertices(j, :);
        projection = dot(vertex - faceCenter, normalVector);
        if projection > 0
            posConcaveSide = [posConcaveSide; vertex]; %#ok<AGROW>
        elseif projection < 0
            negConcaveSide = [negConcaveSide; vertex]; %#ok<AGROW>
        end
    end
end

%% Obstacle-point classification (kept as-is; comments translated)
function obstacle_points1 = classify_points_in_obstacle(points, obstacle_params)

    obstacle_points1 = [];
    tolerance = 1e-5;

    for p_idx = 1:size(points,1)
        point = points(p_idx, :);

        is_inside_or_on_surface = false;

        for i = 1:length(obstacle_params)
            start = obstacle_params(i).start;
            sz    = obstacle_params(i).size;

            min_bounds = start;
            max_bounds = start + sz;

            if (point(1) >= min_bounds(1) - tolerance && point(1) <= max_bounds(1) + tolerance) && ...
               (point(2) >= min_bounds(2) - tolerance && point(2) <= max_bounds(2) + tolerance) && ...
               (point(3) >= min_bounds(3) - tolerance && point(3) <= max_bounds(3) + tolerance)

                % On surface
                if (abs(point(1) - min_bounds(1)) <= tolerance || abs(point(1) - max_bounds(1)) <= tolerance || ...
                    abs(point(2) - min_bounds(2)) <= tolerance || abs(point(2) - max_bounds(2)) <= tolerance || ...
                    abs(point(3) - min_bounds(3)) <= tolerance || abs(point(3) - max_bounds(3)) <= tolerance)

                    is_inside_or_on_surface = true;

                % Strictly inside
                elseif (point(1) > min_bounds(1) + tolerance && point(1) < max_bounds(1) - tolerance && ...
                        point(2) > min_bounds(2) + tolerance && point(2) < max_bounds(2) - tolerance && ...
                        point(3) > min_bounds(3) + tolerance && point(3) < max_bounds(3) - tolerance)

                    is_inside_or_on_surface = true;
                end
            end

            if is_inside_or_on_surface
                break;
            end
        end

        if is_inside_or_on_surface
            obstacle_points1 = [obstacle_points1; point]; %#ok<AGROW>
        end
    end

    obstacle_points1 = unique(obstacle_points1, 'rows');
end

%% Coplanarity / normals / sorting (logic unchanged; translate messages)
function normal = computeNormal(face)
    if size(face, 1) < 3
        error('A face must contain at least three points to compute the normal.');
    end
    normal = cross(face(2,:) - face(1,:), face(3,:) - face(1,:));
    normal = normal / max(norm(normal), 1e-12);
end

function is_coplanar = isCoplanar(normal1, normal2, face1, face2)
    if abs(dot(normal1, normal2)) < 0.999
        is_coplanar = false;
        return;
    end

    d = -dot(normal1, face1(1, :));
    distances = abs(dot(face2, repmat(normal1, size(face2, 1), 1), 2) + d);
    is_coplanar = all(distances < 1e-10);
end

function sorted_vertices = sortCoplanarVertices3D(coplanar_vertices)

    if size(coplanar_vertices, 1) < 3
        error('At least three points are required to determine a consistent ordering.');
    end

    normal = computeNormal(coplanar_vertices);

    [~, maxIdx] = max(abs(normal));
    switch maxIdx
        case 1
            projPoints = coplanar_vertices(:, [2, 3]);
        case 2
            projPoints = coplanar_vertices(:, [1, 3]);
        otherwise
            projPoints = coplanar_vertices(:, [1, 2]);
    end

    k = convhull(projPoints(:, 1), projPoints(:, 2));
    sorted_vertices = coplanar_vertices(k, :);
end

%% =======================================================================
%% Branch B: AABB box-difference decomposition (stable for multi-obstacle)
function [minP, maxP] = cellToAABB(cellFaces)
    V = [];
    for t = 2:numel(cellFaces)
        V = [V; cellFaces{t}]; %#ok<AGROW>
    end
    minP = min(V, [], 1);
    maxP = max(V, [], 1);
end

function flag = aabbIntersect(minA, maxA, minB, maxB)
    flag = ~( ...
        maxA(1) <= minB(1) || minA(1) >= maxB(1) || ...
        maxA(2) <= minB(2) || minA(2) >= maxB(2) || ...
        maxA(3) <= minB(3) || minA(3) >= maxB(3) );
end

function pieces = subtractBoxByBox(pMin, pMax, oMin, oMax)
% Return convex box pieces for pBox \ (pBox ∩ oBox).
% Each piece is a struct('min', ..., 'max', ...).

    pieces = {};

    iMin = max(pMin, oMin);
    iMax = min(pMax, oMax);

    % No intersection -> keep the original box
    if any(iMin >= iMax)
        pieces = {struct('min', pMin, 'max', pMax)};
        return;
    end

    % Fully covered -> removed
    if all(iMin <= pMin + 1e-15) && all(iMax >= pMax - 1e-15)
        pieces = {};
        return;
    end

    % Left
    if pMin(1) < iMin(1)
        pieces{end+1} = struct('min', [pMin(1), pMin(2), pMin(3)], ...
                              'max', [iMin(1), pMax(2), pMax(3)]);
    end
    % Right
    if iMax(1) < pMax(1)
        pieces{end+1} = struct('min', [iMax(1), pMin(2), pMin(3)], ...
                              'max', [pMax(1), pMax(2), pMax(3)]);
    end

    xMin = max(pMin(1), iMin(1));
    xMax = min(pMax(1), iMax(1));

    % Front
    if pMin(2) < iMin(2)
        pieces{end+1} = struct('min', [xMin, pMin(2), pMin(3)], ...
                              'max', [xMax, iMin(2), pMax(3)]);
    end
    % Back
    if iMax(2) < pMax(2)
        pieces{end+1} = struct('min', [xMin, iMax(2), pMin(3)], ...
                              'max', [xMax, pMax(2), pMax(3)]);
    end

    yMin = max(pMin(2), iMin(2));
    yMax = min(pMax(2), iMax(2));

    % Bottom
    if pMin(3) < iMin(3)
        pieces{end+1} = struct('min', [xMin, yMin, pMin(3)], ...
                              'max', [xMax, yMax, iMin(3)]);
    end
    % Top
    if iMax(3) < pMax(3)
        pieces{end+1} = struct('min', [xMin, yMin, iMax(3)], ...
                              'max', [xMax, yMax, pMax(3)]);
    end

    % Filter degenerate boxes
    keep = true(1, numel(pieces));
    for k = 1:numel(pieces)
        if any(pieces{k}.max - pieces{k}.min <= 1e-12)
            keep(k) = false;
        end
    end
    pieces = pieces(keep);
end

function cellStruct = boxToCell(minP, maxP)
% Convert an AABB box to your {center; 6 faces} cell format.

    x0=minP(1); y0=minP(2); z0=minP(3);
    x1=maxP(1); y1=maxP(2); z1=maxP(3);

    v000 = [x0 y0 z0]; v100 = [x1 y0 z0];
    v110 = [x1 y1 z0]; v010 = [x0 y1 z0];
    v001 = [x0 y0 z1]; v101 = [x1 y0 z1];
    v111 = [x1 y1 z1]; v011 = [x0 y1 z1];

    f1 = [v000; v100; v110; v010];
    f2 = [v001; v101; v111; v011];
    f3 = [v000; v100; v101; v001];
    f4 = [v010; v110; v111; v011];
    f5 = [v000; v010; v011; v001];
    f6 = [v100; v110; v111; v101];

    center = (minP + maxP)/2;
    cellStruct = [{center}; {f1}; {f2}; {f3}; {f4}; {f5}; {f6}];
end


function total_volume = calculate_volume(node)

    % Leaf node: return signed volume directly
    if isempty(node.Children)
        total_volume = node.Volume * (2 * node.IsPositive - 1);  % positive part -> +, negative part -> -
        return;
    end

    % Non-leaf: accumulate children volumes
    total_volume = node.Volume * (2 * node.IsPositive - 1);
    for i = 1:length(node.Children)
        total_volume = total_volume + calculate_volume(node.Children{i});
    end
end


