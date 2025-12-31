function [path_points, path_faces] = chainpath( ...
    new_faces2, start_point, end_point, obstacle_params, cfg)

    if nargin < 5 || isempty(cfg)
        cfg.useDetourCost = false;
    end
    if nargin < 4
        obstacle_params = [];
    end

    % ===== cost weights =====
    alpha = 10;     % heuristic weight
    beta  = 1;    % detour (topological) cost weight

    % ===== find start / goal blocks =====
    start_faces_idx = find_closest_face(new_faces2, start_point);
    end_face_idx    = find_closest_face(new_faces2, end_point);

    start_face_idx = start_faces_idx(1);
    excluded_faces = start_faces_idx(2:end);

    hfun = @(p, q) norm(p - q);

    % ===== build projection basis for detour computation =====
    [u, v] = buildOrthonormalBasis(end_point - start_point);

    % ===== open list initialization =====
    openList = start_face_idx;
    openList_path = {start_face_idx};

    start_center = get_center(new_faces2{start_face_idx});
    goal_center  = get_center(new_faces2{end_face_idx(1)});

    g0 = 0;
    h0 = hfun(start_center, goal_center);

    if cfg.useDetourCost
        k0 = zeros(1, numel(obstacle_params));  % winding numbers
        topo0 = sum(abs(k0));
    else
        k0 = [];
        topo0 = 0;
    end

    f0 = g0 + alpha * h0 + beta * topo0;

    openList_cost = [g0, h0, topo0, f0];
    openList_k    = {k0};

    visited = false(1, length(new_faces2));
    visited(start_face_idx) = true;
    for i = 1:numel(excluded_faces)
        visited(excluded_faces(i)) = true;
    end

    % ===== main search loop =====
    while ~isempty(openList)

        [~, min_idx] = min(openList_cost(:,4));
        current_idx = openList(min_idx);

        current_center = get_center(new_faces2{current_idx});
        current_g = openList_cost(min_idx,1);
        current_k = openList_k{min_idx};

        % ----- goal reached -----
        if any(current_idx == end_face_idx)
            path_faces = openList_path{min_idx};
            path_points = cellfun(@(ii) get_center(new_faces2{ii}), ...
                                   num2cell(path_faces), ...
                                   'UniformOutput', false);
            path_points = vertcat(path_points{:});
            return;
        end

        visited(current_idx) = true;
        neighbors = get_neighbors(current_idx, new_faces2, visited);

        for t = 1:length(neighbors)
            neighbor_idx = neighbors(t);
            neighbor_center = get_center(new_faces2{neighbor_idx});

            % ----- geometric cost -----
            stepCost = norm(current_center - neighbor_center);
            g = current_g + stepCost;

            % ----- heuristic -----
            h_val = hfun(neighbor_center, goal_center);

            % ----- detour / winding update -----
            if cfg.useDetourCost && ~isempty(obstacle_params)
                dk = detourIncrement3D(current_center, neighbor_center, ...
                                       obstacle_params, u, v);
                k_new = current_k + dk;
                topo = sum(abs(k_new));
            else
                k_new = current_k;
                topo = 0;
            end

            f = g + alpha * h_val + beta * topo;

            % ----- open list update -----
            [in_open, idx_open] = ismember(neighbor_idx, openList);
            if in_open
                if f < openList_cost(idx_open,4)
                    openList_cost(idx_open,:) = [g, h_val, topo, f];
                    openList_path{idx_open} = ...
                        [openList_path{min_idx}, neighbor_idx];
                    openList_k{idx_open} = k_new;
                end
            else
                openList(end+1) = neighbor_idx; %#ok<AGROW>
                openList_cost(end+1,:) = [g, h_val, topo, f]; %#ok<AGROW>
                openList_path{end+1} = ...
                    [openList_path{min_idx}, neighbor_idx]; %#ok<AGROW>
                openList_k{end+1} = k_new; %#ok<AGROW>
            end
        end

        % ----- remove current node from open list -----
        openList(min_idx) = [];
        openList_cost(min_idx,:) = [];
        openList_path(min_idx) = [];
        openList_k(min_idx) = [];
    end

    % No feasible path found
    path_points = [];
    path_faces  = [];
end



function [u, v] = buildOrthonormalBasis(dir3)

    if norm(dir3) < 1e-12
        dir3 = [1 0 0];
    end
    w = dir3 / norm(dir3);

    % Choose a non-collinear reference vector
    if abs(dot(w, [1 0 0])) < 0.9
        a = [1 0 0];
    else
        a = [0 1 0];
    end

    u = cross(w, a);
    u = u / max(norm(u), 1e-12);
    v = cross(w, u);
    v = v / max(norm(v), 1e-12);
end


function dk = detourIncrement3D(P, Q, obstacle_params, u, v)


    nObs = numel(obstacle_params);
    dk = zeros(1, nObs);

    for i = 1:nObs
        c = obstacle_params(i).start + obstacle_params(i).size / 2;

        p2 = project2D(P - c, u, v);
        q2 = project2D(Q - c, u, v);

        dk(i) = windingIncrementRayCross(p2, q2);
    end
end


function p2 = project2D(x3, u, v)
% Project a 3D vector onto the 2D (u, v) coordinate system.
    p2 = [dot(x3, u), dot(x3, v)];
end


function inc = windingIncrementRayCross(p, q)

    inc = 0;

    px = p(1); py = p(2);
    qx = q(1); qy = q(2);

    % Ignore numerical jitter near the origin
    eps0 = 1e-12;
    if (abs(py) < eps0 && abs(qy) < eps0)
        return;
    end

    % Upward / downward crossing test
    if (py <= 0 && qy > 0)       % upward crossing
        if isLeft(p, q, [0 0]) > 0
            inc = +1;
        end
    elseif (py > 0 && qy <= 0)   % downward crossing
        if isLeft(p, q, [0 0]) < 0
            inc = -1;
        end
    end
end


function val = isLeft(p0, p1, p2)
% Test whether p2 lies to the left (>0) or right (<0)
% of the directed line segment p0 -> p1 in 2D.
    val = (p1(1) - p0(1)) * (p2(2) - p0(2)) - ...
          (p2(1) - p0(1)) * (p1(2) - p0(2));
end

function closest_face_idx = find_closest_face(all_faces, point)
    num_faces = length(all_faces);
    centers = zeros(num_faces, 3);
    for i = 1:num_faces
        centers(i, :) = get_center(all_faces{i});
    end

    distances = vecnorm(centers - point, 2, 2);
    min_distance = min(distances);

    threshold = 1e-6;
    closest_face_idx = find(abs(distances - min_distance) < threshold);
end


function neighbors = get_neighbors(current_idx, all_faces, visited)
    neighbors = [];
    current_face = all_faces{current_idx};

    for i = 1:length(all_faces)
        if ~visited(i) && is_adjacent(current_face, all_faces{i})
            neighbors = [neighbors, i]; %#ok<AGROW>
        end
    end
end


function is_adj = is_adjacent(face1, face2)
    vertices1 = [];
    vertices2 = [];

    for i = 2:length(face1)
        vertices1 = [vertices1; cell2mat(face1(i))]; %#ok<AGROW>
    end
    for j = 2:length(face2)
        vertices2 = [vertices2; cell2mat(face2(j))]; %#ok<AGROW>
    end

    unique_vertices1 = unique(vertices1, 'rows', 'stable');
    unique_vertices2 = unique(vertices2, 'rows', 'stable');

    tolerance = 1e-5;
    common_count = 0;

    for i = 1:size(unique_vertices1, 1)
        for j = 1:size(unique_vertices2, 1)
            if norm(unique_vertices1(i, :) - unique_vertices2(j, :)) < tolerance
                common_count = common_count + 1;
                if common_count >= 2
                    is_adj = true;
                    return;
                end
            end
        end
    end

    is_adj = false;
end


function center = get_center(face)
% Return the center of a convex block.
    center = face{1};
end



