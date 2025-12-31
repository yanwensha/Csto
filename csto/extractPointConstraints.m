function point_constraints = extractPointConstraints(new_faces3)

    num_faces = length(new_faces3);
    if num_faces < 2
        point_constraints = {};
        return;
    end

    point_constraints = cell(num_faces - 1, 1);

    for i = 1:num_faces-1
        poly1 = new_faces3{i};
        poly2 = new_faces3{i+1};

        [shared_faces, shared_lines] = findSharedGeometryByVertices(poly1, poly2);

        if ~isempty(shared_faces)
            % Shared face: output an ordered polygon (>=3 points)
            point_constraints{i} = orderCoplanarPolygon(shared_faces{1});
        elseif ~isempty(shared_lines)
            % Shared edge: output the two endpoints
            point_constraints{i} = shared_lines{1};
        else
            error(['Invalid use of extractPointConstraints: ' ...
                   'No shared geometry between path polytopes %d and %d; cannot build constraints.'], i, i+1);
        end
    end
end


function [shared_faces, shared_lines] = findSharedGeometryByVertices(poly1, poly2)

    faces1 = poly1(2:end);
    faces2 = poly2(2:end);

    % Key: a slightly looser tolerance than 1e-8 for robustness
    tol = 1e-6;

    shared_faces = {};
    shared_lines = {};

    % Tolerance-based unique for each face
    faces1u = cell(size(faces1));
    faces2u = cell(size(faces2));
    for i = 1:numel(faces1)
        faces1u{i} = unique_tol(faces1{i}, tol);
    end
    for j = 1:numel(faces2)
        faces2u{j} = unique_tol(faces2{j}, tol);
    end

    % Choose the face with the maximum number of common points as the shared face
    % (to reduce false matches)
    bestFace = [];
    bestK = 0;
    bestLine = [];

    for i = 1:numel(faces1u)
        V1 = faces1u{i};
        for j = 1:numel(faces2u)
            V2 = faces2u{j};

            common = intersect_tol(V1, V2, tol);
            k = size(common, 1);

            if k >= 3
                if k > bestK
                    bestK = k;
                    bestFace = common;
                end
            elseif k == 2
                if isempty(bestLine)
                    bestLine = common;
                end
            end
        end
    end

    if ~isempty(bestFace)
        shared_faces{1} = unique_tol(bestFace, tol);
        return;
    end

    if ~isempty(bestLine)
        shared_lines{1} = unique_tol(bestLine, tol);
    end
end


function common = intersect_tol(A, B, tol)
% Tolerance-based row intersection:
% return matched points between A and B (via quantization + ismember).
    if isempty(A) || isempty(B)
        common = zeros(0,3);
        return;
    end

    Aq = round(A./tol).*tol;
    Bq = round(B./tol).*tol;

    [lia, ~] = ismember(Aq, Bq, 'rows');
    common = A(lia, :);

    common = unique_tol(common, tol);
end


function P = unique_tol(P, tol)
% Tolerance-based unique (more robust than direct unique).
    if isempty(P), return; end
    Q = round(P./tol).*tol;
    [~, ia] = unique(Q, 'rows', 'stable');
    P = P(ia,:);
end


function ordered = orderCoplanarPolygon(vertices)

    tol = 1e-6;
    vertices = unique_tol(vertices, tol);

    if size(vertices,1) <= 2
        ordered = vertices;
        return;
    end

    % If exactly 3 points, ordering is not critical
    if size(vertices,1) == 3
        ordered = vertices;
        return;
    end

    % Estimate a normal by finding three non-collinear points
    p1 = vertices(1,:);
    n = [0 0 0];
    for k = 2:size(vertices,1)-1
        p2 = vertices(k,:);
        p3 = vertices(k+1,:);
        n0 = cross(p2 - p1, p3 - p1);
        if norm(n0) > 1e-12
            n = n0;
            break;
        end
    end

    if norm(n) < 1e-12
        % Nearly collinear: return as-is (edge-like constraint is more reasonable)
        ordered = vertices;
        return;
    end
    n = n / norm(n);

    % Project onto plane by dropping the largest normal component
    [~, idx] = max(abs(n));
    switch idx
        case 1 % drop x -> use yz
            P2 = vertices(:, [2 3]);
        case 2 % drop y -> use xz
            P2 = vertices(:, [1 3]);
        otherwise % drop z -> use xy
            P2 = vertices(:, [1 2]);
    end

    % 2D convex hull
    try
        k = convhull(P2(:,1), P2(:,2));
    catch
        ordered = vertices;
        return;
    end
    ordered = vertices(k, :);

    % convhull returns a closed loop; remove duplicated endpoint
    if size(ordered,1) >= 2 && all(abs(ordered(1,:) - ordered(end,:)) < 1e-12)
        ordered(end,:) = [];
    end
end



