function [new_faces1, intersection_points] = Intersection_detection(all_faces, obstacle_faces, obstacle_params)

    tol = 1e-9;

    new_faces1 = cell(1, numel(all_faces));
    intersection_points = cell(1, numel(all_faces));

    for i = 1:numel(all_faces)

        % Generator point of current cell (your original logic: all_faces{i}{1}(1,:))
        gen_point = all_faces{i}{1}(1, :); %#ok<NASGU>

        % Keep the cell by default
        new_faces1{i} = all_faces{i};
        all_intersections = [];

        % Compute AABB of current cell from all its vertices
        [cellMin, cellMax] = getCellAABB(all_faces{i});

        % Iterate over obstacles
        removeCell = false;
        for obs_index = 1:numel(obstacle_params)

            obs_start = obstacle_params(obs_index).start;
            obs_size  = obstacle_params(obs_index).size;
            obsMin = obs_start;
            obsMax = obs_start + obs_size;

            % (1) If the entire cell is fully contained in the obstacle: remove the cell
            if all(cellMin >= obsMin - tol) && all(cellMax <= obsMax + tol)
                new_faces1{i} = {};
                removeCell = true;
                break;
            end

            % (2) If AABBs do not intersect: no face intersection is possible -> skip (fast)
            if ~aabbIntersect(cellMin, cellMax, obsMin, obsMax)
                continue;
            end

            % (3) AABBs intersect: perform face-face intersection tests
            for j = 2:numel(all_faces{i})      % j=2..end are the faces of the cell (your setting)
                faceA = all_faces{i}{j};       % 4x3 quad
                for k = 1:numel(obstacle_faces{obs_index})
                    faceB = obstacle_faces{obs_index}{k};  % 4x3 quad

                    pts = quadQuadIntersectionPoints(faceA, faceB, tol);
                    if ~isempty(pts)
                        all_intersections = [all_intersections; pts]; %#ok<AGROW>
                    end
                end
            end
        end

        if ~removeCell
            intersection_points{i} = uniquePoints(all_intersections, 1e-8);
        else
            intersection_points{i} = [];
        end
    end

    % Remove empty cells in new_faces1 while keeping intersection_points aligned
    keep = ~cellfun(@isempty, new_faces1);
    new_faces1 = new_faces1(keep);
    intersection_points = intersection_points(keep);
end

%% ---------- AABB utilities ----------
function [minP, maxP] = getCellAABB(cellFaces)
% cellFaces: {genInfo, face1, face2, ...}
    V = [];
    for t = 2:numel(cellFaces)
        V = [V; cellFaces{t}]; %#ok<AGROW>
    end
    minP = min(V, [], 1);
    maxP = max(V, [], 1);
end

function flag = aabbIntersect(minA, maxA, minB, maxB)
% Axis-aligned bounding box intersection:
% If separated along any axis, they do not intersect.
    flag = ~( maxA(1) <= minB(1) || minA(1) >= maxB(1) || ...
              maxA(2) <= minB(2) || minA(2) >= maxB(2) || ...
              maxA(3) <= minB(3) || minA(3) >= maxB(3) );
end

%% ---------- Quad-Quad intersection (numeric, robust) ----------
function pts = quadQuadIntersectionPoints(a1, b1, tol)

    pts = [];

    % Plane normals
    A=a1(1,:); B=a1(2,:); C=a1(3,:); D=a1(4,:); %#ok<NASGU>
    E=b1(1,:); F=b1(2,:); G=b1(3,:); H=b1(4,:); %#ok<NASGU>

    n1 = cross(B-A, C-A);
    n2 = cross(F-E, G-E);

    n1n = norm(n1); n2n = norm(n2);
    if n1n < tol || n2n < tol
        return;
    end
    n1 = n1 / n1n;
    n2 = n2 / n2n;

    d = cross(n1, n2);
    dn = norm(d);

    % Parallel case
    if dn < 1e-10
        % Check coplanarity: distance from E to plane(a1)
        dist = abs(dot(n1, E - A));
        if dist > 1e-8
            return; % parallel but not coplanar: no intersection
        end
        % Coplanar: compute 2D polygon intersection
        pts = coplanarQuadIntersection(a1, b1, tol);
        return;
    end

    d = d / dn;

    rhs = [dot(n1, A); dot(n2, E); 0];
    M = [n1; n2; d];
    if abs(det(M)) < 1e-12
        % Degenerate: use least squares
        p0 = (M\rhs).';
    else
        p0 = (M\rhs).';
    end

    % Intersect the line with quad edges (numerical line-segment intersection)
    cand = [];
    cand = [cand; lineSegmentIntersections(p0, d, a1, tol)]; %#ok<AGROW>
    cand = [cand; lineSegmentIntersections(p0, d, b1, tol)]; %#ok<AGROW>

    if isempty(cand)
        return;
    end

    % Deduplicate + filter: must lie inside both quads
    cand = uniquePoints(cand, 1e-8);

    valid = [];
    for i = 1:size(cand,1)
        P = cand(i,:);
        if pointInQuad3D(P, a1, tol) && pointInQuad3D(P, b1, tol)
            valid = [valid; P]; %#ok<AGROW>
        end
    end

    pts = uniquePoints(valid, 1e-8);
end

function cand = lineSegmentIntersections(p0, d, quad, tol)
% Compute intersections between line L(s)=p0+s*d and the 4 edges of a quad.
    cand = [];
    for i = 1:4
        a = quad(i,:);
        b = quad(mod(i,4)+1,:);
        e = (b - a);

        % Solve p0 + s d = a + t e
        % => [d, -e] [s;t] = a - p0
        A = [d(:), -e(:)];
        rhs = (a - p0).';
        if rcond(A.'*A) < 1e-14
            continue;
        end
        st = (A.'*A)\(A.'*rhs);  % least squares
        s = st(1); t = st(2);

        P = p0 + s*d;

        % Residual check (ensure it is a valid intersection)
        Q = a + t*e;
        if norm(P - Q) > 1e-7
            continue;
        end

        % Segment parameter t in [0,1]
        if t >= -1e-8 && t <= 1+1e-8
            cand = [cand; 0.5*(P+Q)]; %#ok<AGROW>
        end
    end
end

%% ---------- Coplanar quad intersection via 2D projection ----------
function pts3 = coplanarQuadIntersection(a1, b1, tol)

    pts3 = [];

    % Build local coordinate system (u,v) on the plane of a1
    A=a1(1,:); B=a1(2,:); C=a1(3,:);
    n = cross(B-A, C-A);
    nn = norm(n);
    if nn < tol, return; end
    n = n/nn;

    u = (B-A); u = u / max(norm(u), tol);
    v = cross(n, u); v = v / max(norm(v), tol);

    % 3D -> 2D
    P1 = project2D(a1, A, u, v);
    P2 = project2D(b1, A, u, v);

    % polyshape intersection (requires ordered 2D polygons)
    try
        ps1 = polyshape(P1(:,1), P1(:,2));
        ps2 = polyshape(P2(:,1), P2(:,2));
        psi = intersect(ps1, ps2);
        if isempty(psi.Vertices)
            return;
        end
        V2 = psi.Vertices; % Nx2
    catch
        % If polyshape is unavailable or degeneracy occurs, return empty
        return;
    end

    % 2D -> 3D
    for i = 1:size(V2,1)
        pts3 = [pts3; A + V2(i,1)*u + V2(i,2)*v]; %#ok<AGROW>
    end
    pts3 = uniquePoints(pts3, 1e-8);
end

function P2 = project2D(P3, origin, u, v)
    P2 = zeros(size(P3,1), 2);
    for i = 1:size(P3,1)
        w = P3(i,:) - origin;
        P2(i,1) = dot(w, u);
        P2(i,2) = dot(w, v);
    end
end

%% ---------- Point in quad (3D) ----------
function inside = pointInQuad3D(P, quad, tol)
% Split quad into two triangles and test using barycentric coordinates.
    A = quad(1,:); B = quad(2,:); C = quad(3,:); D = quad(4,:);
    inside = pointInTri3D(P, A, B, C, tol) || pointInTri3D(P, A, C, D, tol);
end

function inside = pointInTri3D(P, A, B, C, tol)
% Use barycentric coordinates after projecting to the plane with maximum normal component.
    n = cross(B-A, C-A);
    [~, idx] = max(abs(n));

    % Choose projection plane
    switch idx
        case 1 % drop x -> use yz
            P2 = P([2 3]); A2=A([2 3]); B2=B([2 3]); C2=C([2 3]);
        case 2 % drop y -> use xz
            P2 = P([1 3]); A2=A([1 3]); B2=B([1 3]); C2=C([1 3]);
        otherwise % drop z -> use xy
            P2 = P([1 2]); A2=A([1 2]); B2=B([1 2]); C2=C([1 2]);
    end

    v0 = C2 - A2;
    v1 = B2 - A2;
    v2 = P2 - A2;

    den = v0(1)*v1(2) - v1(1)*v0(2);
    if abs(den) < tol
        inside = false;
        return;
    end
    a = (v2(1)*v1(2) - v1(1)*v2(2)) / den;
    b = (v0(1)*v2(2) - v2(1)*v0(2)) / den;

    inside = (a >= -1e-8) && (b >= -1e-8) && (a + b <= 1 + 1e-8);
end

%% ---------- Unique points with tolerance ----------
function U = uniquePoints(P, epsv)
    if isempty(P)
        U = P;
        return;
    end
    % Quantize by epsv for deduplication
    Q = round(P/epsv)*epsv;
    [~, ia] = unique(Q, 'rows', 'stable');
    U = P(ia,:);
end


