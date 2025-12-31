function [merged_faces, removed_faces] = processConcaveFaces(faces_cells)

    concaveFaces = faces_cells{2};  % Concave face set

    % Merge all convex vertices into one array
    convexVertices = [];
    for i = 1:length(faces_cells{1})
        convexVertices = [convexVertices; cell2mat(faces_cells{1}{i})]; %#ok<AGROW>
    end
    convexVertices = unique(convexVertices, 'rows');

    % Merge all concave vertices into one array
    concaveVertices = [];
    for i = 1:length(faces_cells{2})
        concaveVertices = [concaveVertices; cell2mat(faces_cells{2}{i})]; %#ok<AGROW>
    end
    concaveVertices = unique(concaveVertices, 'rows');

    merged_faces  = {};  % Store generated convex hull face sets
    removed_faces = {};  % Store removed concave faces

    % Iterate over each concave face
    for i = 1:length(concaveFaces)
        faceVertices = cell2mat(concaveFaces{i});  % Convert current face vertices to a matrix

        % Compute face normal vector (using the first three vertices)
        normalVector = cross(faceVertices(2,:) - faceVertices(1,:), ...
                             faceVertices(3,:) - faceVertices(1,:));

        nrm = norm(normalVector);
        if nrm < 1e-12
            % Degenerate face -> mark as removed
            removed_faces{end+1} = faceVertices; %#ok<AGROW>
            continue;
        end
        normalVector = normalVector / nrm;

        % Face centroid
        faceCenter = mean(faceVertices, 1);

        % Classify vertices by the sign of projection onto the normal direction
        [posConvexSide, negConvexSide, posConcaveSide, negConcaveSide] = ...
            classifyVerticesBySide(faceCenter, convexVertices, concaveVertices, normalVector);

        % One side contains only convex vertices, while the other side contains concave vertices
        if isempty(posConcaveSide) && ~isempty(posConvexSide) && ~isempty(negConcaveSide)
            % Combine current face vertices and positive-side convex vertices
            combinedVertices = [faceVertices; posConvexSide];
        elseif isempty(negConcaveSide) && ~isempty(negConvexSide) && ~isempty(posConcaveSide)
            % Combine current face vertices and negative-side convex vertices
            combinedVertices = [faceVertices; negConvexSide];
        else
            % Condition 2:
            % If neither condition is met, remove the face
            removed_faces{end+1} = faceVertices; %#ok<AGROW>
            continue;
        end

        % Deduplicate merged vertices and compute convex hull
        combinedVertices = unique(combinedVertices, 'rows');
        try
            K = convhulln(combinedVertices);
        catch
            removed_faces{end+1} = faceVertices; %#ok<AGROW>
            continue;
        end

        % Store all faces of the convex hull for the current concave face
        current_convex_hull = {};
        for j = 1:size(K, 1)
            face = combinedVertices(K(j, :), :);
            current_convex_hull{end+1} = face; %#ok<AGROW>
            % plotFace(face); % Enable if visualization is needed
        end

        % Append the current convex hull face set
        merged_faces{end+1} = current_convex_hull; %#ok<AGROW>
    end
end

function [posConvexSide, negConvexSide, posConcaveSide, negConcaveSide] = ...
    classifyVerticesBySide(faceCenter, convexVertices, concaveVertices, normalVector)

    posConvexSide  = [];
    negConvexSide  = [];
    posConcaveSide = [];
    negConcaveSide = [];

    % Classify convex vertices
    for j = 1:size(convexVertices, 1)
        vertex = convexVertices(j, :);
        projection = dot(vertex - faceCenter, normalVector);
        if projection > 0
            posConvexSide = [posConvexSide; vertex]; %#ok<AGROW>
        elseif projection < 0
            negConvexSide = [negConvexSide; vertex]; %#ok<AGROW>
        end
    end

    % Classify concave vertices
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

function plotFace(face)

    fill3(face(:,1), face(:,2), face(:,3), [0.6, 0.8, 1], ...
          'FaceAlpha', 0.5, 'EdgeColor', 'k');
    hold on;
    axis equal;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
end
