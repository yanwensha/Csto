clear; close all; clc
program_start = tic;  % Total runtime timer

%% ========== Parameters ==========
total_size   = 1.0;
n_divisions  = 10;
cube_size    = total_size / n_divisions;

start_point  = [0 0 0];
end_point    = [1 1 1];

%% ========== Space partition (only once) ==========
all_faces = createCubeCells(n_divisions, cube_size);

%% ========== Random obstacle environment (one-shot) ==========

obstacle_params = obstacle_data();
num_obstacles   = length(obstacle_params);
obstacle_faces  = cell(1, num_obstacles);

for i = 1:num_obstacles
    obstacle_faces{i} = createObstacleFaces( ...
        obstacle_params(i).start, obstacle_params(i).size);
end

%% ========== Intersection detection + convex decomposition ==========
[new_faces1, intersection_points] = ...
    Intersection_detection(all_faces, obstacle_faces, obstacle_params);

new_faces2 = convexDecompose(new_faces1, intersection_points, obstacle_params);

%% ========== Chain search (detour cost ON) ==========
cfg.useDetourCost = true;

[path_pts, path_faces] = ...
    chainpath(new_faces2, start_point, end_point, obstacle_params, cfg);

if isempty(path_faces)
    error('No feasible chain found (detour cost ON).');
end

% Build the solved polytope chain as a cell array (for plotting)
poly_chain = cell(1, length(path_faces));
for i = 1:length(path_faces)
    poly_chain{i} = new_faces2{path_faces(i)};
end

%% ========== Convex optimization on the chain ==========
[optimal_path, final_length, t_opt] = ...
    optimize_on_chain(new_faces2, path_faces, start_point, end_point);
total_time = toc(program_start);

%% ========== Output ==========
fprintf('\nFinal optimized path length (with detour cost): %.6f\n', final_length);
fprintf('Optimization time: %.4f s\n', t_opt);
fprintf('Total runtime: %.4f s\n', total_time);

%% ========== Visualization: polytope chain + obstacles + optimized path ==========
figure('Name','Polytope Chain + Obstacles + Optimized Path');
hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Solved Polytope Chain and Optimized Path (Detour Cost ON)');
view(3);

% Plot polytope chain (wireframe)
plot_faces(poly_chain);

% Plot obstacles (semi-transparent)
drawCubeObstacles(obstacle_params, [1 0 0], 0.25);

% Plot optimized path
plot3(optimal_path(:,1), optimal_path(:,2), optimal_path(:,3), '-o', ...
      'LineWidth', 2, 'MarkerSize', 5);
scatter3(start_point(1), start_point(2), start_point(3), 80, 'filled');
scatter3(end_point(1), end_point(2), end_point(3), 80, 'filled');

legend({'Polytope chain (wireframe)','Obstacles','Optimized path','Start','Goal'}, ...
       'Location','best');

hold off;

%% ========================================================================
%  Local functions
% ========================================================================

function cells = createCubeCells(n_divisions, cube_size)
    % Create small cube cells (center + 6 faces), stored in a cell array.
    cells = {};

    for i = 0:(n_divisions - 1)
        for j = 0:(n_divisions - 1)
            for k = 0:(n_divisions - 1)

                center = [(i + 0.5) * cube_size, (j + 0.5) * cube_size, (k + 0.5) * cube_size];

                x_min = i * cube_size;       x_max = (i + 1) * cube_size;
                y_min = j * cube_size;       y_max = (j + 1) * cube_size;
                z_min = k * cube_size;       z_max = (k + 1) * cube_size;

                faces = {
                    [x_min, y_min, z_min; x_max, y_min, z_min; x_max, y_max, z_min; x_min, y_max, z_min], % bottom
                    [x_min, y_min, z_max; x_max, y_min, z_max; x_max, y_max, z_max; x_min, y_max, z_max], % top
                    [x_min, y_min, z_min; x_min, y_min, z_max; x_min, y_max, z_max; x_min, y_max, z_min], % left
                    [x_max, y_min, z_min; x_max, y_min, z_max; x_max, y_max, z_max; x_max, y_max, z_min], % right
                    [x_min, y_min, z_min; x_max, y_min, z_min; x_max, y_min, z_max; x_min, y_min, z_max], % front
                    [x_min, y_max, z_min; x_max, y_max, z_min; x_max, y_max, z_max; x_min, y_max, z_max]  % back
                };

                cube_data = cell(1, 7);
                cube_data{1} = center;
                for f = 1:6
                    cube_data{f+1} = faces{f};
                end

                cells{end+1} = cube_data;
            end
        end
    end
end

function obstacle_faces = createObstacleFaces(start_point, size)
    % Create faces of one axis-aligned box obstacle (6 quad faces).
    x_min = start_point(1); x_max = start_point(1) + size(1);
    y_min = start_point(2); y_max = start_point(2) + size(2);
    z_min = start_point(3); z_max = start_point(3) + size(3);

    obstacle_faces = {
        [x_min, y_min, z_min; x_max, y_min, z_min; x_max, y_max, z_min; x_min, y_max, z_min], % bottom
        [x_min, y_min, z_max; x_max, y_min, z_max; x_max, y_max, z_max; x_min, y_max, z_max], % top
        [x_min, y_min, z_min; x_min, y_min, z_max; x_min, y_max, z_max; x_min, y_max, z_min], % left
        [x_max, y_min, z_min; x_max, y_min, z_max; x_max, y_max, z_max; x_max, y_max, z_min], % right
        [x_min, y_min, z_min; x_max, y_min, z_min; x_max, y_min, z_max; x_min, y_min, z_max], % front
        [x_min, y_max, z_min; x_max, y_max, z_min; x_max, y_max, z_max; x_min, y_max, z_max]  % back
    };
end

function L = compute_path_length(P)
    % Polyline length.
    L = 0;
    for i = 1:size(P,1)-1
        L = L + norm(P(i+1,:) - P(i,:));
    end
end

function [optimal_path, total_distance, t_opt] = ...
    optimize_on_chain(new_faces2, path_faces, start_point, end_point)

    % Optimize waypoint positions constrained to the polytope chain.
    tic;

    % Extract polytope chain (subset of cells)
    new_faces3 = cell(1, length(path_faces));
    for i = 1:length(path_faces)
        new_faces3{i} = new_faces2{path_faces(i)};
    end

    % This function must be provided by your project
    point_constraints = extractPointConstraints(new_faces3);

    % Build optimization problem
    P = optimvar('P', length(point_constraints)+2, 3);
    prob = optimproblem;

    prob.Constraints.start = P(1,:) == start_point;
    prob.Constraints.goal  = P(end,:) == end_point;

    for i = 1:length(point_constraints)
        geom = point_constraints{i};
        k = size(geom,1);

        alpha = optimvar(sprintf('a_%d', i), k, ...
            'LowerBound', 0, 'UpperBound', 1);

        prob.Constraints.(sprintf('c_%d', i)) = P(i+1,:) == alpha' * geom;
        prob.Constraints.(sprintf('w_%d', i)) = sum(alpha) == 1;
    end

    % Objective: sum of squared segment lengths (keeps your original logic)
    prob.Objective = sum(sum((P(2:end,:) - P(1:end-1,:)).^2, 2));

    sol = solve(prob);

    optimal_path   = sol.P;
    total_distance = compute_path_length(optimal_path);
    t_opt          = toc;
end

function plot_faces(cells_poly)
    % Plot a set of convex polytopes as wireframes.
    % Each poly: {center, face1, face2, ..., face6?} but can be variable faces.

    for i = 1:length(cells_poly)
        poly = cells_poly{i};
        if isempty(poly), continue; end

        for j = 2:length(poly)  % skip the center
            V = poly{j};
            if size(V,1) < 3, continue; end

            % Close the loop
            plot3(V([1:end, 1], 1), V([1:end, 1], 2), V([1:end, 1], 3), ...
                  'k-', 'LineWidth', 1);
        end
    end
end

function drawCubeObstacles(obstacle_params, color, alphaVal)
    % Draw all axis-aligned box obstacles with patch.
    for i = 1:length(obstacle_params)
        origin = obstacle_params(i).start;
        sz     = obstacle_params(i).size;
        drawCube(origin, sz, color, alphaVal);
    end
end

function drawCube(origin, sz, color, alphaVal)
    % origin: [x,y,z]
    % sz:     [lx,ly,lz]
    [X, Y, Z] = ndgrid([0 1]);

    pts = [X(:), Y(:), Z(:)] .* sz + origin;

    faces = [
        1 2 4 3;
        1 2 6 5;
        2 4 8 6;
        3 4 8 7;
        1 3 7 5;
        5 6 8 7
    ];

    patch('Vertices', pts, 'Faces', faces, ...
          'FaceColor', color, 'FaceAlpha', alphaVal, 'EdgeColor', 'k');
end

