function obstacle_params = obstacle_data()

    numObstacles   = 70;
    cubeSize       = 0.1;     % Edge length of each cube
    minCenterDist  = 0.213;   % Minimum distance between obstacle centers
    spaceMin       = 0.0;
    spaceMax       = 1.0;

    % Preallocate obstacle structure array
    obstacle_params(numObstacles) = struct('start', [], 'size', []);

    % Maximum number of random trials for placing each obstacle
    % (larger value is needed due to the strict distance constraints)
    maxTrials = 200000;
    count = 0;

    minCenterDist2 = minCenterDist^2;

    for i = 1:numObstacles
        success = false;
        trial = 0;

        while ~success
            trial = trial + 1;
            if trial > maxTrials
                error(['Failed to generate obstacle %d: ' ...
                       'constraints are too strict or space is too crowded ' ...
                       '(maxTrials = %d).'], i, maxTrials);
            end

            % Random start position in continuous space (stay within bounds)
            start = spaceMin + (spaceMax - cubeSize - spaceMin) * rand(1,3);
            minP  = start;
            maxP  = start + cubeSize;

            % Center of the new cube
            cNew = start + 0.5 * cubeSize;

            % Check against existing obstacles:
            % (1) no AABB intersection
            % (2) center distance larger than minCenterDist
            ok = true;
            for j = 1:count
                otherMin = obstacle_params(j).start;
                otherMax = otherMin + obstacle_params(j).size;
                cOld     = otherMin + 0.5 * obstacle_params(j).size;

                % Condition 1: AABB non-intersection
                if aabbIntersect(minP, maxP, otherMin, otherMax)
                    ok = false;
                    break;
                end

                % Condition 2: center-to-center distance constraint
                d = cNew - cOld;
                if (d * d.') <= minCenterDist2   % Equivalent to norm(d)^2 <= minCenterDist^2
                    ok = false;
                    break;
                end
            end

            if ok
                count = count + 1;
                obstacle_params(i).start = start;
                obstacle_params(i).size  = [cubeSize, cubeSize, cubeSize];
                success = true;
            end
        end
    end
end

function flag = aabbIntersect(minA, maxA, minB, maxB)
    % Return true if two axis-aligned bounding boxes intersect.
    % If they are separated along any axis, they do not intersect.

    flag = ~( ...
        maxA(1) <= minB(1) || minA(1) >= maxB(1) || ...
        maxA(2) <= minB(2) || minA(2) >= maxB(2) || ...
        maxA(3) <= minB(3) || minA(3) >= maxB(3) );
end


