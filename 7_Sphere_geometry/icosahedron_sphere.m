% Geodesic triangulation of a (partial) sphere from a subdivided icosahedron.
function [vertices, faces] = icosahedron_sphere(R, theta_range, phi_range, subdivisions)
    % R - sphere radius
    % theta_range - azimuth range, rad, e.g. [0, 2*pi]
    % phi_range - polar angle range (from +Z), rad, e.g. [0, pi/2]
    % subdivisions - number of subdivision iterations

    % Golden ratio
    t = (1 + sqrt(5)) / 2;

    % Initial icosahedron vertices
    vertices = [
        -1,  t,  0;  1,  t,  0; -1, -t,  0;  1, -t,  0;
         0, -1,  t;  0,  1,  t;  0, -1, -t;  0,  1, -t;
         t,  0, -1;  t,  0,  1; -t,  0, -1; -t,  0,  1
    ];

    % Normalize vertices onto the unit sphere
    vertices = vertices ./ sqrt(sum(vertices.^2, 2));

    % Initial icosahedron faces
    faces = [
        1, 12, 6; 1, 6, 2; 1, 2, 8; 1, 8, 11; 1, 11, 12;
        2, 6, 10; 6, 12, 5; 12, 11, 3; 11, 8, 7; 8, 2, 9;
        4, 10, 5; 4, 5, 3; 4, 3, 7; 4, 7, 9; 4, 9, 10;
        5, 10, 6; 3, 5, 12; 7, 3, 11; 9, 7, 8; 10, 9, 2
    ];

    % Subdivision: each triangle splits into 4, new vertices are projected
    % back onto the sphere
    for s = 1:subdivisions
        new_faces = [];
        edge_midpoints = containers.Map;

        for f = 1:size(faces, 1)
            v1 = faces(f, 1); v2 = faces(f, 2); v3 = faces(f, 3);

            % Midpoints of the three edges (shared midpoints are reused
            % via the map, so the mesh stays watertight)
            e1 = sort([v1, v2]); e2 = sort([v2, v3]); e3 = sort([v3, v1]);
            e1_key = sprintf('%d-%d', e1(1), e1(2));
            e2_key = sprintf('%d-%d', e2(1), e2(2));
            e3_key = sprintf('%d-%d', e3(1), e3(2));

            % Midpoint of e1
            if ~isKey(edge_midpoints, e1_key)
                mid = mean(vertices([v1, v2], :));
                mid = mid / norm(mid); % project onto the sphere
                vertices = [vertices; mid];
                edge_midpoints(e1_key) = size(vertices, 1);
            end
            m1 = edge_midpoints(e1_key);

            % Midpoint of e2
            if ~isKey(edge_midpoints, e2_key)
                mid = mean(vertices([v2, v3], :));
                mid = mid / norm(mid);
                vertices = [vertices; mid];
                edge_midpoints(e2_key) = size(vertices, 1);
            end
            m2 = edge_midpoints(e2_key);

            % Midpoint of e3
            if ~isKey(edge_midpoints, e3_key)
                mid = mean(vertices([v3, v1], :));
                mid = mid / norm(mid);
                vertices = [vertices; mid];
                edge_midpoints(e3_key) = size(vertices, 1);
            end
            m3 = edge_midpoints(e3_key);

            % Four new triangles
            new_faces = [new_faces;
                v1, m1, m3;
                v2, m2, m1;
                v3, m3, m2;
                m1, m2, m3];
        end
        faces = new_faces;
    end

    % Scale vertices to radius R
    vertices = vertices * R;

    % Filter vertices by the angular ranges
    [theta, elev, ~] = cart2sph(vertices(:,1), vertices(:,2), vertices(:,3));
    theta = mod(theta, 2*pi); % wrap theta to [0, 2*pi]
    phi = pi/2 - elev;        % cart2sph returns elevation; convert to polar angle from +Z

    % Keep vertices inside the requested ranges
    valid_vertices = (theta >= theta_range(1)) & (theta <= theta_range(2)) & ...
                     (phi >= phi_range(1)) & (phi <= phi_range(2));

    valid_idx = find(valid_vertices);
    vertex_map = zeros(size(vertices, 1), 1);
    vertex_map(valid_vertices) = 1:sum(valid_vertices);

    vertices = vertices(valid_vertices, :);

    % Filter faces
    valid_faces = false(size(faces, 1), 1);
    for f = 1:size(faces, 1)
        if all(valid_vertices(faces(f, :)))
            valid_faces(f) = true;
        end
    end
    faces = faces(valid_faces, :);

    % Reindex faces
    faces = vertex_map(faces);

    % Drop any faces with zero indices (if present)
    faces = faces(all(faces > 0, 2), :);

    % Visualization
    figure;
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'c', 'EdgeColor', 'k');
    axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Icosahedron-based triangulation of a partial sphere');
    grid on;
end

% Usage example:
% [vertices, faces] = icosahedron_sphere(1, [0, 2*pi], [0, pi/2], 2);
