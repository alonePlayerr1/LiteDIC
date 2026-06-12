% Demo: geodesic triangulation of the upper hemisphere.
% Icosahedron -> each triangle subdivided into 4 -> new vertices projected
% onto the sphere (2 subdivision levels). See icosahedron_sphere.m.
[vertices, faces] = icosahedron_sphere(1, [0, 2*pi], [0, pi/2], 2);
