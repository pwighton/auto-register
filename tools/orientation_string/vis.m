function vis(R)

   planes.vertices = [ 0, -1,  1; % sag
                       0,  1,  1;
                       0,  1, -1;
                       0, -1, -1;
                      -1,  0,  1; % cor
                       1,  0,  1;
                       1,  0, -1;
                      -1,  0, -1;
                      -1,  1,  0; % tra
                       1,  1,  0;
                       1, -1,  0;
                      -1, -1,  0;
                       ];

   planes.faces = [  1,  2,  3; % sag
                     1,  3,  4;
                     5,  6,  7; % cord
                     5,  7,  8;
                     9, 10, 11; % tra
                     9, 11, 12;
                  ];

   planes.cdata = [linspace(1, length(planes.faces), length(planes.faces) - 1), 2 * length(planes.faces)]';

   show_planes(planes); hold on;

   rot_planes = planes;
   rot_planes.vertices = (R * rot_planes.vertices')';

   rot_planes.cdata = [1, linspace(length(planes.faces) + 1, 2 * length(planes.faces), length(planes.faces) - 1)]';

   show_planes(rot_planes);

end

function show_planes(planes)

  patch('vertices', planes.vertices, ...
        'faces', planes.faces, ...
        'facevertexcdata', planes.cdata, ...
        'facecolor', 'flat', ...
        'edgecolor', 'none');

  axis on;
  xlabel('x');
  ylabel('y');
  zlabel('z');

  axis vis3d;
  axis equal;
  lighting gouraud
  material dull
  camlight;
  cameratoolbar;
end
