%% Making Surface Plots From Scatter Data
% How do you turn a collection of XYZ triplets into a surface plot? This is
% the most frequently asked 3D plotting question that I got when I was in
% Tech Support.

%%
% The problem is that the data is made up of individual (x,y,z)
% measurements. It isn't laid out on a rectilinear grid, which is what the
% SURF command expects. A simple plot command isn't very useful.
data = load("data.txt");
plot3(data(:,1),data(:,3),data(:,2),'.-')

%% Little triangles
% The solution is to use Delaunay triangulation. Let's look at some
% info about the "tri" variable.

tri = delaunay(data(:,1),data(:,3));
plot(data(:,1),data(:,3),'.')

%%
% How many triangles are there?

[r,c] = size(tri);
disp(r)

%% Plot it with TRISURF

h = trisurf(tri, data(:,1),data(:,3),data(:,2));
axis vis3d

%% Clean it up

axis off
l = light('Position',[-50 -15 29])
set(gca,'CameraPosition',[208 -50 7687])
lighting phong
shading interp
colorbar EastOutside
