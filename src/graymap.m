% Creates a gray scale colormap
% and sets the system colormap to it.
map = [0:255;0:255;0:255];
map = map';
map = map/255;
colormap(map);