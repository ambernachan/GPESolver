% get phi, xx, yy, zz (meshgrids)

selection = phi(:);
indices = find(selection < max(selection)/1000);
selection(indices) = [];

markers = selection ./ max(selection) .* 40; 
markers = ceil(markers) + 1;

xx = xx(:);
xx(indices) = [];
yy = yy(:);
yy(indices) = [];
zz = zz(:);
zz(indices) = [];

fig = scatter3(xx, yy, zz, markers, selection, 'filled');
xlabel('x'); ylabel('y'); zlabel('z');
xlim([-10 10]); ylim([-10 10]); zlim([-10 10]);
fig.MarkerFaceAlpha = 0.5;
