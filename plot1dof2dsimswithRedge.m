load('something\workspace_fittingdata.mat')
path = 'something';
PHIS = make_1d_cutouts_from_phi(phi_dyn.phisq{1});
phix = PHIS{1}; phiy = PHIS{2}';
geom = phi_dyn.return_geometry();
xx = geom.X(1,:);
save([path 'workspace_fittingdata'],'PHIS', 'phix', 'phiy', 'xx', 'path', 'geom', '-append')
Redge = (16*S)^(1/4);
dx = geom.dx;
area = dx*sum(phix);
phixN = phix/area;
save([path 'workspace_fittingdata'],'Redge', 'dx', 'area', 'phixN', '-append')

plot(xx,phixN)
%labels & limits
xlabel('x'); ylabel('|\phi|_x^2');
xlim([-1 1]); ylim([0 max(phixN)*1.1]);
%drawing Redge lines
hold on;
line([-Redge -Redge], [0 max(phixN)*1.1], 'Color','red','LineStyle','--'); line([Redge Redge], [0 max(phixN)*1.1],'Color','red','LineStyle','--');
%legend
legend(sprintf('S=1e-7,N=257^2\ngrid(-%d,%d)',max(xx),max(xx)), sprintf('R_{edge} (dx=%.4f)',dx));
