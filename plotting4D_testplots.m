
clear all; 

% I have: 
% r1 = [1 2 3 4 5 6]; 
% p  = [2 4 5 8 9 7]; 
% r2 = [10 45 1 0 7 9]; 
% v  = 5 .*r1 .^2 + 2 .*p .^2 + r2 .^2;

%% Method 1

func = @(x,y,z)(sin(2*x.^2+2*y.^2).*cos(2*z.^2));
[x,y,z] = meshgrid(-1:0.1:1, -1:0.1:1, -1:0.1:1);
v = func(x,y,z);
figure
p1 = patch(isosurface(x,y,z,v,0.5));
hold on
p2 = patch(isosurface(x,y,z,v,0.8));
isonormals(x,y,z,v,p1);
isonormals(x,y,z,v,p2);
p1.FaceColor = 'red';
p2.FaceColor = 'green';
p1.EdgeColor = 'none';
p2.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud
alpha(0.3);

%% Method 2
% 
% r1 = [1 2 3 4 5 6]'; 
% p = [2 4 5 8 9 7]'; 
% r2 = [10 15 1 0 7 9]'; 
% 
% func = @(r1,p,r2)( 5 .*r1 .^2 + 2 .*p .^2 + r2 .^2 );
% % [r1,p,r2] = meshgrid(-1:0.1:1, -1:0.1:1, -1:0.1:1);
% v = func(r1,p,r2);
% 
% scatter3(r1, p, r2, 40, v, 'filled')
% 
% for x=0:length(p)
%     p=circshift(p,1);
%     for w=0:length(r2)
%         r2=circshift(r2,1);
%         for u =1:1:length(r1)
%             v(u)=5*r1(u).^2 + 2*p(u).^2 + r2(u).^2;
%         end
%         
%         %scatter3(X, Y, Z, circleSize, circleColor)
%         scatter3(r1, p, r2, 40, v, 'filled')
%         xlabel('r1') 
%         ylabel('p')
%         zlabel('r2')
%         hold on
%     end
% end
% 
% cb = colorbar;
  
%% Method 3

func = @(x,y,z)(sin(2*x.^2+2*y.^2).*cos(2*z.^2));
[x,y,z] = meshgrid(-1:0.1:1, -1:0.1:1, -1:0.1:1);
v = func(x,y,z);

