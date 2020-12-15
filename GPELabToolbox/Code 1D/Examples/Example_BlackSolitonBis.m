%%% This file is an example of how to use GPELab (FFT version)

%% Ground state of a Gross-Pitaevskii equation with quadratic potential and cubic nonlinearity in 1D


clear all;

%% Setting the method and geometry
Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Relaxation';
Deltat = 1e-3;
Stop_time = 1;
Stop_crit = [];
Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -20;
xmax = 20;
Nx = 2^11+1;
Geometry1D = Geometry1D_Var1d(xmin,xmax,Nx);

%% Setting the initial data
Beta = 1;

L_h = sqrt(2);
c = sqrt(Beta/2);
v = 0.1*c;
Xi = sqrt(1-(v/c)^2);
X_0= 2;
X = Geometry1D.X ; 
eps = 0.4;
Phi_stock = zeros(Nx-2,floor(Stop_time/Deltat)/10);
Phi = Xi*tanh(Xi*X(2:Nx-1)/L_h)+1i*(v/c);

%% Setting informations and outputs
Lpl_vect = (1/Geometry1D.dx^2)*ones(Nx-2,1);
Lpl_mat = spdiags([-Lpl_vect,2*Lpl_vect,-Lpl_vect],[-1,0,1],Nx-2,Nx-2);
Lpl_mat(1,1) = (1/2)*Lpl_mat(1,1);
Lpl_mat(Nx-2,Nx-2) = (1/2)*Lpl_mat(Nx-2,Nx-2);
Id_mat = speye(Nx-2,Nx-2);
NL_mat = spdiags(Beta*abs(Phi).^2,0,Nx-2,Nx-2);

for k = 1:floor(Stop_time/Deltat)
    NL_mat = spdiags(Beta*abs(Phi).^2,0,Nx-2,Nx-2);
    Phi = (Id_mat+Method.Deltat*1i*(Lpl_mat+NL_mat))\Phi;
    plot(X(2:Nx-1),abs(Phi).^2)
    drawnow
    if (mod(k,10) == 0)
    Phi_stock(:,k/10) = Phi; 
    end
end

Time = [10*Deltat:10*Deltat:Stop_time];

surf(Time,X(2:Nx-1),abs(Phi_stock),'EdgeColor','none')
view(2)
ylabel('x')
xlabel('Time')

