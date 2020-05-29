function Dispersion3D = Dispersion_SpinOrbit3d(Kappa)
Dispersion3D = cell(3);
Dispersion3D{1,1} = @(fftx, ffty, fftz) 0.5 * (fftx.^2 + ffty.^2 + fftz.^2);
Dispersion3D{2,2} = @(fftx, ffty, fftz) 0.5 * (fftx.^2 + ffty.^2 + fftz.^2);
Dispersion3D{3,3} = @(fftx, ffty, fftz) 0.5 * (fftx.^2 + ffty.^2 + fftz.^2);

Dispersion3D{1,2} = @(fftx, ffty, fftz) Kappa * (fftx - 1i * ffty);
Dispersion3D{2,1} = @(fftx, ffty, fftz) Kappa * (fftx + 1i * ffty);

Dispersion3D{1,3} = @(fftx, ffty, fftz) 0;
Dispersion3D{3,1} = @(fftx, ffty, fftz) 0;
Dispersion3D{2,3} = @(fftx, ffty, fftz) 0;
Dispersion3D{3,2} = @(fftx, ffty, fftz) 0;

%Dispersion3D{1,3} = @(fftx, ffty, fftz) Kappa * (fftx - 1i * ffty);
%Dispersion3D{3,1} = @(fftx, ffty, fftz) Kappa * (fftx + 1i * ffty);
%Dispersion3D{2,3} = @(fftx, ffty, fftz) Kappa * (fftx - 1i * ffty);
%Dispersion3D{3,2} = @(fftx, ffty, fftz) Kappa * (fftx + 1i * ffty);

end 