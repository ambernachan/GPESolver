function Dispersion2D = Dispersion_SpinOrbit2d(Kappa)
Dispersion2D = cell(2);
Dispersion2D{1,1} = @(fftx, ffty) 0.5 * (fftx^2 + ffty.^2);
Dispersion2D{1,2} = @(fftx, ffty) Kappa * (fftx - 1i * ffty);
Dispersion2D{2,1} = @(fftx, ffty) Kappa * (fftx + 1i * ffty);
Dispersion2D{2,2} = @(fftx, ffty) 0.5 * (fftx^2 + ffty.^2);
end 