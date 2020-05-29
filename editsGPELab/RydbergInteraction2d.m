function [Dipolar2d] = RydbergInteraction2d(R0, phi, fftx, ffty)
q = R0*sqrt(fftx.^2+ffty.^2);
const = R0*(2*pi^2/3);
K = const*(exp(-q/2)./q).*(exp(-q/2)-2*sin(pi/6-sqrt(3)/2*q));
K(isinf(K)) = const;
Dipolar2d = (1/2)*real(ifft2(K.*fft2(abs(phi).^2)))';
end