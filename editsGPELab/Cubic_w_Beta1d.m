function [CubicNonlinearity] = Cubic_w_Beta1d(Method, Beta)

CubicNonlinearity = cell(Method.Ncomponents);
for n = 1:Method.Ncomponents
        CubicNonlinearity{n,n} =  @(Phi,X) Beta * abs(Phi{n}).^2;
end

end
