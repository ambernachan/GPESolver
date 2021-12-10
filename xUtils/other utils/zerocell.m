function [c] = zerocell(dims)
%     if dims > 3
% %         error('This function does not support >3 dimensions yet.')
%     end
    
    c = cell(dims,dims);
    for i=1:dims
        if dims > 1
            for j=1:dims
                if dims > 2
                    c{i,j} = @(Phi, X, Y, Z) 0;
                else
                    c{i,j} = @(Phi, X, Y) 0;
                end
            end
        else
            c{i} = @(Phi, X) 0;
        end
    end
    
end