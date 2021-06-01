function [c] = zerocell(dims)
%     if dims > 3
% %         error('This function does not support >3 dimensions yet.')
%     end
    
    c = cell(dims,dims);
    for i=1:dims
        if dims > 1
            for j=1:dims
                c{i,j} = 0;
            end
        else
            c{i} = 0;
        end
    end
    
end