function [rowindex, colindex] = findMaxRowAndCol(image)

    n = size(image,1);
    m = size(image,2);
    
    total = 0;
    index = 0;
    
    %{
    for i = 1:n
        previous = total;
        iprev = index;
        total = sum(image(i,:));
        index = i;
        if total < previous
            total = previous;
            index = iprev;
        end
    end
    %}
    if index == 0
        im1 = sum(image,1);
        [val,col] = max(im1);
        index = col;
    end
    
    colindex = index;
    
    total = 0;
    index = 0;
    
    %{
    for j = 1:m
        previous = total;
        iprev = index;
        total = sum(image(:,j));
        index = j;
        if total < previous
            total = previous;
            index = iprev;
        end
    end
    %}
    
    if index == 0
        if index == 0
            im1 = sum(image,2);
            [val,row] = max(im1);
            index = row;
        end
    end
    
    rowindex = index;
end