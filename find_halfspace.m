function l = find_halfspace(data,x);
% returns the number of the last sample in the ASCENDING array data that is
% &lt; x using a simple half space search
 
    upper = numel(data); lower = 1;
    [data,x,l,u] = find_halfspace_rec(data,x,lower,upper);
end

function [data,x,lower,upper]= find_halfspace_rec(data,x,lower,upper);
    mp=floor((lower+upper)/2);
    if mp~=lower
        if data(mp)<x
            [data,x,lower,upper]= find_halfspace_rec(data,x,mp,upper);
        else
            [data,x,lower,upper]= find_halfspace_rec(data,x,lower,mp);
        end
    end
end

