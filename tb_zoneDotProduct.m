function [ x ] = tb_zoneDotProduct( z , performAverageFirst )
%tb_zoneDotProduct Generates dot product of a group of zones,
% perform average first generates an average of the signal.  Without it,
% the code is cascading the dot products.
     [~, o, n] = size(z);
     x = zeros(o);
    if (performAverageFirst == 0)
        for j=1:1:n
            data = z(:,:,j);
            data( ~any(data,2), : ) = [];  %rows
            [m, ~] = size(data);
            for k=1:1:m
                if(k==1 && j==1)
                    x = data(k,:)/max(data(k,:));
                else
                    x = x.*(data(k,:)/max(data(k,:)));
                end
            end
        end
    else
        count = 0;
        for j=1:1:n
            data = z(:,:,j);
            data( ~any(data,2), : ) = [];  %rows
            [m, ~] = size(data);
            for k=1:1:m
                count = count + 1;
                if(k==1 && j==1)
                    x = data(k,:)/max(data(k,:));
                else
                    x = x + ( data(k,:)/max(data(k,:)) );
                end
            end
        end  
        x = x / count;
    end

end

