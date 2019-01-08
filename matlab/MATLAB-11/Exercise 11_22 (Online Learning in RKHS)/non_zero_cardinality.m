

function r = non_zero_cardinality(vector)

r=0;
N = length(vector);
for i=1:N
    if vector(i)~=0
        r = r+1;
    end;
end;