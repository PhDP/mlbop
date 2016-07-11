function   vec = sum_mean_val2(input, K)


N=length(input);
for i=1:N
    vec(i)=mean(input(max(1,i-K+1):i));
end;