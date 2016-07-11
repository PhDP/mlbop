%-----------------------------------------------------------------
%  Exercise 5.23
%  Decision Feedback Equalization 
%-----------------------------------------------------------------





rseed	= 1;
rand('state', rseed);
N_T = 500;

N=10000; %The length of the data sequence
L = 31;
Lb = 10;
Lf = 21;

%the impulse responses of the linear channel
h = [0.04 -0.05 0.07 -0.21 0.72 0.36 0.21 0.03 0.07]';

Lh = length(h);
mse = zeros(N,1);
u = zeros(N,1);

N0 = 250; %After N0 samples we assume convergence
errors = 0;

%Run N_T experiments
for i = 1:N_T
    
    %Generate a set of 10000 random ±1 values (BPSK)
    s = sign(2*rand(N,1)-1);
    
    %Direct this sequence into a linear channel
    for n=1:N
        u(n) = s(n:-1:max(1,n-Lh+1))'*h(1:min(Lh,n));
    end;
    
    %add noise
    u=awgn(u,11,'measured');


    
    % Adaptive Decision Feedback Equalizer 
    % Estimate the original signal (s) 
    % from u (the output of the linear channel)
    w = zeros(Lf+Lb,1);
    e = zeros(N,1);
    d = zeros(N,1);
    mu = 0.025;
    %filter
    for n=Lf+Lb:N
        if n<=N0
            v_n = [u(n:-1:n-Lf+1); s(n-Lf:-1:n-Lf-Lb+1)];
            d_hat = w'*v_n;
            d(n) = s(n-Lf+1);
        else
            v_n = [u(n:-1:n-Lf+1); d(n-1:-1:n-Lb)];
            d_hat = w'*v_n;
            d(n) = sign(d_hat);
            if d(n) ~= s(n-Lf+1)
                errors = errors + 1;
            end;
        end;
        e(n) = d(n) - d_hat;
        w = w + mu*v_n*e(n);
        e(n) = d_hat - s(n-Lf+1);
    end;
    % end DFE
    
    %compute the mean square error over all experiments
    mse = mse+e.^2;
end;

%plot results
mse = mse/N_T;
plot(10*log10(mse),'r');
xlabel('n');
ylabel('MSE');
percentage = errors/(N-N0)/N_T;
disp(['Errors: ', num2str(percentage*100,4),'%']);
    