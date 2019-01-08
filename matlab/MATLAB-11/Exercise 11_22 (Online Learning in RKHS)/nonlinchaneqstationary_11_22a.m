%-----------------------------------------------------------------
% Exercise 11.22
% Nonlinear Channel Equalization in RKHS
% Part A: Stationary Case
%-----------------------------------------------------------------


clear;
rng(0);

N=5000;
%input signal s
%s = abs(s)>=0.5;

%noise n



NUMBER_OF_TESTS = 1000;

TOTAL_KNLMS = 0;
TOTAL_APSM = 0;
TOTAL_NORMA = 0;
    
TOTAL_EXPANSION_KLMS = 0;
TOTAL_EXPANSION_KAPSM = 0;
TOTAL_EXPANSION_NORMA = 0;


p=5;
h=[-0.9   0.6  -0.7  0.2   0.1]';
residue(1,h)

for t=1:NUMBER_OF_TESTS

    var_signal=0.8;
    s = var_signal*randn(1,N);

    
    cur_snr=15;
    var_noise = sqrt( var_signal^2/(10^(cur_snr/10)) );


    %gaussian noise n
    noise = var_noise*(randn(1,N));
    
    
    SNR_num=10*log10(var(s)/var(noise));
    SNR=10*log10(var_signal^2/var_noise^2);
    %disp( [''] );
    %disp( [''] );
    disp('-------------------------------------------------------------------');
    disp(['TEST NUMBER : ', num2str(t)]);


    
    %test 1 linear channel and then non linear channel without memory
    for n=1:p
        t(n)=s(n);
    end;
    for n=p+1:N
        xi_n = (s(n:-1:n-p+1))';
        t(n)= transpose(h)*xi_n + noise(n);
        x(n)= t(n) + (0.15)*t(n)^2 + 0.03*t(n)^3  + noise(n);
        r(n)=x(n);
    end;
    

    
    

    L=5;
    D=2;

    %construct z and d for the LMS
    N0=-D+L;
    N1=N-L+1;
    for n=-D+L:N-D
        %z(n)=[r(n+D-L+1) r(n+D-L+2)     r(n+D)]
        z(n-N0+1,:)=r(n+D-L+1:n+D);
        d(n-N0+1)=s(n);
        r2(n)=r(n+D);
    end;

    %disp('*******************************');

    mu=1/2;
    sigma = 5;
    Q_size=5;
    sparse_params = [Q_size];
    %tic;
    %disp('NKLMS started');
    [a1, centers1, e_KNLMS, expansion_size_klms]=QKernel_NLMS(z,d,N1,mu,1,0,'gaus_f',[sigma], sparse_params);
    %disp('NKLMS ended');
    %toc;
    
    
    mu=1/4;
    lambda=0.01;
    sigma = 5;
    sparse_params = [80]; %size of window
    %tic;
    %disp('NORMA L2 started');
    [a2, centers2, e_NORMA, expansion_size_norma]=NORMA_L2(z,d,N1,mu,lambda,0,'gaus_f',[sigma], sparse_params);
    %disp('NORMA L2 ended');
    %toc;
    
    
    
    %APSM
    %tic;
    %disp('APSM started');
    epsilon = 10^(-5);
    Delta=10000;
    thresh1 = 0.1;
    thresh2 = 0.1;
    sigma=5;
    huber_sigma = 2;
    Q=5;
    [a4,centers4, e_APSM,expansion_size_KAPSM] = QKernel_APSM(z,d,N1,epsilon,Q,'l2',[huber_sigma],'gaus_f',[sigma],[Delta, Q_size],0);
    expansion_size_APSM=0;
    %[a1,centers,e1] = Kernel_APSM(z,d,N1,epsilon,Q,'l2',[huber_sigma],'gaus_f',[sigma],2,[thresh1 thresh2],1);
    %disp('APSM ended');
    %toc;


    
   
    
    
    %disp('*******************************');

    
    
    
    
    TOTAL_KNLMS = TOTAL_KNLMS + abs(e_KNLMS).^2;
    TOTAL_NORMA = TOTAL_NORMA + abs(e_NORMA).^2;
    TOTAL_APSM = TOTAL_APSM + abs(e_APSM).^2;
    
    TOTAL_EXPANSION_KLMS = TOTAL_EXPANSION_KLMS + expansion_size_klms;
    TOTAL_EXPANSION_NORMA = TOTAL_EXPANSION_NORMA + expansion_size_norma;
    TOTAL_EXPANSION_KAPSM = TOTAL_EXPANSION_KAPSM + expansion_size_KAPSM;
end;

MEAN_KNLMS = TOTAL_KNLMS / NUMBER_OF_TESTS;
MEAN_NORMA = TOTAL_NORMA / NUMBER_OF_TESTS;
MEAN_APSM = TOTAL_APSM / NUMBER_OF_TESTS;



MEAN_EXPANSION_KLMS = TOTAL_EXPANSION_KLMS / NUMBER_OF_TESTS;
MEAN_EXPANSION_NORMA = TOTAL_EXPANSION_NORMA / NUMBER_OF_TESTS;
MEAN_EXPANSION_KAPSM = TOTAL_EXPANSION_KAPSM / NUMBER_OF_TESTS;


figure(1);
hold on;
ylabel('MSE');
xlabel('n');
box    = [1, N, -8, 0];
axis(box);



%print 1
K=10;
plot(10*log10(sum_mean_val2(MEAN_NORMA(1:N1),K)),'k','LineWidth',1);
plot(10*log10(sum_mean_val2(MEAN_KNLMS(1:N1),K)),'r','LineWidth',1);
plot(10*log10(sum_mean_val2(MEAN_APSM(1:N1),K)),'Color', 0.5*[1 1 1],'LineWidth',1);


axes_handle=get(gcf,'CurrentAxes');
set(axes_handle,'YGrid','on');

h = legend('NORMA', 'QKNLMS', 'QAPSM' ,  'Location', 'NorthEast');
set(h,'Interpreter','none');
title(['Non linear channel Equalization']);






figure(2);
hold on;
ylabel('Expansion size (M)');
xlabel('n');

%line([1,N], [1, N], 'Color', 'k', 'LineWidth',2);
plot(MEAN_EXPANSION_NORMA,'k', 'LineWidth', 2);
plot(MEAN_EXPANSION_KLMS,'b', 'LineWidth', 2);
plot(MEAN_EXPANSION_KAPSM,'g', 'LineWidth', 2);

axes_handle=get(gcf,'CurrentAxes');
set(axes_handle,'YGrid','on');
h = legend( 'NORMA L_2', 'QKNLMS',  'QAPSM', 'Location', 'NorthEast');
title(['Evolution of the Expansion`s size']);
    
    
       