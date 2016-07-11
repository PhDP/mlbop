%========================================================
% Case Study 14.11: Change Point Detection
% The code has been kindly provided by Taylan Cemgil
%========================================================

% CMPE58N_CP_POISS        Gibbs sampler for the Coal Mining Data

% Change History :
% Date        Time        Prog    Note
% 24-Mar-2009     7:53 PM    ATC    Created under MATLAB 7.7.0

% ATC = Ali Taylan Cemgil,
% Department of Computer Engineering, Bogazici University
% e-mail : taylan.cemgil@boun.edu.tr

%randn('seed', 2);
rand('seed', 1);

%Coal mining data
x=[4 5 4 1 0 4 3 4 0 6 3 3 4 0 2 6 3 3 5 4 5 3 1 4 4 1 5 5 3 4 2 5 2 2 ...
3 4 2 1 3 2 2 1 1 1 1 3 0 0 1 0 1 1 0 0 3 1 0 3 2 2 0 1 1 1 0 1 0 1 0 ...
0 0 2 1 0 0 0 1 1 0 2 3 3 1 1 2 1 1 1 1 2 4 2 0 0 0 1 4 0 0 0 1 0 0 0 ...
0 0 1 0 0 1 0 1];
M = length(x);

mm = ceil(M*rand);

a = 2; b = 1;

figure(1)
clf
EP = 1000;

BURN_IN = 200;
verbose = 0;

m = 10;
lam = zeros(2,1);

subplot(211);
set(gca, 'xlim', [0 EP])

subplot(212);
set(gca, 'xlim', [0 EP])

gibbs.m = zeros(1, EP-BURN_IN);
gibbs.lambda = zeros(2, EP-BURN_IN);


for e=1:EP,
  lam(1) = gamrnd(a + sum(x(1:m)), 1/(m+b) );
  lam(2) = gamrnd(a + sum(x(m+1:end)), 1/(M-m+b) );


  lp = zeros(1, M);
  for i=1:M,
    lp(i) = sum(x(1:i)).*log(lam(1)) - i*lam(1) + ...
sum(x((i+1):M)).*log(lam(2)) - (M-i)*lam(2) ;
  end;

  m = randgen(exp(lp-max(lp)));

  if verbose,
    subplot(211); line([e e], [lam(1) lam(2)], 'marker', '.', 'lines','n' );
    subplot(212); line([e], [m],  'marker', 'o'  );
    drawnow;
  end;

  if e>BURN_IN,
    gibbs.lambda(:, e-BURN_IN) = lam;
    gibbs.m(1, e-BURN_IN) = m;
  end;


end;

figure(1)
clf
plot(gibbs.lambda(1,:),gibbs.lambda(2,:),'.');
set(gca, 'xlim', [0 4], 'ylim', [0 4]);
axis square; xlabel('l1'); ylabel('l2');
%subplot(221);
%hist(gibbs.lambda(1,:),100);

figure(2)
subplot(211)
idx = 1850 + (1:length(x));
stem(idx, x);
ylabel('xj')
set(gca, 'xlim', [min(idx) max(idx)+1])


subplot(212);
hs = accumarray(gibbs.m', 1, [length(x) 1]);
bar(idx, hs/sum(hs))
%set(gca, 'xlim', [0 M+1])
set(gca, 'xlim', [min(idx) max(idx)+1])
ylabel('p(tau|x)')
xlabel('j');