%=======================================================================
% Exercise 8.39
% The script  generates the data and the network. 
% The  CTA APSM is compared to
% the CTA LMS, the ATC LMS and the non-cooperative LMS algorithms.
%======================================================================

function distributed_test_APSM

 
close all
err = {};
leg = {};
N=2000;

resultsfolder = '/results_timevar/';
%rng(123,'twister')
plots = 0;

%% Make unknown vector
L=60;
 
rng(12,'twister')
 

nodes = 10;
connections = 32;

 

% for node = 1:nodes
% rndstreams{node}=RandStream('mrg32k3a','Seed',node);
% end
%inputvec = @(node,t) inputvector(rndstreams,node,t,L);


%% Make the network
distributed =1;
if distributed
    ruleweights = 'Metropolis';
    %ruleweights = 'noncoperation';
    [Nmat,C,connected] = makenetwork(nodes,connections,ruleweights,100,0);
    if ~connected
        disp('I did not manage to build a connected network')
        return
    end
else
    C=eye(nodes);
    connections=0;
end

totalrep=100;
for inod=1:nodes, errlms{inod} = zeros(1,N); end
totalerr1=zeros(1,N);
totalerr2=zeros(1,N);
totalerr3=zeros(1,N);
totalerr4=zeros(1,N);
ensaver=1;
fullmatinput=0; if fullmatinput==1, X=cell(1,nodes); end
namelms = ['LMS','L',num2str(L),  '_Nd',num2str(nodes),'_Cn',num2str(connections)];
%% Start ensemble averages
for rep=1:totalrep
    namelms = ['LMS','_r',num2str(rep),'_L',num2str(L), '_Nd',num2str(nodes),'_Cn',num2str(connections)];

rng(rep,'twister')
 h1=randn(L,1);  
h2=h1;
changepoint = round(N/2);
h={h1,h2};
h =  @(i) h{(i>changepoint)+1};
%h = {h1,h2,changepoint};
        for j=1:nodes
        X{j}=randn(L,N);    
        end
    if fullmatinput==0
   inputvec = @(node,k) X{node}(:,k);
    else
        inputvec = @(node,i) X{node}(:,i);
    end


%% Tune the noise variance in order to have a certain SNR in the first section (which has \norm{x}^2 input energy) 
%sigmaxsqr = 1/12*(1-0)^2; %variance of uniform distribution
sigmaxsqr = 100^2/L; %The input vectors have normalized norm equal to 100
SNR_high = 25;
SNR_low =20;
SNR=linspace(SNR_low,SNR_high,nodes);
noisevec=cell(nodes,1);
for inode = 1:nodes
    noisevar(inode) = (var(h1)*sigmaxsqr)/10^(SNR(inode)/10);
    st = RandStream('mlfg6331_64','Seed',nodes+inode);
    noisevec{inode} = single(randn(st,N,1)*sqrt(noisevar(inode)));
end


normh1=norm(h1); 
normh2=normh1;
 % Determine how to compute the error
errfun = @(i,h_hat) comperr(h_hat,i,changepoint,h1,normh1,h2,normh2);

  
   
    for i=1:nodes
       if ~isa(h,'function_handle')
            y_clean{i} = X{i}'*h;
        else
            y_clean{i} =  X{i}'*h(i);
        end
        y_noise{i} = y_clean{i} + noisevec{inod}(i);
      
    end
        
playLMS=1;
if playLMS==1
 
 mu = 0.01;

LMS_param.mu = 0.03;
LMS_param.normalized = 0.0;
LMS_param.gamma = 0; %0.25*10^(-3);
LMS_param.L = L; LMS_param.N = N;
LMS_param.h = h; LMS_param.C = C;
LMS_param.nodes = nodes;
LMS_param.dataexchange = 0; %if it excanges data with its neighbors

APSM_param.mu =1;APSM_param.L = L;
 APSM_param.q = 20;  
APSM_param.epsilon = 0.5; APSM_param.N = N;
APSM_param.h = h; APSM_param.C = C;
APSM_param.nodes = nodes;
 

if isa(mu,'function_handle')
    mutag = ['_mudmn',num2str(g0),'-',num2str(mu_n),'-',num2str(c)];
else
    mutag = ['_mu',num2str(mu)];
end

end
[err1, x]=LMS_distrib_ATC(inputvec,y_noise,LMS_param,errfun)
[err2, x]=LMS_distrib_CTA(inputvec,y_noise,LMS_param,errfun)
errmeanlms1=zeros(1,N);
errmeanlms2=zeros(1,N);
errmeanlms3=zeros(1,N);
errmeanlms4=zeros(1,N);
LMS_param.C=eye(nodes); %%%Here is the nondistributed case
[err3, x]=LMS_distrib_ATC(inputvec,y_noise,LMS_param,errfun)


errmeanapsm=zeros(1,N);
[err4, x]=APSM_distrib(inputvec,y_noise,APSM_param,errfun)

for inod=1:nodes errmeanlms1=errmeanlms1+err1{inod}; end
for inod=1:nodes errmeanlms2=errmeanlms2+err2{inod}; end
for inod=1:nodes errmeanlms3=errmeanlms3+err3{inod}; end
for inod=1:nodes errmeanlms4=errmeanlms4+err4{inod}; end

errmeanlms1 = errmeanlms1/nodes;
errmeanlms2 = errmeanlms2/nodes;
errmeanlms3 = errmeanlms3/nodes;
errmeanlms4 = errmeanlms4/nodes;


totalerr1=totalerr1+errmeanlms1;
totalerr2=totalerr2+errmeanlms2;
totalerr3=totalerr3+errmeanlms3;
totalerr4=totalerr4+errmeanlms4;
end
% if playLMS==1
 
  
 totalerr1 = totalerr1/totalrep;
  totalerr2 = totalerr2/totalrep;
  totalerr3 = totalerr3/totalrep;
 totalerr4 = totalerr4/totalrep;
  %the average of all the nodes.
%end
 plot(totalerr1,'r');hold on
  plot(totalerr2,'g');hold on
  plot(totalerr3,'k');hold on 
   plot(totalerr4,'m') 
end


function vec = inputvector(rndstreams,stream,substream,L)
sstream = rndstreams{stream};
RandStream.setGlobalStream(sstream);
sstream.Substream=substream;
vec = randn(L,1);
end

 function vec = inputvector2(stream,node,t,L)
stream.Substream=node*(t-1)+node;
vec = rand(L,1);
 end
   
 function vec = inputvector3(node,t,L)
 sstream = RandStream('swb2712','Seed',node*(t-1)+node);
 set(sstream,'FullPrecision',false)
 vec = rand(sstream,L,1)-0.5;
 end

function vec = inputvector4(rndstreams,node,L)
vec = rand(rndstreams{node},L,1,'single')-0.5;
vec = vec/norm(vec)*100; % So the variance of vec equals to 100^2/L
% vec = randn(rndstreams{node},L,1)*sqrt(0.1);
end

function vec = inputvector5(X,node,i)
vec = X{node}(:,i);
end

function err = comperr(h_hat,i,change_i,h1,normh1,h2,normh2)

if i<=change_i
    % err = 20*log10(norm(h1-h_hat)/normh1);
     err = 20*log10(norm(h1-h_hat));
else
    % err = 20*log10(norm(h2-h_hat)/normh2);
    err = 20*log10(norm(h2-h_hat));
end

end