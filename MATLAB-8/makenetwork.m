%=====================================================================
% Exercise 8.39
% This script generates the network for the distributed adaptive
% learning experiments.
%=====================================================================

function [N,A,connected] = makenetwork(nodes,connections,ruleweights,trials,visualize)

% Builds a network with ``nodes'' nodes and ``connections'' connections.
% ruleweights: 'Metropolis', 'uniform', 'noncoperation'
% tryconnectivity: Runs many realizations in order to get a connected
% network. However, it can still fail to do so.
% visualize = 1 to plot the graph.
%Returns:
% matrix N, with ones and zeros in connected nonconnected edges
% matrix C, with metropolis rule entries
% logical connected: True for connected -- False for non connected

if nargin==4
    visualize = 0;
elseif nargin==3
    visualize = 0;
    trials = 1;
end

Lowtriangular = tril(ones(nodes,nodes),-1);
Lowtriangpos=find(Lowtriangular);
if connections>length(Lowtriangpos), connections=length(Lowtriangpos); end

connected = false;
tr = 1;
while ~connected && tr<=trials
    if mod(tr,1000)==0, disp(['this is the ',num2str(tr),' try to find a connected']), drawnow, end
        
    N = zeros(nodes,nodes);
    
    K = randperm(length(Lowtriangpos),connections);
    N(Lowtriangpos(K)) = 1;
    N=N+N'+eye(nodes); % make it symmetric with full diagonal
    N=N~=0;
    [p,q,r,s] = dmperm(N);
    if length(r)==2,
        connected = true;
    end
    tr = tr+1;
end

A=zeros(size(N));
N_k = sum(N);

if strcmpi(ruleweights,'metropolis')
    for jj=1:nodes
       for ii=1:nodes
          if N(ii,jj)~=0 && ii~=jj
              A(ii,jj)=1/max([N_k(ii),N_k(jj)]);
          end    
       end
       A(jj,jj) = 1 - (sum(A(:,jj))-A(jj,jj));
    end
elseif strcmpi(ruleweights,'uniform')
    for jj=1:nodes
        A(:,jj) = N(:,jj)/N_k(jj);
    end
elseif strcmpi(ruleweights,'noncoperation')
    A=eye(nodes);
end

if visualize
xy = [randperm(nodes);randperm(nodes)]';
figure;gplot(N,xy)
end