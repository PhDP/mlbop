%-----------------------------------------------------------------
%  MATLAB code for Exercise 7.24
%  Protein folding prediction
%  Requires: Statistics Toolbox (Machine Learning >> CART)
%-----------------------------------------------------------------

fprintf('DEMO: CART decision trees for protein folding prediction\nDataset: UCSD-MKL / 27-class subset of the PDB-40D SCOP\n\n');

clear all;
load 'ppdemo_data.mat';			% imported dataset from PDB-40D SCOP

%...  CART: classification mode  ...
% Using: 'ppdata_trn_comp' = feature vectors (composition)
         'ppdata_trn_id'   = assigned class (1...27)
% See: 'help fit' for details on training parameters  
fprintf('Creating classification tree...\n');
treeC=ClassificationTree.fit(ppdata_trn_comp,ppdata_trn_id,'prior','empirical','ScoreTransform','none','SplitCriterion','gdi','ScoreTransform','symmetric','MinParent',5);
disp(treeC);

% evaluate performance with the testing subset, show results (total success rate)
evalC=treeC.predict(ppdata_tst_comp);
accC=length(find(evalC==ppdata_tst_id))/length(ppdata_tst_id);		% sum correct predictions over all classes
fprintf('Final classification accuracy (1-vs-all):  %g%%\n\n',accC*100);

% display the CART classifier (full tree plot)
view(treeC,'mode','graph');
