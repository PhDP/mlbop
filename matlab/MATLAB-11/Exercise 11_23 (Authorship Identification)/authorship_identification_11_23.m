%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 11.23
% Authorship identification.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ** IMPORTANT **: You must have installed the Statistics toolkit from the
% main Matlab installation. 
% You must also have installed the TMG toolkit 
% (available at http://scgroup20.ceid.upatras.gr:8000/tmg/ at the time
% of writing)
% This script also uses LibSVM for Matlab
% ** IMPORTANT ** : Prior to the use, run "make" to compile LibSVM for Matlab.

% The data set is available via the UCI Machine Learning Repository, http://archive.ics.uci.edu/ml
% See, also, description in the text.

% IMPORTANT: Set the following variable to the full path of the directory 
% containing this m file. E.g. if the file is '/full/path/to/authorship_identification_11_23.m'
% then you should write '/full/path/to/'

baseDir = '/full/path/to/'; % Update this!!!
cd(baseDir);

% N-gram graphs (NGGs)

% Initialize libraries for N-gram graphs
% Built NGGSVMKernel.jar from:
% https://github.com/ggianna/ngg2svm/
% and downloaded (JInsect.jar, OpenJGraph.jar) from
% http://sourceforge.net/projects/jinsect/
% at the time of writing.
javaaddpath JInsect.jar;
javaaddpath OpenJGraph.jar;
javaaddpath NGGSVMKernel.jar;

% Extract similarity matrices from training data
% based on NGG similarity (calculated through the Java libraries)
a=javaArray('java.lang.String', 2);
a(1)=java.lang.String(strcat(baseDir, 'data/train/'));
a(2)=java.lang.String(strcat(baseDir, 'nggkernel.txt'));
javaMethod('main','nggsvmkernel.Main',  a); % This may take some time (e.g. ~3-5 mins)

% Add path to the libsvm toolbox
addpath(strcat(baseDir, 'libsvm-3.20/matlab'));

% load precomputed NGGs kernel 
[nggLabels,nggData] = libsvmread(strcat(baseDir, 'nggkernel.txt'));
% perform cross-validation and return results
libsvmcv(nggLabels, nggData, '-t 4') % Once again you may need to wait for
                                     % some time...


% Bag of words

% ** WARNING ** : You must have installed the TMG toolkit
opt.delimiter = 'none_delimiter';
% Load data from directory
% NOTE: You may need to confirm loading all directories at this point
% by typing in "yes" when asked.
docMat = tmg(strcat(baseDir, 'data/train/'), opt);
% return to base directory, because TMG does not return
cd(baseDir);
% Initialize class labels vector
aClasses=[1*ones(50,1); 2*ones(50,1); 3*ones(50,1);4*ones(50,1)];
% run cross-validation on bag-of-words data on RBF kernel
libsvmcv(aClasses, docMat)


% NOTE: Data taken from https://archive.ics.uci.edu/ml/datasets/Reuter_50_50
% as indicated in the book.