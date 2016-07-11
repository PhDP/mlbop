%----------------------------------------------------------------
% Exercise 10.14
% Time-frequency analysis of echolocation signals transmitted by bats.
% The Time-Frequency toolbox (http://tftb.nongnu.org/) need to be
% downloaded and added to the matlab path
% addpath('tftb/mfiles').
%----------------------------------------------------------------

function Matlab10_14


        
        addpath 'C:\Users\Usuario\Google Drive\documents\DSP_mathima\tftb\mfiles'
        
        % The NESTA toolbox need to be downloaded and added to the matlab path (http://statweb.stanford.edu/~candes/nesta/nesta.html).
        %Hadamard function (distributed with NESTA) might also need compiling
        addpath(genpath('NESTA_v1.1'))
        
        % The Large Time-Frequency Analysis Toolbox
        % (http://ltfat.sourceforge.net/download.php) need to be downloaded
        addpath('ltfat'); ltfatstart
        
        [x_exact,FS,Nbits]=wavread([pwd,filesep,'mixed_calls']);
        N = length(x_exact);
        topfreq=0.15; Nw=2^10+1; Nf=4*1024; Ntpl=4*512;
        
        
        % The Hadamard code used below requires inputs of 2^k
        % So, padd the signal with zeros at the end
        k = nextpow2(N); x_exact = [x_exact; zeros(2^k-N,1) ]; N_short = N; N = 2^k; 
        M = round(N/8);  % number of measurements

        seed = 1234;     % make it reproducible
        randn('state',seed); rand('state',seed);

        % The following formulation is taken from NESTA examples. See NESTA
        % documentation for similar examples.
        
        % --- Randomly permute the columns ---
        % When x is a column vector it essentially permutes its components.
        % However, when the Hadamard transform of the permuted x is taken (see
        % operator A below) it is equivalent to have the columns of the Hadamard
        % matrix permuted. In this way we can use measurement vectors with +/- ones
        % and on the same time, we take the measurements with a fast Hadamard transform.
        perm = randperm(N);      % pick a permuation
        permute_cols = @(x) x(perm,:);
        Sperm.type = '()'; Sperm.subs{1} = perm; Sperm.subs{2} = ':';
        i_permute_cols = @(x) subsasgn( zeros(size(x)),Sperm,x);
        
        % --- Randomly subsample the rows
        ROWS = randperm(N); ROWS = ROWS(1:M);   % pick some rows
        downsample = @(x) x(ROWS,:);
        S.type = '()'; S.subs{1} = ROWS; S.subs{2} = ':';
        upsample = @(x) subsasgn( zeros(N,size(x,2)),S,x);
        
        % -- And make the operator
        sqrtN = 1/sqrt(N);
        A = @(x) sqrtN*downsample( hadamard( permute_cols(x) ) );
        At = @(x) sqrtN*i_permute_cols( hadamard( upsample(x) ) );
        
       
      
       
%% Asuming the singal is sparse in Gabor - moderate redundancy        
 % Construction of GAbor window and its dual       
 Time_samples=2^7;
 Freq_shift=2^6;
 gabor.a = N / Time_samples;  %  time shift
 gabor.M = N / Freq_shift; % frequency bands
 
 % For this Gabor system, the optimally concentrated Gaussian is given by
 tf_ratio = gabor.a/Freq_shift;
 
 % Compute the window and the canonical dual window
 gabor.w = pgauss(N,tf_ratio*4);
 gabor.wd=gabdual(gabor.w,gabor.a,gabor.M);
 
           
yy = dgt(ones(size(x_exact)),gabor.w,gabor.a,gabor.M); 
gabor.fpoints=size(yy,1);
gabor.tpoints=size(yy,2);
clear yy
psi_A = @(x) op_igabor(x,gabor,N); 
psiT_A = @(y) op_gabor(y,gabor);            

c = psiT_A(x_exact);
cc=sort(abs(c).^2,'descend');
tot_en = sum(cc);
tot99=tot_en*99/100;
cc_cum = cumsum(cc);
Sparsity = length(cc)-length(find(cc_cum>tot99));
Sparsityratio = Sparsity/length(cc);


c_real=dgtreal(x_exact,gabor.w,gabor.a,gabor.M);
gabor.fpoints=size(c_real,1); gabor.tpoints=size(c_real,2);
gabor.f = linspace(0 , 1, gabor.fpoints);
gabor.t = (0:gabor.tpoints-1)*gabor.a  ;
ff=round(gabor.fpoints*topfreq/0.5);
nam=['sparsity ratio',num2str(Sparsityratio)];
plotTFRsimple(abs(c_real(1:ff,:)),gabor.t,gabor.f(1:ff),FS,0.2,nam)

        
        %% RECONSTRUCT via NESTA, using "ANALYSIS"
        opts = [];
        opts.U = psiT_A;
        opts.Ut = psi_A;
        opts.Verbose = 10;
        
        b = A(x_exact);     % the data
        delta = 1e-11;
        muf = 1e-7;
        disp('------- Reconstruction via Analysis -----------');
        tic
        [xk,niter,resid,outData] = NESTA(A,At,b,muf,delta,opts);
        time.analysis = toc;
        
%         % RECONSTRUCT via NESTA, using "SYNTHESIS"
%         
        opts = [];
        opts.Verbose = 10;
        AA = @(x) A( psi_A( x ) );
        AAt =@(y) psiT_A( At( y ) );
        
        delta = 1e-11;
        mu = 1e-7;
        disp('------- Reconstruction via Synthesis -----------');
        tic
        [coeff_k,niter,resid,outData] = NESTA(AA,AAt,b,muf,delta,opts);
        xk2 = psi_A(coeff_k);
        time.synthesis = toc;
     %   coeff_k=zeros(size(xk));

        %% Figure 10.18(a)
        figure; clf;
        %subplot(2,2,1);
        coeff = psiT_A( xk );
        semilogy( sort(abs(coeff),'descend') )
        hold all
        semilogy( sort(abs(coeff_k),'descend') );
        semilogy( sort(abs(psiT_A(x_exact)),'descend') );
        title('Dictionary coefficients of recovered signal');
        xlabel('Sorted coefficients');
        ylabel('Magnitude');
        legend('recovered via analysis','recovered via synthesis',...
            'original signal');
        
  
%  Compare the original bat signal and the reconstructed one. 
figure
plot( (1:N)/FS, x_exact,'k');
hold on
plot( (1:N)/FS, xk, '--r' );
%plot( (1:N)/FS, xk2, '--g' );


% Figure 10.18 (b)
[s,t,f]=STFT(xk,Nw,2*Nf,Ntpl);
[l,c]=size(s);
%topfreq=0.3;
ff=round(Nf*topfreq/0.5);
nam=['SPEC, W=',num2str(Nw)];
plotTFRsimple(s(1:ff,:),t,f(1:ff),FS,0.2,nam)

% Figure 10.18 (c)
[s,t,f]=STFT(xk2,Nw,2*Nf,Ntpl);
[l,c]=size(s);
%topfreq=0.3;
ff=round(Nf*topfreq/0.5);
nam=['SPEC, W=',num2str(Nw)];
plotTFRsimple(s(1:ff,:),t,f(1:ff),FS,0.2,nam)


%% Asuming the singal is sparse in Gabor - High redundancy        
 % Construction of GAbor window and its dual       
 Time_samples=2^8;
 Freq_shift=2^5;
 gabor.a = N / Time_samples;  %  time shift
 gabor.M = N / Freq_shift; % frequency bands
 
 % For this Gabor system, the optimally concentrated Gaussian is given by
 tf_ratio = gabor.a/Freq_shift;
 
 % Compute the window and the canonical dual window
 gabor.w = pgauss(N,tf_ratio*4);
 gabor.wd=gabdual(gabor.w,gabor.a,gabor.M);
 
           
yy = dgt(ones(size(x_exact)),gabor.w,gabor.a,gabor.M); 
gabor.fpoints=size(yy,1);
gabor.tpoints=size(yy,2);
clear yy
psi_A = @(x) op_igabor(x,gabor,N); 
psiT_A = @(y) op_gabor(y,gabor);            

c = psiT_A(x_exact);
cc=sort(abs(c).^2,'descend');
tot_en = sum(cc);
tot99=tot_en*99/100;
cc_cum = cumsum(cc);
Sparsity = length(cc)-length(find(cc_cum>tot99));
Sparsityratio = Sparsity/length(cc);


        
        %% RECONSTRUCT via NESTA, using "ANALYSIS"
        opts = [];
        opts.U = psiT_A;
        opts.Ut = psi_A;
        opts.Verbose = 10;
        
        b = A(x_exact);     % the data
        delta = 1e-11;
        muf = 1e-7;
        disp('------- Reconstruction via Analysis -----------');
        tic
        [xk,niter,resid,outData] = NESTA(A,At,b,muf,delta,opts);
        time.analysis = toc;
        

[s,t,f]=STFT(xk,Nw,2*Nf,Ntpl);
[l,c]=size(s);
ff=round(Nf*topfreq/0.5);
nam=['SPEC, W=',num2str(Nw)];
plotTFRsimple(s(1:ff,:),t,f(1:ff),FS,0.1,nam)




end

function [s,t,f]=STFT(signal,Nw,Nf,Nt)
%Short Time Fourier Transform (STFT).

%--- Input ---
%Signal: Vector containing the signal samples.
%Nw: window length
%Nf: number of frequency points
%Nt: number of time points

% --- Output ---
% s: Matrix of the time frequency representation.
% t: Time instances (corresponding to the columns of s)
% F: Vector of normalized frequencies (corresponding to the lines of s)

ww=tftb_window(Nw,'Gauss');
[s,t,f]=tfrstft(signal,1:round(length(signal)/Nt):length(signal),Nf,ww,1);
s=abs(s(1:Nf/2,:));
f=f(1:Nf/2);
end


function plotTFRsimple(s,t,f,Fs,dbval,nam)
% Plots a time frequency representation.
% s: Matrix of the time frequency representation.
% t: Time instances (corresponding to the columns of s)
% F: Vector of normalized frequencies (corresponding to the lines of s)
% Fs: Sampling Frequency 
% dbval: Value determining the initial dynamic range. After the plot, the
%        dynamic range can be modified with slider existing in the produced
%        figure. Reducing this value enhances the visibility of the lower energy 
%        components of the time frequency representations.
% nam: The name that is given to the produced figure.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The script is written by Yannis Kopsinis & Elias Aboutanios for demonstration 
% purposes. Therefore, it might not be suitable in general signal cases and
% it is not optimised with respect to speed.

% For questions and/or corrections email at kopsinis@ieee.org.
% Feel free to use and modify the specific script. However, in case of
% publication the citation of the following paper will be greatly
% appreciated:

% Y. Kopsinis, E. Aboutanios, D. A. Waters, S. McLaughlin, "Instantaneous
% Frequency - based Time-Frequency representation of multi-harmonic
% signals," Journal of Acoustical Society of America (JASA),
% vol. 127(2), Feb 2010.

data.s=s;
data.t=t;
data.f=f;
data.Fs=Fs;
if nargin==5
fh = figure('Position',[250 250 350 350]);
else
    fh = figure('Position',[250 250 350 350],'name',nam);
end
ax1 = axes(...    % Axes for plotting the selected plot
    'Parent', fh, ...
    'Units', 'normalized', ...
    'HandleVisibility','callback', ...
    'Position',[0.11 0.13 0.80 0.67]);
% sh = uicontrol(fh,'Style','slider',...
%     'Max',60,'Min',0,'Value',dbval,...
%     'SliderStep',[0.01 0.1],...
%     'units','normalized',...
%     'Position',[0.11 0.85 0.70 0.05],...
%     'Callback',@slider_callback);
% eth = uicontrol(fh,'Style','edit',...
%     'String',num2str(get(sh,'Value')),...
%     'units','normalized',...
%     'Position',[0.81 0.85 0.1 0.05],...
%     'Callback',@edittext_callback);
% number_errors = 0;

% Set edit text UserData property to slider structure.

data.sliderval=dbval;

thres=data.sliderval;
maxi=max(max(s));
mini=max(min(min(s)),maxi*thres/100.0);
indmin=find(s<mini);
s(indmin)=mini*ones(1,length(indmin));
indmax=find(s>maxi);
s(indmax)=maxi*ones(1,length(indmax));

imagesc(t/Fs,f*Fs,10*log10(s),'parent',ax1)
data.typeofplot='imagesc';
axis('xy')
clmp=colormap('hot');
%clmp=colormap('pink');

colormap(flipud(clmp));
set(fh,'UserData',data)
ylabel('Frequency (KHz)')
xlabel('Time')

end

function y = op_gabor(x,gabor)


      y = dgt(x,gabor.w,gabor.a,gabor.M); 
      y = y(:);
end

function y = op_igabor(x,gabor,N)

      x = reshape(x, gabor.fpoints, gabor.tpoints);
      y = real(idgt(x,gabor.wd,gabor.a,N)); 
end
