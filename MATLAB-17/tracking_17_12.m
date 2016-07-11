%-----------------------------------------------------------------
%  Exercise 15.11
%  Visual Tracking
%  Use particle filter for tracking a circle moving according to 
%  a random uniform model in the image.
%
%  Includes functions for circle plotting, likelihood estimation 
%  and particle sampling.
%
%  Writing to a video file can be activated by activating the 
%  related commented code
%-----------------------------------------------------------------

function tracking_15_11

clear

DIM = 300;%image size (DIM x DIM)

CYCLES = 120;%the number of cycles

r=10;%the radius of the circle to track

C=5;%coefficient for the noise covariance
    %bigger C => higher dispersion of particles
    
%initial position of the circle at the image center
i0 = floor(DIM/2);
j0 = floor(DIM/2);

%mean and variance for the noise model 
M = [0 0]';
V = C*[2 0.5; 0.5 2];

N = 50; %number of particles

%initialize particles
particle =  mvnrnd(M,V,N) + repmat(  [i0 j0], N, 1);

% initialize particle weights to the same value
w = ones(1,N)/N; 

%activate commented 2 lines to write into a video file
%writerObj = VideoWriter('particlefilter.avi');
%open(writerObj);


for k=1:CYCLES

    % the circle moves randomly following a uniform
    % motion model in the interval [-10 10] for each dimension
    eta = floor(20*rand(2,1)-10);  
    i0 = eta(1) + i0;
    j0 = eta(2) + j0;

    %the real circle image
    target_image = plotCircle(r,i0,j0, DIM); 
    
    %add noise to image
    %target_image = imnoise(target_image,'salt & pepper', 0.04);
    
    %get the likelihood for each particle
    for c=1:N
        w(c) = calculateObservationLikelihood( i0, j0, particle(c,1), particle(c,2) );
    end

    %weight normalization
    if sum(w)==0
        break;
    end
    w = w/sum(w);
    
    %create image for plotting 
    particle_image = zeros(DIM,DIM);

    %Draw all particles into an RGB image with red color
    for c=1:N
        particle_image = particle_image + plotCircle(r,particle(c,1),particle(c,2), DIM);
    end
    
    %do resampling based on the weights
    for c=1:N
        s = performSampling(w);
        M1 = [particle(s,1); particle(s,2)]; 
        y(c,:) = M1 + mvnrnd(M,V,1)'; %generate observations
    end
    
    particle = y;

    rgb_image = ones(DIM, DIM, 3);
    
    %paint red the particles
    rgb_image(:,:,2) = rgb_image(:,:,1) .* double(1 - particle_image);
    rgb_image(:,:,3) = rgb_image(:,:,2) ;
    
    %paint black the target
    rgb_image(:,:,1) = rgb_image(:,:,1).*double(1-target_image);
    rgb_image(:,:,2) = rgb_image(:,:,2).*double(1-target_image);
    rgb_image(:,:,3) = rgb_image(:,:,2);

    imshow(rgb_image);
    pause(0.25) %delay for display purposes
        
    %activate commented line to write into a video file
    %writeVideo(writerObj,rgb_image);
end

%activate commented line to write into a video file
%close(writerObj);
end


%Plot a white circle onto an image with black background
%r the circle radius
%i0, j0 the circle center
%DIM the size of the rectangular image
function outimage = plotCircle(r,i0,j0, DIM)
t = 0:0.1:2*pi;

i = round (r * cos(t) + i0);
j = round (r * sin(t) + j0);

outimage = zeros(DIM,DIM);

valid = (i>0 & i<DIM & j>0 & j<DIM);

i = i(valid);
j = j(valid);

for c=1:size(i,2)
    outimage(i(c),j(c))=1;
end
end

%Calculate likelihood p(y_n|x_n) of a particle based on the distance from the real center
% i0, j0 (x_n) the circle center coordinates
% i, j (y_n) the particle coordinates
function f = calculateObservationLikelihood(i0, j0, i, j)

a = 2; %decay factor

d = sqrt( (i-i0)^2 + (j-j0)^2 ); 

f = exp(-a * d);
end

%Do sampling given the weight function w
function out = performSampling(w)

x = rand;
acc = 0;
i=1;

while 1
    
    acc = acc + w(i);
    if acc > x
        break;
    end
    i=i+1;
    
end
out=i;

end
