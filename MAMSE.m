function [se,sd,en]=MAMSE(img,angV,scales,m,r)
%MAXIMUM ANGULAR MULTISCALESAMPLEENTROPY
%
% se = Maxumum Angular Sample entropy for each scales
% sd = difference between max and min sample entropy
% en = Matrix with all data. Columns are angles, lines are scales, cells
% are mean sample entropy over all image profiles.
% This code was implemented by José Garcia Vivas Miranda on 01 May 2022.
% Contact: vivas@ufba.br

t = zeros(scales, 1);   % store the sample entropy values
en=t;                   % store the mean sample entropy values
szO = length(img);
sz=floor((szO/2)*1.41);        % The size of the máximum square within the rotated image
col=1;

for ang=angV
    disp(['calc angle' num2str(ang)]);
    imgR=imrotate(img,ang,'nearest','crop');  % Rotate original images
    ul=szO/2-sz/2;
    rc = [ul,ul,sz-1,sz-1];
    imgR=imcrop(imgR,rc);
    n=length(imgR);
    
    for tau = 1:scales      % Loop for scales
        
        for j = 1:n         % vertical profiles
            
            signal = imgR(:,j)';
            
            [ t(tau, j)] = multiscaleSampleEntropy( signal, m, r, tau );           
        end
    end
    en(:,col)=mean(t, 2); 
    col=col+1;
end
se=nanmax(en,[],2)';
sm=nanmin(en,[],2)';
sd=se-sm;
end


function [ e, A, B ] = multiscaleSampleEntropy( x, m, r, tau )
%MULTISCALESAMPLEENTROPY
%
% Based on "Multiscale entropy analysis of biological signals"
% By Madalena Costa, Ary L. Goldberger, and C.-K. Peng
% Published on 18 February 2005 in Phys. Rev. E 71, 021906.
% And the code implemented by John Malik on 26 April 2017.
% Contact: john.malik@duke.edu
switch nargin
    case 1
        m = 2;
        r = 0.15;
        tau = 1;
    case 2
        r = 0.15;
        tau = 1;
    case 3
        tau = 1;
end
% coarse signal
y = mean(buffer(x(:), tau), 1);
% (m+1)-element sequences
X = buffer(y, m + 1, m, 'nodelay')';
% matching (m+1)-element sequences
A = sum(pdist(X, 'chebychev') < r * nanstd(x, 1));
% matching m-element sequences
X = X(:, 1:m);
B = sum(pdist(X, 'chebychev') < r * nanstd(x, 1));
% take log
if A == 0 || B == 0
    e = NaN;
    return
end
e = log(B / A);
end
