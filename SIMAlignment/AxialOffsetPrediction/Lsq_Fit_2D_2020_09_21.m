function FitData = Lsq_Fit_2D_2020_09_21(PeakImg,r_max,c_max)
%Performs least square fitting to a 2D stack of a bead
%Determines centroid position in units of pixels

A=size(PeakImg);
[xx,yy]=meshgrid((-A(2)+1)/2:(A(2)-1)/2,(-A(1)+1)/2:(A(1)-1)/2);
PeakData=[PeakImg(:),xx(:),yy(:)];
xo(1) = max(PeakImg(:));  %N(guess) = max value within the window
xo(2) = c_max;  %xo(guess) = x pixel value of peak pixel in the window
xo(3) = r_max;  %yo(guess) = y pixel value of peak pixel in the window
xo(4) = 1.5; %sigx(guess) = x sigma value
xo(5) = 1.5; %sigy(guess) = y sigma value
% xo(8)= 0;  %b(guess) = background photons per pixel

%find the initial best fit with unadjusted weights:
options = optimset('Display', 'off', 'MaxIter', 10000, 'MaxFunEvals', 10000, 'TolX', 1e-6, 'TolFun', 1e-6);
[x, resnorm] = fminsearch(@(x) objFun_Gaussian_Fit_2D_2020_09_21(x,PeakData), xo,options);

FitData=x;
% 
F=Gaussian_Fit_2D_2020_09_21(x,PeakData,A);
figure
subplot(2,1,1)
%imagesc(squeeze(max(reshape(PeakData(:,1),A),[],3)))
imagesc(reshape(PeakData(:,1),A))
axis equal
title('XY projection of the CC-peak')
subplot(2,1,2)
imagesc(F)
axis equal
title('XY projection of the fitted CC-peak')
