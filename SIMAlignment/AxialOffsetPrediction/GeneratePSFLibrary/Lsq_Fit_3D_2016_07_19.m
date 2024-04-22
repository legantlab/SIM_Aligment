function FitData = Lsq_Fit_3D_2016_07_19_complex(PeakImg)
%Performs least square fitting to a 3D stack of a bead
%Determines centroid position in units of pixels

[nx,ny,nz]=size(PeakImg);
%Generate pixel coordinates for Gaussian fitting
[xx,yy,zz]=meshgrid(1:nx,1:ny,1:nz);
PeakData=[PeakImg(:),xx(:),yy(:),zz(:)];

%Use the value and location of the max voxel as initial guesses
[maxval,ind]=max(PeakImg(:));
[cx,cy,cz]=ind2sub([nx,ny,nz],ind);

xo(1) = maxval;  %N(guess) = max value within the window
xo(2) = cx;  %xo(guess) = x pixel value of peak pixel in the window
xo(3) = cy;  %yo(guess) = y pixel value of peak pixel in the window
xo(4) = cz;  %zo(guess) = z pixel value of peak pixel in the window
xo(5) = 1.5; %sigx(guess) = x sigma value
xo(6) = 1.5; %sigy(guess) = y sigma value
xo(7) = 1.5; %sigz(guess) = z sigma value
% xo(8)= 0;  %b(guess) = background photons per pixel

%find the initial best fit with unadjusted weights:
options = optimset('Display', 'off', 'MaxIter', 10000, 'MaxFunEvals', 10000, 'TolX', 1e-6, 'TolFun', 1e-6);
[x, resnorm] = fminsearch(@(x) objFun_Gaussian_Fit_3D_2016_07_19(x,PeakData), xo,options);

FitData=x;
% 
F=Gaussian_Fit_3D_2016_07_19(x,PeakData,[nx,ny,nz]);
% figure
% subplot(2,2,1)
% imagesc(squeeze(max(reshape(PeakData(:,1),[nx,ny,nz]),[],3)))
% axis equal
% title('XY projection of the input image')
% subplot(2,2,2)
% imagesc(squeeze(max(reshape(PeakData(:,1),[nx,ny,nz]),[],1))')
% axis equal
% title('XZ projection of the input image')
% subplot(2,2,3)
% imagesc(squeeze(max(F,[],3)))
% axis equal
% title('XY projection of the fit image')
% subplot(2,2,4)
% imagesc(squeeze(max(F,[],1))')
% axis equal
% title('XZ projection of the fit image')
