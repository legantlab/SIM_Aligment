function [shiftedImage]=ifourierShift3D(inputImage,shifts)
%Performs sub-pixel real-space shifts of a 3D image by applying a phase
%ramp in Fourier space

%Inputs:
%inputImage is a 3D array of either real or complex values
%shifts is a 1 x 3 vector of shift values in real-space pixels

%Set up vectors based on image size
[nx,ny,nz]=size(inputImage);
kxx=[-ceil((nx-1)/2):floor((nx-1)/2)];
kyy=[-ceil((ny-1)/2):floor((ny-1)/2)];
kzz=[-ceil((nz-1)/2):floor((nz-1)/2)];

%Set up arrays for frequency space pixels - note the switch in meshgrid
%inputs between x and y ordering. 
[kxx_arr,kyy_arr,kzz_arr]=meshgrid(kyy/ny,kxx/nx,kzz/nz);

%FFT image, apply phase ramp, and inverse transform
Fimage=fftshift(ifftn((inputImage)));
Fimage_shift=Fimage.*exp(1i*2*pi*(kxx_arr.*shifts(1)+kyy_arr.*shifts(2)+kzz_arr.*shifts(3)));
shiftedImage=fftn(ifftshift(Fimage_shift));

end