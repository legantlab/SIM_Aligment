function [OTF_extent]=calc_OTF_extent(NA,nimm,wvl,nz,ny,nx,dx,dz,shift)
%Compute OTF extent by computing the overlap of shifted toroids

%NA is the numerical aperature of the objective lens

%nimm is the refractive index of the immersion media

%wvl is the wavelength in microns

%nx, ny, and nz are the extent of the volume in which the PSF is calculated in microns

%dx == dy, and dz are the x/y and z voxel dimensions of the PSF grid in microns

%shift is an 3 x m position vectors specifying the locations of the shifted
%copies in z-pixels

dkx=1/(nx*dx); %Pixel dimensions in inverse space
dky=1/(ny*dx);
dkz=1/(nz*dz);

alpha=asin(NA/nimm); %Half angle of acceptance

kmag=nimm/wvl; %Magnitude of the wave vector in microns

kx=([-(nx-1)/2:(nx-1)/2])*dkx; %Populate k-space arrs
ky=([-(ny-1)/2:(ny-1)/2])*dky;
kz=([-(nz-1)/2:(nz-1)/2])*dkz;
[kx_arr,ky_arr,kz_arr]=meshgrid(kx,ky,kz);

R=kmag*sin(alpha); %radius of toroidal circle of rotation about the kz axis
a=kmag; %radius of toroidal cross-secton in the kz plane
kz_offset=kmag*cos(alpha);

OTF_extent=zeros(nx,ny,nz);

[num_positions,~]=size(shift);

for ii=1:num_positions
    %Note that the extra factor of dkz/2 was added to account for potential
    %pixel discretization errors
T1=((sqrt((kx_arr-shift(ii,1)*dkx).^2+(ky_arr-shift(ii,2)*dky).^2)-R).^2/a^2+(kz_arr-kz_offset-shift(ii,3)*dkz).^2/a^2)<=(1+dkz/2);
T2=((sqrt((kx_arr-shift(ii,1)*dkx).^2+(ky_arr-shift(ii,2)*dky).^2)-R).^2/a^2+(kz_arr+kz_offset-shift(ii,3)*dkz).^2/a^2)<=(1+dkz/2);
OTF_extent=OTF_extent|(T1&T2);
end