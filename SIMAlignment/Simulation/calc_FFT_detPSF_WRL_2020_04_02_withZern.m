function [PSFa,PSFi,OTF,kx,ky,kz]=calc_FFT_detPSF_WRL_2020_04_02_withZern(NA,F,nimm,wvl,nx,ny,nz,dx,dy,dz,offset,tip_tilt,piston,zernNum,zernWeight)
%[PSFa,OTF]=calc_FFT_PSF(1.27,1.33,1,512,512,512,.01,.01);

%Use FFT's to Calculate the point spread function and optical transfer function for a
%back pupil function using the Debye approximation

%NA is the numerical aperature of the objective lens

%F is the focal length of the objective lens in microns

%nimm is the refractive index of the immersion media

%wvl is the wavelength in microns

%nx, ny, and nz are the extent of the volume in which the PSF is calculated in microns

%dx == dy, and dz are the x/y and z voxel dimensions of the PSF grid in microns

%offset is a two element vector specifying the x and y offset of the collimated excitation beam 
%at the back pupil with respect to the optical axis of the objective (in microns)

%tip_tilt is a two element vector specifying the tip/tilt amplitude of the
%collimated excitation beam at the back pupil.  The input values are angles
%in degrees.

%piston is the relative path length difference between the opposing lenses

%zernNum is a 1x2 vector of Zernike polynomial aberrations to be added to
%the pupil function.

%zernWeight is a 1x2 vector with the weights of the Zernike functions - weight describes the max value of the aberration in wavelengths.

%Incorporates vectorial effects.
dkx=1/(nx*dx); %Pixel dimensions in inverse space
dky=1/(ny*dy);
dkz=1/(nz*dz);

alpha=asin(NA/nimm); %Half angle of acceptance

kmag=nimm/wvl; %Magnitude of the wave vector in microns

kx_offset=offset(1)/(wvl*F); %Convert from offset at physical back pupil dimensions to offset in k-space
ky_offset=offset(2)/(wvl*F); %Convert from offset at physical back pupil dimensions to offset in k-space

kx=([-(nx-1)/2:(nx-1)/2])*dkx; %Populate k-space arrs
ky=([-(ny-1)/2:(ny-1)/2])*dky;
kz=([-(nz-1)/2:(nz-1)/2])*dkz;
[kx_arr,ky_arr]=meshgrid(kx,ky);
kx_arr=single(kx_arr);
ky_arr=single(ky_arr);
kr_arr=sqrt(kx_arr.^2+ky_arr.^2);
kz_arr=(kr_arr<kmag).*sqrt(kmag.^2-kr_arr.^2); %kz(kx,ky)

theta_arr=(kr_arr<kmag).*asin(kr_arr/kmag); %Incident angle
phi_arr=atan2(ky_arr,kx_arr); %Azimuthal angle

%Construct the pupil function
krmax=NA/wvl;
pupil=1./sqrt(kz_arr/kmag);
pupil(isinf(pupil))=0;
kr_arr2=sqrt((kx_arr+kx_offset).^2+(ky_arr+ky_offset).^2);
pupil_mask=kr_arr2<(kmag*sin(alpha)); %Mask the back pupil according to the objective NA
pupil=pupil.*pupil_mask;

phase_ramp=exp(1i*2*pi/wvl*(kx_arr*wvl*F*tand(tip_tilt(1))+ky_arr*wvl*F*tand(tip_tilt(2))));

pupil=pupil.*phase_ramp.*exp(1i*2*pi);
%%
%Compute Zernike aberrations for the given NA and PSF dimensions
kr_norm=kr_arr2/krmax;
is_in_circle=kr_norm<1;

N = []; M = [];
for n = 0:9
N = [N n*ones(1,n+1)];
M = [M -n:2:n];
end 
Z = zernfun(N,M,kr_norm(is_in_circle),phi_arr(is_in_circle));

%Apply Zernike aberrations to the pupil function
for pp=1:length(zernNum)
pupil(is_in_circle)=pupil(is_in_circle).*exp(1i*2*pi*zernWeight(pp).*Z(:,zernNum(pp)+1));
end

%%
%Introduce vectorial factors
Px1=cos(theta_arr).*cos(phi_arr).^2+sin(phi_arr).^2;
Px2=(cos(theta_arr)-1).*sin(phi_arr).*cos(phi_arr);
Py1=(cos(theta_arr)-1).*sin(phi_arr).*cos(phi_arr);
Py2=cos(theta_arr).*sin(phi_arr).^2+cos(phi_arr).^2;    
Pz1=sin(theta_arr).*cos(phi_arr);
Pz2=sin(theta_arr).*sin(phi_arr);

PSFa=zeros(6,ny,nx,nz,'single');
PSFi=zeros(ny,nx,nz,'single');

t=1;
for z=[-(nz-1)/2:(nz-1)/2]
    kz_shell=exp(1i*2*pi*kz_arr*z*dz);
    PSFa(1,:,:,t)=fftshift(ifft2(ifftshift(pupil.*kz_shell.*Px1))); %X component of the electric field due to x polarized light at back pupil
    PSFa(2,:,:,t)=fftshift(ifft2(ifftshift(pupil.*kz_shell.*Py1))); %Y component of the electric field due to x polarized light at back pupil
    PSFa(3,:,:,t)=fftshift(ifft2(ifftshift(pupil.*kz_shell.*Pz1))); %Z component of the electric field due to x polarized light at back pupil
    PSFa(4,:,:,t)=fftshift(ifft2(ifftshift(pupil.*kz_shell.*Px2))); %X component of the electric field due to y polarized light at back pupil
    PSFa(5,:,:,t)=fftshift(ifft2(ifftshift(pupil.*kz_shell.*Py2))); %Y component of the electric field due to y polarized light at back pupil
    PSFa(6,:,:,t)=fftshift(ifft2(ifftshift(pupil.*kz_shell.*Pz2))); %Z component of the electric field due to y polarized light at back pupil
    t;
t=t+1;
end
    
for ii=1:6
    PSFi=PSFi+squeeze(PSFa(ii,:,:,:).*conj(PSFa(ii,:,:,:)));
end
    
OTF=fftshift(fftn(ifftshift(PSFi)));
