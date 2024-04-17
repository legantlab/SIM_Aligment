function [PSFa,PSFi,OTF,xx,yy,zz,kx,ky,kz,PW_scaled,pupil]=calc_FFT_extPSF_HEX_2017_03_28(NA,F,nimm,wvl,nx,ny,nz,dx,dy,dz,offset,tip_tilt,piston,input_pol,phase)
%[PSFa,OTF]=calc_FFT_PSF(1.27,1.33,1,512,512,512,.01,.01);

%Use FFT's to calculate the point spread function and optical transfer function for a
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

%input_pol is either 'x' or 'y' corresponding to linear polarization of the
%pupil function

%phase - is a lateral shift applied to the illumination pattern at the sample given as a
%2 pi fraction of the highest lateral frequency in the illumination patter.  Useful for specifying phase steps for strucured illumination microscopy.

alpha_max=asin(NA(1)/nimm); %max angular extent of the input mask
alpha_min=asin(NA(2)/nimm); %min angular extent of the input mask
alpha_cent=asin(NA(3)/nimm); %cone angle for input illumination beams

%Generate wavevectors for hexagonal lattice
PW=zeros(6,3);

for ii=1:6
    theta=2*pi/6*(ii-1);
    PW(ii,1:2)=[0,1]*[cos(theta),-sin(theta);sin(theta),cos(theta)];
end

%Restrict wavevectors to a cone defined by the central NA
PW(:,1:2) = sin(alpha_cent) .* PW(:,1:2);
PW(:,3) = cos(alpha_cent);

% %Normalize wavevectors to unit length
% PW(1,:)=PW(1,:)/norm(PW(1,:));
% PW(2,:)=PW(2,:)/norm(PW(2,:));
% PW(3,:)=PW(3,:)/norm(PW(3,:));

[nvec,~]=size(PW);

%Incorporates vectorial effects.
dkx=1/(nx*dx); %Pixel dimensions in inverse space
dky=1/(ny*dy);
dkz=1/(nz*dz);

alpha_vec=asin(PW(:,1)); %an array of the lateral (x) angles of the excitation pupil plane
kmag=nimm/wvl; %Magnitude of the wave vector in microns
PW_scaled=PW*kmag;


kx_offset=offset(1)/(wvl*F); %Convert from offset at physical back pupil dimensions to offset in k-space
ky_offset=offset(2)/(wvl*F); %Convert from offset at physical back pupil dimensions to offset in k-space

xx=([-(nx-1)/2:(nx-1)/2])*dx; %Populate real-space vectors
yy=([-(ny-1)/2:(ny-1)/2])*dy;
zz=([-(nz-1)/2:(nz-1)/2])*dz;

kx=([-(nx-1)/2:(nx-1)/2])*dkx; %Populate k-space arrs
ky=([-(ny-1)/2:(ny-1)/2])*dky;
kz=([-(nz-1)/2:(nz-1)/2])*dkz;
[kx_arr,ky_arr]=meshgrid(ky,kx); %Note shift here to account for meshgrid swapping dimensions
kr_arr=sqrt(kx_arr.^2+ky_arr.^2);
kz_arr=(kr_arr<kmag).*sqrt(kmag.^2-kr_arr.^2); %kz(kx,ky)

theta_arr=(kr_arr<kmag).*asin(kr_arr/kmag); %Incident angle
phi_arr=atan2(ky_arr,kx_arr); %Azimuthal angle

%Construct the pupil function
pupil=1./sqrt(kz_arr/kmag);
pupil(isinf(pupil))=0;
kr_arr=sqrt((kx_arr+kx_offset).^2+(ky_arr+ky_offset).^2);
pupil_mask=0*pupil;
for ii=1:nvec
pupil_mask=pupil_mask|(kr_arr<(kmag*sin(alpha_max))&kr_arr>=(kmag*sin(alpha_min))&kx_arr>=(kmag*(sin(alpha_vec(ii)))-dkx/2)&kx_arr<=(kmag*(sin(alpha_vec(ii)))+dkx/2)); %Mask the back pupil according lattice vectors and mask NA
end
pupil=pupil.*pupil_mask;

phase_ramp=exp(1i*2*pi*F*(kx_arr*sind(tip_tilt(1))+ky_arr*sind(tip_tilt(2))));

deltas=find(sum(pupil_mask,1)); %Find the columns of the pupil where the hex spots are located
deltas=deltas(deltas>(ny+1)/2); %Just look at half of the image
delta_max_period=min(deltas); %Find the maximum period - will be used to compute the phase steps

%Apply lateral phase offset based on input parameter "phase"
tilt_step=asind(((phase)/(2*pi))/(F*ky(delta_max_period)));
phase_step=exp(1i*2*pi*F*(kx_arr*sind(tilt_step)));
pupil=pupil.*phase_ramp.*exp(1i*2*pi*piston).*phase_step;

%Introduce vectorial factors
Px1=cos(theta_arr).*cos(phi_arr).^2+sin(phi_arr).^2;
Px2=(cos(theta_arr)-1).*sin(phi_arr).*cos(phi_arr);
Py1=(cos(theta_arr)-1).*sin(phi_arr).*cos(phi_arr);
Py2=cos(theta_arr).*sin(phi_arr).^2+cos(phi_arr).^2;
Pz1=sin(theta_arr).*cos(phi_arr);
Pz2=sin(theta_arr).*sin(phi_arr);

t=1;
for z=[-(nz-1)/2:(nz-1)/2]
    kz_shell=exp(1i*2*pi*kz_arr*z*dz);
    if strcmp(input_pol(1),'x')
        PSFa(1,:,:,t)=fftshift(ifft2(ifftshift(pupil.*kz_shell.*Px1))); %X component of the electric field due to x polarized light at back pupil
        PSFa(2,:,:,t)=fftshift(ifft2(ifftshift(pupil.*kz_shell.*Py1))); %Y component of the electric field due to x polarized light at back pupil
        PSFa(3,:,:,t)=fftshift(ifft2(ifftshift(pupil.*kz_shell.*Pz1))); %Z component of the electric field due to x polarized light at back pupil
    end
    if strcmp(input_pol(1),'y')
        PSFa(1,:,:,t)=fftshift(ifft2(ifftshift(pupil.*kz_shell.*Px2))); %X component of the electric field due to y polarized light at back pupil
        PSFa(2,:,:,t)=fftshift(ifft2(ifftshift(pupil.*kz_shell.*Py2))); %Y component of the electric field due to y polarized light at back pupil
        PSFa(3,:,:,t)=fftshift(ifft2(ifftshift(pupil.*kz_shell.*Pz2))); %Z component of the electric field due to y polarized light at back pupil
    end
    t;
    t=t+1;
end

if strcmp(input_pol(1),'x')
    PSFi_x=zeros(nx,ny,nz);
    for ii=1:3
        PSFi_x=PSFi_x+squeeze(PSFa(ii,:,:,:).*conj(PSFa(ii,:,:,:)));  %Component due to x-polarized light at back pupil
    end
    PSFi(1,:,:,:)=PSFi_x;
end
if strcmp(input_pol(1),'y')
    PSFi_y=zeros(nx,ny,nz);
    for ii=1:3
        PSFi_y=PSFi_y+squeeze(PSFa(ii,:,:,:).*conj(PSFa(ii,:,:,:)));  %Component due to y-polarized light at back pupil
    end
    PSFi(1,:,:,:)=PSFi_y;
end
if (length(input_pol)==2 && strcmp(input_pol(2),'y'))
    PSFi_y=zeros(nx,ny,nz);
    for ii=4:6
        PSFi_y=PSFi_y+squeeze(PSFa(ii,:,:,:).*conj(PSFa(ii,:,:,:)));  %Component due to y-polarized light at back pupil
    end
    PSFi(2,:,:,:)=PSFi_y;
end

if strcmp(input_pol(1),'x')
    OTF(1,:,:,:)=fftshift(fftn(ifftshift(PSFi_x)));
end
if strcmp(input_pol(1),'y')
    OTF(1,:,:,:)=fftshift(fftn(ifftshift(PSFi_y)));
end
if (length(input_pol)==2 && strcmp(input_pol(2),'y'))
    OTF(2,:,:,:)=fftshift(fftn(ifftshift(PSFi_y)));
end
        

