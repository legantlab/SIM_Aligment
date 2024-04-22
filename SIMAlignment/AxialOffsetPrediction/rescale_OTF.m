function O_scaled = rescale_OTF(O,pxl_dim_PSF,data_dim,pxl_dim_data)
% This function rescale pixel size in the calibrated OTF to match the acquired image
% O: calibrated OTF
% pxl_dim_PSF: pixel size of calibrated PSF
% data_dim: size of acquired image in [y,x,z]
% pxl_dim_data: pixel size of acquired image

[ny_PSF,nx_PSF,nz_PSF,norders] = size(O);
ny_data = data_dim(1);
nx_data = data_dim(2);
nz_data = data_dim(3);
dk_PSF=1./([nx_PSF,ny_PSF,nz_PSF].*pxl_dim_PSF);
dk_data=1./([nx_data,ny_data,nz_data].*pxl_dim_data);


OTF_xx=[-ceil((nx_PSF-1)/2):floor((nx_PSF-1)/2)]*dk_PSF(1);
OTF_yy=[-ceil((ny_PSF-1)/2):floor((ny_PSF-1)/2)]*dk_PSF(2);
OTF_zz=[-ceil((nz_PSF-1)/2):floor((nz_PSF-1)/2)]*dk_PSF(3);

map_xx=[-ceil((nx_data-1)/2):floor((nx_data-1)/2)]*dk_data(1);
map_yy=[-ceil((ny_data-1)/2):floor((ny_data-1)/2)]*dk_data(2);
map_zz=[-ceil((nz_data-1)/2):floor((nz_data-1)/2)]*dk_data(3);

[OTF_xx_arr,OTF_yy_arr,OTF_zz_arr]=meshgrid(OTF_xx,OTF_yy,OTF_zz);
[map_xx_arr,map_yy_arr,map_zz_arr]=meshgrid(map_xx,map_yy,map_zz);
O_scaled=zeros(size(map_xx_arr));
mask_scaled=zeros(size(map_xx_arr));

for ii=1:norders
    O_scaled(:,:,:,ii) = interp3(OTF_xx_arr,OTF_yy_arr,OTF_zz_arr,O(:,:,:,ii),map_xx_arr,map_yy_arr,map_zz_arr,'spline',0+0i);
end
