function O = seperate_sample(PSF_filename, n_phase)
% This function separate different frequency components of acquired image
% PSF_filename: path to the image
% n_phase: number of phases

%% read image stack
%PSF_filename = 'Z:\Yu\2020_05_21(LLS_SIM)\NAp55p50\RAW_z_stack_latticeSIM_z100nm_dphi_234nm_6_CamB_ch0_CAM1_stack0000_560nm_0000000msec_0075656857msecAbs_000x_000y_000z_0000t.tif';
img_num = length(imfinfo(PSF_filename));
img_stack = imread(PSF_filename,1);
for i = 2:img_num
    img_temp = imread(PSF_filename,i);
    img_stack = cat(3,img_stack,img_temp);
end

%% subtract background
background = mode(img_stack,'all')*1.5;
nz = img_num/n_phase;
nx = size(img_temp,2);
ny = size(img_temp,1);

%% apply shift to PSF to move beads to the center    
for i = 1:n_phase
    Dr_temp = img_stack(:,:,i:n_phase:end) - background;
    Dr_temp(Dr_temp<0)=0;
    [Dr_temp_p,~] = perdecomp_3D(Dr_temp);
    Dk_temp = fftshift(fftn(fftshift(Dr_temp_p)));
    Dk_shift = Dk_temp;
    eval(['Dk_' num2str(i-1) '=Dk_shift;']);
end

%% make seperation
[sep_matrix]=make_forward_separation_matrix(5,5);
inv_sep_matrix=pinv(sep_matrix);

O=zeros(ny,nx,nz,n_phase);

for ii=0:n_phase-1
    %ii
    eval(['O(:,:,:,1)=O(:,:,:,1)+inv_sep_matrix(1,ii+1)*Dk_',num2str(ii),';']); %m2 order
    eval(['O(:,:,:,2)=O(:,:,:,2)+inv_sep_matrix(2,ii+1)*Dk_',num2str(ii),';']); %m1 order
    eval(['O(:,:,:,3)=O(:,:,:,3)+inv_sep_matrix(3,ii+1)*Dk_',num2str(ii),';']); %0th order
    eval(['O(:,:,:,4)=O(:,:,:,4)+inv_sep_matrix(4,ii+1)*Dk_',num2str(ii),';']); %p1 order
    eval(['O(:,:,:,5)=O(:,:,:,5)+inv_sep_matrix(5,ii+1)*Dk_',num2str(ii),';']); %p2 order
end


end



