clear;clc;

addpath('scripts_matlab');

subjects_file_list='subjects.txt'; %txt file with subject IDs (one per line)
mask_file='GM_mask_3mm.nii.gz'; %grey matter mask in the resolution of your data
data_path='data3mm'; % a folder with processed fmri files (one per individual)
out_path='results/weighted_degree_voxel'; %output folder for the results
apply_fisher=0;


%%%%%%%%%%%%%%%%%%%%%% Code starts here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Read subjects ids from the text file
fid = fopen(subjects_file_list);
tline = fgetl(fid);
subjects={}; num_subjects=0;
while ischar(tline)
    num_subjects=num_subjects+1;
    subjects{num_subjects}=tline;
    tline = fgetl(fid);
end
fclose(fid);

%Read mask and identify the voxels of gm
[hdr mask]=read_nii(mask_file);
indx_mask=find(reshape(mask,[],1)>0);
num_rois=length(indx_mask);

%For each subject compute the weighted degree
mkdir(out_path); 
wd_all = zeros(1,num_rois,num_subjects);
for s=1:num_subjects
   fprintf('[%.3d/%.3d] %s\n',s,num_subjects,subjects{s});

   %Read the fmri time series convert to 2d and apply grey matter mask
   [hdri img]=read_nii([ data_path '/' subjects{s} '.nii.gz']);
   img=reshape(img,[],size(img,4));
   img=img(indx_mask,:);

  %Compute functional connectivity matrix
  fc=corr(img');
  fc(1:num_rois+1:end)=0; %diagonal to zero
  if (apply_fisher==1)
      fc=0.5*log((1+fc)./(1-fc));
  end

  clear img fcs;

  %Compute weighted degree
  wd=nansum(fc);

  %Save weighted degree map
  img2=zeros(size(mask));
  img2(indx_mask)=wd;
  save_nii(hdr,img2,[ out_path  '/' subjects{s} '.nii.gz' ]);

  wd_all(:,:,s) = wd;
  clear fc wd img2;
  
end

save([out_path '/wd_voxelLevel.mat'],'wd_all');