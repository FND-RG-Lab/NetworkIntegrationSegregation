clear;clc;
addpath( 'scripts_matlab');

subjects_file_list='subjects.txt'; %txt file with subject IDs (1 per line)
mask_file='GM_mask_3mm.nii.gz'; %grey matter mask in the resolution of your data
atlas_file='Yeo7_3mm.nii.gz'; %network parcellation map
data_path='data3mm'; % a folder with processed fmri files (one per individual)
out_path='results/intgr_segr_3mm/'; %output folder
apply_fisher=0;
block_sz=10000;


%%%%%%%%%%%%%%% Code starts here %%%%%%%%%%%%%%%%%%

%Read subjects filenames from txt file
fid = fopen(subjects_file_list);
tline = fgetl(fid);
subjects={}; num_subjects=0;
while ischar(tline)
    num_subjects=num_subjects+1;
    subjects{num_subjects}=tline;
    tline = fgetl(fid);
end
fclose(fid);

%Read mask
[hdr mask]=read_nii(mask_file);

%Read atlas
[hdra atlas]=read_nii(atlas_file);
atlas=atlas.*(mask>0); %apply mask to atlas
atlas=reshape(atlas,[],1); %reshape to 1d
indx_mask=find(atlas>0); %find gm voxels
num_voxels=length(indx_mask);
atlas=atlas(indx_mask);
num_rois_atlas=max(atlas(:));

%Create output directories
if exist([out_path '/segregation'])==0
   mkdir([out_path '/segregation']);
end

if exist([out_path '/integration'])==0
   mkdir([out_path '/integration']);
end

%For each file
files=strsplit(ls([data_path '/*.nii.gz'])); files(end)=[];
seg_all = zeros(1,num_voxels,num_subjects);
int_all = zeros(1,num_voxels,num_subjects);
for i=1:num_subjects %length(files)

  %read fmri
  fprintf('[%.3d/%.3d] %s\n',i,num_subjects,subjects{i});
  [hdri fmri]=read_nii([ data_path '/' subjects{i} '.nii.gz']);
  fmri=reshape(fmri,[],size(fmri,4)); % convert fmri to 2d matrix
  fmri=fmri(indx_mask,:); %apply mask

  %Comptue integration, segregation
  segregation=zeros(1,num_voxels);
  integration=zeros(1,num_voxels);

  %for each roi
  for a=1:num_rois_atlas
      fprintf(' Network %d\n',a);
      sgr=find(atlas==a);
      intr=setdiff(1:num_voxels,sgr);
      num_blocks=ceil(length(sgr)./block_sz);
      for b=1:num_blocks
          sv=((b-1)*block_sz)+1;
          ev=sv+block_sz-1;
          if ev>length(sgr)
              ev=length(sgr);
          end

          fc=corr(fmri(sgr(sv:ev),:)',fmri');
          fc(find(fc<0))=0;
          if (apply_fisher==1)
            fc= 0.5 * log( ( 1 + fc)./( 1 - fc) ); %fisher
          end

          segregation(sgr(sv:ev))=nansum(fc(:,sgr)')-1;
          integration(sgr(sv:ev))=nansum(fc(:,intr)');
          clear fc;
      end
      clear sgr intr;
  end

  int_all(:,:,i) = integration;
  seg_all(:,:,i) = segregation;

  img=zeros(size(mask));
  img(indx_mask)=segregation;
  save_nii(hdr,img,[out_path '/segregation/' subjects{i} '.nii.gz']);

  img=zeros(size(mask));
  img(indx_mask)=integration;
  save_nii(hdr,img,[out_path '/integration/' subjects{i} '.nii.gz']);

end

save([out_path '/integration_voxelLevel.mat'],'int_all');
save([out_path '/segregation_voxelLevel.mat'],'seg_all');