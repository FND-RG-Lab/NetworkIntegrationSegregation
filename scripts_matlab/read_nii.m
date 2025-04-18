function [hdr img]= read_nii(filename)
% Adaptation from https://es.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
% Suports Nifti (*.nii or *.nii.gz) and analyze data (*.hdr/*.img)
%
% TODO: xform and reorient to RPI (take a look http://cosmomvpa.org/matlab/cosmo_fmri_dataset.html
%                                              https://github.com/CoSMoMVPA/CoSMoMVPA/blob/master/mvpa/cosmo_fmri_reorient.m)

 isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
 if isOctave
    warning('off', 'Octave:possible-matlab-short-circuit-operator');
    confirm_recursive_rmdir(0);
 end

% Check if file exist
 if ~exist(filename,'file')
    error(['File doesn''t exist: ',filename]);
 end 

 %If compressed unzip the file
 [pathstr,name,ext] = fileparts(filename); 
 if strcmp(ext,'.gz')
	tmpDir = tempname;
    mkdir(tmpDir);
    if isOctave
	  copyfile(filename,tmpDir);
      filename = gunzip([tmpDir '/' name ext ], tmpDir);
    else
      filename = gunzip(filename, tmpDir);
    end
    filename = char(filename);
 elseif strcmp(ext,'.img')
    filename=[pathstr,'/',name,'.hdr'];
 end
 
 %Check machine and if is valid nifti image
 machine = 'ieee-le';
 fid = fopen(filename,'r',machine);
 if (fid < 0)
    error(sprintf('Cannot open file %s.',filename));
 else
    fseek(fid,0,'bof');
    testval=fread(fid,1,'int32');
	fclose(fid);
    if testval ~= 348
        switch machine,
           case 'ieee-le', machine = 'ieee-be';
           case 'ieee-be', machine = 'ieee-le';
        end

        fid = fopen(filename,'r',machine);
        if fid < 0,
            error(sprintf('Cannot open file %s.',filename));
        else
            fseek(fid,0,'bof');
            if fread(fid,1,'int32') ~= 348
               error(sprintf('File format is not valid',filename));
            end
            fclose(fid);
        end
     end
 end
 
 %Read header
 hdr = read_hdr(filename,machine);
 
 %Read image
 isnifti=1;
 if strcmp(ext,'.hdr')
    filename=[pathstr,'/',name,'.img'];
	isnifti=0;
 end
 
 [img,hdr] = read_image(hdr,filename,machine,isnifti);  
   
 
 %Perform some of sform/qform transform
 %[hdr,img] = xform_nii(hdr,img,filename);
 
 
 %Remove temp folder if gz file
 if exist('tmpDir', 'var')
   rmdir(tmpDir,'s');
 end
 
 return
 

 
 
function [img,hdr]=read_image(hdr,filename,machine,isnifti)
   fid = fopen(filename,'r',machine);

   if fid < 0,
      msg = sprintf('Cannot open file %s.',fn);
      error(msg);
   end

   %  Set bitpix according to datatype
   %
   %  /*Acceptable values for datatype are*/ 
   %
   %     0 None                     (Unknown bit per voxel) % DT_NONE, DT_UNKNOWN 
   %     1 Binary                         (ubit1, bitpix=1) % DT_BINARY 
   %     2 Unsigned char         (uchar or uint8, bitpix=8) % DT_UINT8, NIFTI_TYPE_UINT8 
   %     4 Signed short                  (int16, bitpix=16) % DT_INT16, NIFTI_TYPE_INT16 
   %     8 Signed integer                (int32, bitpix=32) % DT_INT32, NIFTI_TYPE_INT32 
   %    16 Floating point    (single or float32, bitpix=32) % DT_FLOAT32, NIFTI_TYPE_FLOAT32 
   %    32 Complex, 2 float32      (Use float32, bitpix=64) % DT_COMPLEX64, NIFTI_TYPE_COMPLEX64
   %    64 Double precision  (double or float64, bitpix=64) % DT_FLOAT64, NIFTI_TYPE_FLOAT64 
   %   128 uint8 RGB                 (Use uint8, bitpix=24) % DT_RGB24, NIFTI_TYPE_RGB24 
   %   256 Signed char            (schar or int8, bitpix=8) % DT_INT8, NIFTI_TYPE_INT8 
   %   511 Single RGB              (Use float32, bitpix=96) % DT_RGB96, NIFTI_TYPE_RGB96
   %   512 Unsigned short               (uint16, bitpix=16) % DT_UNINT16, NIFTI_TYPE_UNINT16 
   %   768 Unsigned integer             (uint32, bitpix=32) % DT_UNINT32, NIFTI_TYPE_UNINT32 
   %  1024 Signed long long              (int64, bitpix=64) % DT_INT64, NIFTI_TYPE_INT64
   %  1280 Unsigned long long           (uint64, bitpix=64) % DT_UINT64, NIFTI_TYPE_UINT64 
   %  1536 Long double, float128  (Unsupported, bitpix=128) % DT_FLOAT128, NIFTI_TYPE_FLOAT128 
   %  1792 Complex128, 2 float64  (Use float64, bitpix=128) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
   %  2048 Complex256, 2 float128 (Unsupported, bitpix=256) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
   %
   
   switch hdr.image_dimension.datatype
   case   1,
      hdr.image_dimension.bitpix = 1;  precision = 'ubit1';
   case   2,
      hdr.image_dimension.bitpix = 8;  precision = 'uint8';
   case   4,
      hdr.image_dimension.bitpix = 16; precision = 'int16';
   case   8,
      hdr.image_dimension.bitpix = 32; precision = 'int32';
   case  16,
      hdr.image_dimension.bitpix = 32; precision = 'float32';
   case  32,
      hdr.image_dimension.bitpix = 64; precision = 'float32';
   case  64,
      hdr.image_dimension.bitpix = 64; precision = 'float64';
   case 128,
      hdr.image_dimension.bitpix = 24; precision = 'uint8';
   case 256 
      hdr.image_dimension.bitpix = 8;  precision = 'int8';
   case 511 
      hdr.image_dimension.bitpix = 96; precision = 'float32';
   case 512 
      hdr.image_dimension.bitpix = 16; precision = 'uint16';
   case 768 
      hdr.image_dimension.bitpix = 32; precision = 'uint32';
   case 1024
      hdr.image_dimension.bitpix = 64; precision = 'int64';
   case 1280
      hdr.image_dimension.bitpix = 64; precision = 'uint64';
   case 1792,
      hdr.image_dimension.bitpix = 128; precision = 'float64';
   otherwise
      error('This datatype is not supported'); 
   end
   hdr.image_dimension.dim(find(hdr.image_dimension.dim < 1)) = 1;
   

   % move pointer to the start of image block
   if isnifti
      fseek(fid, hdr.image_dimension.vox_offset, 'bof');
   else
      fseek(fid, 0, 'bof');
   end


   % For each frame, precision of value will be read in img_siz times, where img_siz is only the dimension size of an image, not the byte storage size of an image.
      
   img_siz = prod(hdr.image_dimension.dim(2:8));

   %  For complex float32 or complex float64, voxel values include [real, imag]
   if hdr.image_dimension.datatype == 32 | hdr.image_dimension.datatype == 1792
      img_siz = img_siz * 2;
   end
	 
   % For RGB24, voxel values include 3 separate color planes
   if hdr.image_dimension.datatype == 128 | hdr.image_dimension.datatype == 511
  	  img_siz = img_siz * 3;
   end

   img = fread(fid, img_siz, sprintf('*%s',precision));

   d1 = hdr.image_dimension.dim(2);
   d2 = hdr.image_dimension.dim(3);
   d3 = hdr.image_dimension.dim(4);
   d4 = hdr.image_dimension.dim(5);
   d5 = hdr.image_dimension.dim(6);
   d6 = hdr.image_dimension.dim(7);
   d7 = hdr.image_dimension.dim(8);

   img_idx = 1:d4; dim5_idx = 1:d5; dim6_idx = 1:d6; dim7_idx = 1:d7;

 
   %  For complex float32 or complex float64, voxel values 
   if hdr.image_dimension.datatype == 32 | hdr.image_dimension.datatype == 1792
      img = reshape(img, [2, length(img)/2]);
      img = complex(img(1,:)', img(2,:)');
   end

   fclose(fid);

   %  Update the global min and max values 
   
   hdr.image_dimension.glmax = double(max(img(:)));
   hdr.image_dimension.glmin = double(min(img(:)));

   
   % For RGB data
   if hdr.image_dimension.datatype == 128 & hdr.image_dimension.bitpix == 24
      % remove squeeze
      img = (reshape(img, [3 hdr.image_dimension.dim(2:4) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
      img = permute(img, [2 3 4 1 5 6 7 8]);
   elseif hdr.image_dimension.datatype == 511 & hdr.image_dimension.bitpix == 96
      img = double(img(:));
      img = single((img - min(img))/(max(img) - min(img)));
      % remove squeeze
      img = (reshape(img, [3 hdr.image_dimension.dim(2:4) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
      img = permute(img, [2 3 4 1 5 6 7 8]);
   else
      % remove squeeze
      img = (reshape(img, [hdr.image_dimension.dim(2:4) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
   end

   hdr.image_dimension.dim(5) = length(img_idx);
   hdr.image_dimension.dim(6) = length(dim5_idx);
   hdr.image_dimension.dim(7) = length(dim6_idx);
   hdr.image_dimension.dim(8) = length(dim7_idx);
  

 return 
 
 
function [hdr img] = xform_nii(hdr,img,filename)
   tolerance=0.1; %  10%
   
   % if scl_slope field is nonzero, then each voxel value in the dataset should be scaled as: y = scl_slope * x + scl_inter
   if hdr.image_dimension.scl_slope ~= 0 & ismember(hdr.image_dimension.datatype, [2,4,8,16,64,256,512,768]) & (hdr.image_dimension.scl_slope ~= 1 | hdr.image_dimension.scl_inter ~= 0)
      img = hdr.image_dimension.scl_slope * double(img) + hdr.image_dimension.scl_inter;
      if hdr.image_dimension.datatype == 64
         hdr.image_dimension.bitpix = 64;
      else
         img = single(img);
         hdr.image_dimension.datatype = 16;
         hdr.image_dimension.bitpix = 32;
      end
      hdr.image_dimension.glmax = max(double(img(:)));
      hdr.image_dimension.glmin = min(double(img(:)));
      %  set scale to non-use
      hdr.image_dimension.scl_slope = 0;
   end

   % the scaling is ignored if datatype is DT_RGB24.

   %  If datatype is a complex type, then the scaling is to be applied to both the real and imaginary parts.
  
   if hdr.image_dimension.scl_slope ~= 0 & ismember(hdr.image_dimension.datatype, [32,1792])
      img = hdr.image_dimension.scl_slope * double(img) + hdr.image_dimension.scl_inter;
      if hdr.image_dimension.datatype == 32
         img = single(img);
      end
      hdr.image_dimension.glmax = max(double(img(:)));
      hdr.image_dimension.glmin = min(double(img(:)));
      %  set scale to non-use
      hdr.image_dimension.scl_slope = 0;

    end

	
   %  There is no need for this program to transform Analyze data
   if (~strcmp(hdr.data_history.magic, 'n+1')) & (~strcmp(hdr.data_history.magic, 'ni1'))
       [pathstr,name,ext] = fileparts(filename); 
	   if exist([pathstr,'/',name,'.mat'],'file')
	      load([pathstr,'/',name,'.mat']);	% old SPM affine matrix
          R=M(1:3,1:3);
          T=M(1:3,4);
          T=R*ones(3,1)+T;
          M(1:3,4)=T;
          hdr.data_history.qform_code=0;
          hdr.data_history.sform_code=1;
          hdr.data_history.srow_x=M(1,:);
          hdr.data_history.srow_y=M(2,:);
          hdr.data_history.srow_z=M(3,:);
	   else
	      hdr.data_history.rot_orient = [];
          hdr.data_history.flip_orient = [];
		  return;
	   end
   end
   

   [hdr,orient]=change_hdr(hdr,tolerance);

   %  flip and/or rotate image data
   if ~isequal(orient, [1 2 3])

      old_dim = hdr.image_dimension.dim([2:4]);

      %  More than 1 time frame
      if ndims(img) > 3
         pattern = 1:prod(old_dim);
		 pattern = reshape(pattern, old_dim);
      else
         pattern = [];
      end


      % calculate for rotation after flip
      rot_orient = mod(orient + 2, 3) + 1;

      % do flip:
	  flip_orient = orient - rot_orient;

      for i = 1:3
         if flip_orient(i)
            if ~isempty(pattern)
               pattern = flipdim(pattern, i);
            else
               img = flipdim(img, i);
            end
         end
      end

      %  get index of orient (rotate inversely)
      [tmp rot_orient] = sort(rot_orient);
      new_dim = old_dim;
      new_dim = new_dim(rot_orient);
      hdr.image_dimension.dim([2:4]) = new_dim;

      new_pixdim = hdr.image_dimension.pixdim([2:4]);
      new_pixdim = new_pixdim(rot_orient);
      hdr.image_dimension.pixdim([2:4]) = new_pixdim;

      % re-calculate originator
      tmp = hdr.data_history.originator([1:3]);
      tmp = tmp(rot_orient);
      flip_orient = flip_orient(rot_orient);
      for i = 1:3
         if flip_orient(i) & ~isequal(tmp(i), 0)
            tmp(i) = new_dim(i) - tmp(i) + 1;
         end
      end

      hdr.data_history.originator([1:3]) = tmp;
      hdr.data_history.rot_orient = rot_orient;
      hdr.data_history.flip_orient = flip_orient;

      %  do rotation
      if ~isempty(pattern)
         pattern = permute(pattern, rot_orient);
         pattern = pattern(:);

        if hdr.image_dimension.datatype == 32 | hdr.image_dimension.datatype == 1792 | hdr.image_dimension.datatype == 128 | hdr.image_dimension.datatype == 511
            tmp = reshape(img(:,:,:,1), [prod(new_dim) hdr.image_dimension.dim(5:8)]);
            tmp = tmp(pattern, :);
            img(:,:,:,1) = reshape(tmp, [new_dim       hdr.image_dimension.dim(5:8)]);

            tmp = reshape(img(:,:,:,2), [prod(new_dim) hdr.image_dimension.dim(5:8)]);
            tmp = tmp(pattern, :);
            img(:,:,:,2) = reshape(tmp, [new_dim       hdr.image_dimension.dim(5:8)]);

            if hdr.image_dimension.datatype == 128 | hdr.image_dimension.datatype == 511
               tmp = reshape(img(:,:,:,3), [prod(new_dim) hdr.image_dimension.dim(5:8)]);
               tmp = tmp(pattern, :);
               img(:,:,:,3) = reshape(tmp, [new_dim       hdr.image_dimension.dim(5:8)]);
            end

         else
            img = reshape(img, [prod(new_dim) hdr.image_dimension.dim(5:8)]);
            img = img(pattern, :);
            img = reshape(img, [new_dim       hdr.image_dimension.dim(5:8)]);
         end
      else
         if hdr.image_dimension.datatype == 32 | hdr.image_dimension.datatype == 1792 | hdr.image_dimension.datatype == 128 | hdr.image_dimension.datatype == 511
            img(:,:,:,1) = permute(img(:,:,:,1), rot_orient);
            img(:,:,:,2) = permute(img(:,:,:,2), rot_orient);

            if hdr.image_dimension.datatype == 128 | hdr.image_dimension.datatype == 511
               img(:,:,:,3) = permute(img(:,:,:,3), rot_orient);
            end
         else
            img = permute(img, rot_orient);
         end
      end
   else
      hdr.data_history.rot_orient = [];
      hdr.data_history.flip_orient = [];
   end

   return


function [hdr, orient] = change_hdr(hdr, tolerance)

   orient = [1 2 3];
   affine_transform = 1;

   %  NIFTI can have both sform and qform transform. This program will check sform_code prior to qform_code by default.

   useForm='s';
   if hdr.data_history.sform_code > 0
	   useForm='s';
   elseif hdr.data_history.qform_code > 0
	   useForm='q';
   end					
   

   if isequal(useForm,'s')
      R = [hdr.data_history.srow_x(1:3)
           hdr.data_history.srow_y(1:3)
           hdr.data_history.srow_z(1:3)];

      T = [hdr.data_history.srow_x(4)
           hdr.data_history.srow_y(4)
           hdr.data_history.srow_z(4)];

      if det(R) == 0 | ~isequal(R(find(R)), sum(R)')
         old_affine = [ [R;[0 0 0]] [T;1] ];
         R_sort = sort(abs(R(:)));
         R( find( abs(R) < tolerance*min(R_sort(end-2:end)) ) ) = 0;
         new_affine = [ [R;[0 0 0]] [T;1] ];

         if det(R) == 0 | ~isequal(R(find(R)), sum(R)')
            error('   Non-orthogonal rotation or shearing ');
         end
      end

   elseif isequal(useForm,'q')
      b = hdr.data_history.quatern_b;
      c = hdr.data_history.quatern_c;
      d = hdr.data_history.quatern_d;

      if 1.0-(b*b+c*c+d*d) < 0
         if abs(1.0-(b*b+c*c+d*d)) < 1e-5
            a = 0;
         else
            error('Incorrect quaternion values in this NIFTI data.');
         end
      else
         a = sqrt(1.0-(b*b+c*c+d*d));
      end

      qfac = hdr.image_dimension.pixdim(1);
      if qfac==0, qfac = 1; end
      i = hdr.image_dimension.pixdim(2);
      j = hdr.image_dimension.pixdim(3);
      k = qfac * hdr.image_dimension.pixdim(4);

      R = [a*a+b*b-c*c-d*d     2*b*c-2*a*d        2*b*d+2*a*c
           2*b*c+2*a*d         a*a+c*c-b*b-d*d    2*c*d-2*a*b
           2*b*d-2*a*c         2*c*d+2*a*b        a*a+d*d-c*c-b*b];

      T = [hdr.data_history.qoffset_x
           hdr.data_history.qoffset_y
           hdr.data_history.qoffset_z];

      %  qforms are expected to generate rotation matrices R which are
      %  det(R) = 1; we'll make sure that happens.
      %  
      %  now we make the same checks as were done above for sform data
      %  BUT we do it on a transform that is in terms of voxels not mm;
      %  after we figure out the angles and squash them to closest 
      %  rectilinear direction. After that, the voxel sizes are then
      %  added.
      %
      %  This part is modified by Jeff Gunter.
      %
      if det(R) == 0 | ~isequal(R(find(R)), sum(R)')

         %  det(R) == 0 is not a common trigger for this ---
         %  R(find(R)) is a list of non-zero elements in R; if that
         %  is straight (not oblique) then it should be the same as 
         %  columnwise summation. Could just as well have checked the
         %  lengths of R(find(R)) and sum(R)' (which should be 3)
         %
         old_affine = [ [R * diag([i j k]);[0 0 0]] [T;1] ];
         R_sort = sort(abs(R(:)));
         R( find( abs(R) < tolerance*min(R_sort(end-2:end)) ) ) = 0;
         R = R * diag([i j k]);
         new_affine = [ [R;[0 0 0]] [T;1] ];

         if det(R) == 0 | ~isequal(R(find(R)), sum(R)')
            error('   Non-orthogonal rotation or shearing ');
         end

      else
         R = R * diag([i j k]);
      end					% 1st det(R)

   else
      affine_transform = 0;	% no sform or qform transform
   end

   if affine_transform == 1
      voxel_size = abs(sum(R,1));
      inv_R = inv(R);
      originator = inv_R*(-T)+1;
      orient = get_orient(inv_R);  
	  
      %  modify pixdim and originator
      %
      hdr.image_dimension.pixdim(2:4) = voxel_size;
      hdr.data_history.originator(1:3) = originator;

      %  set sform or qform to non-use, because they have been
      %  applied in xform_nii
      %
      hdr.data_history.qform_code = 0;
      hdr.data_history.sform_code = 0;
   end


   %  apply space_unit to pixdim if not 1 (mm)
   
   switch bitand(hdr.image_dimension.xyzt_units, 7)	% mask with 0x07
   case 1
      space_unit = 1e+3;		% meter, m
   case 3
      space_unit = 1e-3;		% micrometer, um
   otherwise
      space_unit = 1;			% millimeter, mm
   end

   if space_unit ~= 1
      hdr.image_dimension.pixdim(2:4) = hdr.image_dimension.pixdim(2:4) * space_unit;
      %  set space_unit of xyzt_units to millimeter, because  voxel_size has been re-scaled
      hdr.image_dimension.xyzt_units = char(bitset(hdr.image_dimension.xyzt_units,1,0));
      hdr.image_dimension.xyzt_units = char(bitset(hdr.image_dimension.xyzt_units,2,1));
      hdr.image_dimension.xyzt_units = char(bitset(hdr.image_dimension.xyzt_units,3,0));
   end

   hdr.image_dimension.pixdim = abs(hdr.image_dimension.pixdim);

   return;					
 

 
function orient = get_orient(R)

   orient = [];

   for i = 1:3
      switch find(R(i,:)) * sign(sum(R(i,:)))
      case 1
         orient = [orient 1];		% Left to Right
      case 2
         orient = [orient 2];		% Posterior to Anterior
      case 3
         orient = [orient 3];		% Inferior to Superior
      case -1
         orient = [orient 4];		% Right to Left
      case -2
         orient = [orient 5];		% Anterior to Posterior
      case -3
         orient = [orient 6];		% Superior to Inferior
      end
   end

   return;					% get_orient
 
 
 
