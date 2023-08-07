function [M0_vox2ras, vox_mm, slice_gap_pct, R] = mrir_measdat_vox2ras(prot, evp, varargin)
%MRIR_MEASDAT_VOX2RAS
%
% M0_vox2ras = mrir_measdat_vox2ras(prot, evp)
% M0_vox2ras = mrir_measdat_vox2ras(prot, evp, [is3D])
% M0_vox2ras = mrir_measdat_vox2ras(prot, evp, ..., [dims])
%
% [M0_vox2ras, vox_mm, slice_gap_pct] = mrir_measdat_vox2ras(prot, evp)
%
% see also MRIR_MEASDAT_DIRCOS.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2012/nov/26
% $Id: mrir_measdat_vox2ras.m,v 1.2 2014/03/03 16:00:56 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  FLAG__is3D = 0;
  if ( nargin >= 3 ),
    FLAG__is3D = varargin{1};
  end;

  if ( evp.NParMeas > 1 ),
    FLAG__is3D = 1;
  end;

  dims = [];
  if ( nargin >= 4 ),
    dims = varargin{2};
    dims(end+1:16) = 1;
  end;


  OPTION__mimic_FreeSurfer      = 1;
  OPTION__force_DICOM_standard  = 0;
  OPTION__flip_phaseencode      = 0;


  %==--------------------------------------------------------------------==%

  [COL_dircos_LPS, LIN_dircos_LPS, IMG_dircos_LPS] = mrir_measdat_vox2ras__dircos(prot);

  % [2014/jul/12: for kuleuven data]
  if ( OPTION__flip_phaseencode ),
    LIN_dircos_LPS = (-1) * LIN_dircos_LPS;
  end;


  if ( ~strcmp(prot.tPatientPosition, 'HFS') ),
    error('not HFS!')
  end;

  % TODO: clarify that siemens orientation is LPS only when patient is HFS
  R_LPS2RAS = [-1 0 0; 0 -1 0; 0 0 +1];

  COL_dircos_RAS = R_LPS2RAS * COL_dircos_LPS;
  LIN_dircos_RAS = R_LPS2RAS * LIN_dircos_LPS;
  IMG_dircos_RAS = R_LPS2RAS * IMG_dircos_LPS;

  % temporary hack for debugging
  %COL_dircos_RAS = -COL_dircos_RAS;
  %LIN_dircos_RAS = -LIN_dircos_RAS;
  %IMG_dircos_RAS = -IMG_dircos_RAS;

  COL_vox_mm = prot.sSliceArray(1).dReadoutFOV/evp.NImageCols;
  LIN_vox_mm = prot.sSliceArray(1).dPhaseFOV/evp.NImageLins;

%  if ( FLAG__no_resample ),
%  COL_vox_mm = prot.sSliceArray(1).dReadoutFOV/evp.NColMeas*prot.flReadoutOSFactor
%  LIN_vox_mm = prot.sSliceArray(1).dPhaseFOV/prot.lPhaseEncodingLines
  IMG_vox_mm = prot.sSliceArray(1).dThickness;

  distance_between_slices = 0;
  slice_gap_pct = 0;
  if ( prot.sSliceArray_lSize > 1 && length(prot.sSliceArray) > 1 ),
    % check for slice gap

    pos_slice1 = [
        prot.sSliceArray(1).sPosition_dSag;
        prot.sSliceArray(1).sPosition_dCor;
        prot.sSliceArray(1).sPosition_dTra;
                 ];
    pos_slice2 = [
        prot.sSliceArray(2).sPosition_dSag;
        prot.sSliceArray(2).sPosition_dCor;
        prot.sSliceArray(2).sPosition_dTra;
                 ];

    distance_between_slices = norm(pos_slice2 - pos_slice1);

    slice_gap_pct = (distance_between_slices-IMG_vox_mm)/IMG_vox_mm;
    slice_gap_pct = round(slice_gap_pct * 1000)/1000;
    if ( slice_gap_pct == 0 ), slice_gap_pct = abs(slice_gap_pct); end;

    % this seems like a bad idea, but it appears to be the convention---so
    % most headers don't distinguish between 1 mm slices with a 25% gap and
    % 1.25 mm slices.
    IMG_vox_mm = distance_between_slices;

    disp(sprintf('==> [%s]: slice gap = %2.2f', mfilename, slice_gap_pct));

  end;

  if ( isempty(evp.NImagePar) ),
    evp.NImagePar = 1;
  end;

  % in 3D volumes the thickness is given for the entire slab, so this has to
  % be corrected here
  if ( FLAG__is3D ),

    % if 3D and prot.sSliceArray_lSize > 1, this means multi-slab, in which
    % case the IMG_vox_mm is not the size in the partition direction. even
    % with negative slice gap, the slab thickness is accurate, so voxel size
    % is slab thickness divided by NImagePar. once negative overlap is
    % accounted for, the total number of slices will be less than
    % evp.NImagePar*prot.sSliceArray_lSize by a factor of:
    %    slice_gap_pct*(prot.sSliceArray_lSize-1)/prot.sSliceArray_lSize

    %    IMG_vox_mm = IMG_vox_mm/prot.lPartitions;

    IMG_vox_mm = prot.sSliceArray(1).dThickness/evp.NImagePar;
  end;


  % rotation component of matrix, with axis of rotation at centroid of image
  % volume
  R = [ [COL_dircos_RAS*COL_vox_mm, LIN_dircos_RAS*LIN_vox_mm, IMG_dircos_RAS*IMG_vox_mm, [0; 0; 0]]
        [0,0,0,1] ];

%  if ( FLAG__no_resample ),
%  COL_vox_count = prot.lBaseResolution;
%  LIN_vox_count = prot.lPhaseEncodingLines;
%  IMG_vox_count = max([prot.lPartitions, prot.sSliceArray_lSize]);
  COL_vox_count = evp.NImageCols;
  LIN_vox_count = evp.NImageLins;
  IMG_vox_count = max([evp.NImagePar, evp.NSlcMeas]);

  % override evp params
  if ( ~isempty(dims) ),
      COL_vox_count = dims(1);
      LIN_vox_count = dims(2);
      IMG_vox_count = max([dims(9), dims(10)]);
  end;

  % index of voxel at center of volume (assume integer since all dims divisible by 2)
  P_center_VOX = [COL_vox_count/2; LIN_vox_count/2; ceil(IMG_vox_count/2); 1];


  if ( FLAG__is3D ),

    % here we subtract because, if the corner voxel is the origin the
    % displacement from its centroid to the position in the center of the
    % volume is not an integer number of voxels but 1/2 voxel less, i.e.,
    % the volume center is *closer* ot the corner voxel centroid than it is
    % to the corner voxel edge.
    if ( OPTION__mimic_FreeSurfer ),
      P_center_VOX = P_center_VOX - [0.0; 0.0; 0.5; 0];
    elseif ( OPTION__force_DICOM_standard ),
      P_center_VOX = P_center_VOX - [0.5; 0.5; 0.5; 0];
    end;
  end;

  % HACK #1: for 3D, out of laziness we can assume that the first slice is
  % the center
  if ( ~FLAG__is3D ),
    P_center_VOX(3) = 0;
  end;


  % displacement vector between corner voxel and volume centroid (if the
  % XYZ_mm origin were at the first voxel, this would be the XYZ_mm
  % coordinate of the volume centroid)
  P_center_RAS = R * P_center_VOX;

  % the position of the voxel at the "center" of the volume
  Q_center_LPS = [prot.sSliceArray(1).sPosition_dSag; ...
                  prot.sSliceArray(1).sPosition_dCor; ...
                  prot.sSliceArray(1).sPosition_dTra];
  % if we don't use HACK #1, this line would be something like:
  %
  %  Q_center_LPS = [prot.sSliceArray(ceil(IMG_vox_count/2)+1).sPosition_dSag; ...
  %                  prot.sSliceArray(ceil(IMG_vox_count/2)+1).sPosition_dCor; ...
  %                  prot.sSliceArray(ceil(IMG_vox_count/2)+1).sPosition_dTra];
  %
  % does this work for cases with an odd numbers of slices?

  % the most important thing is for P_center_VOX and Q_center_LPS to
  % represent the same point in the volume!

  % (TODO: DICOM standard requires the position of the "first voxel
  % transmitted", which it claims is the "upper left hand corner of the
  % image". because the final reconstructed images are reoriented
  % consistently to be in radiological orientation (regardless of the
  % polarity or direction of the RO or PE), check to see whether this
  % position changes in the meas.dat header when the polarity or direction
  % of RO or PE are changed.

  %
  Q_center_RAS = R_LPS2RAS * Q_center_LPS;


  % corner of volume RAS XYZ_mm coordinates that accounts for the shift of
  % the volume such that the corner lands at the location specified in the
  % header. (a.k.a. "ImgPos")
  D_corner_RAS = [Q_center_RAS; 1] - P_center_RAS;

  m = [COL_dircos_RAS, LIN_dircos_RAS, IMG_dircos_RAS];
  vox_mm = [COL_vox_mm; LIN_vox_mm; IMG_vox_mm];
  vox_count = [COL_vox_count; LIN_vox_count; IMG_vox_count];
  FOV_mm = vox_mm .* vox_count;
  %disp(sprintf(' c_ras:  %s', mat2str([D_corner_RAS(1:3) + (m * FOV_mm/2)])));
  %disp(sprintf(' ImgPos: %s', mat2str(D_corner_RAS(1:3))));


  M0_vox2ras = [[ R(1:3,1:3), D_corner_RAS(1:3) ]; [0,0,0,1]];

  % ON HALF-VOXEL SHIFTS:

  % meas.dat: position is corner of voxel

  % siemens DICOM: position is corner of plane and centroid of frame

  % DICOM standard: position is centroid of voxel

  % NOTE DICOM standard re-wrote this section because of ambiguity of
  % position in slice direction and noted potential for half-voxel shift
  % through-slice


  return;


%**************************************************************************%
function [Rdc, Pdc, Sdc] = mrir_measdat_vox2ras__dircos(prot)
%MRIR_MEASDAT_VOX2RAS__DIRCOS
%
% varargout = mrir_measdat_vox2ras(varargin)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2012/nov/17
% $Id: mrir_measdat_vox2ras.m,v 1.2 2014/03/03 16:00:56 jonp Exp $
%**************************************************************************%

  Sdc = [prot.sSliceArray(1).sNormal_dSag;
         prot.sSliceArray(1).sNormal_dCor;
         prot.sSliceArray(1).sNormal_dTra];

  inplanerot = prot.sSliceArray(1).dInPlaneRot;


  %==--------------------------------------------------------------------==%

  orientation = mrir_measdat_vox2ras__principalorientation(Sdc(3),Sdc(1),Sdc(2));

  % in-plane rotation is rotation about the slice-normal direction, so
  % ignoring that for now, there can be rotation either about the PE
  % direction (which by definition does not change the PE vector) or about
  % the RO direction, which we consider first. since PE is orthogonal to RO,
  % the PE direction can be found by

  switch(orientation),

   case 'SMS_TRANSVERSE', % PE dir is COR (A>>P), RO dir is SAG
    Pdc = cross(Sdc/norm(Sdc), [+1;0;0]);
   case 'SMS_CORONAL',    % PE dir is SAG (R>>L), RO dir is TRA
    Pdc = cross(Sdc/norm(Sdc), [0;0;+1]);
   case 'SMS_SAGITTAL',   % PE dir is COR (A>>P), RO dir is TRA
    Pdc = cross(Sdc/norm(Sdc), [0;0;-1]);
  end;

  % direction cosines assumed to be unit normal, and ||axb|| = ||a|| ||b|| sin(theta)
  Pdc = Pdc/norm(Pdc);

  % calculate the readout direction as GR = GS x GP
  Rdc = cross(Sdc, Pdc);

  % direction cosines assumed to be unit normal
  Rdc = Rdc/norm(Rdc);

  % apply the inplane rotation, if its not zero
  if ( inplanerot ~= 0.0 ),
    Pdc(1) = cos(inplanerot)*Pdc(1) - sin(inplanerot)*Rdc(1);
    Pdc(2) = cos(inplanerot)*Pdc(2) - sin(inplanerot)*Rdc(2);
    Pdc(3) = cos(inplanerot)*Pdc(3) - sin(inplanerot)*Rdc(3);

    % recalculate the readout direction
    Rdc = cross(Sdc, Pdc);

  end;

  return;


%**************************************************************************%
function rval = mrir_measdat_vox2ras__principalorientation(ntra, nsag, ncor)
%MRIR_MEASDAT_VOX2RAS__PRINCIPALORIENTATION
%
% varargout = mrir_measdat_vox2ras(varargin)
%
%
% see also MRIR_MEASDAT_VOX2RAS_ALT.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2012/nov/17
% $Id: mrir_measdat_vox2ras.m,v 1.2 2014/03/03 16:00:56 jonp Exp $
%**************************************************************************%

  abs_sag = abs(nsag);
  abs_cor = abs(ncor);
  abs_tra = abs(ntra);

  almost_equal_sc = isAlmostEqual(abs_sag, abs_cor);
  almost_equal_st = isAlmostEqual(abs_sag, abs_tra);
  almost_equal_tc = isAlmostEqual(abs_tra, abs_cor);

  % mainly transverse case, which captures most special cases
  if (   (almost_equal_sc && almost_equal_st) ...
         || (almost_equal_sc && (abs_sag < abs_tra)) ...
         || (almost_equal_st && (abs_sag > abs_cor)) ...
         || (almost_equal_tc && (abs_cor > abs_sag)) ...
         || ((abs_sag > abs_cor) && (abs_sag < abs_tra)) ...
         || ((abs_sag < abs_cor) && (abs_cor < abs_tra)) ...
         || ((abs_sag < abs_tra) && (abs_tra > abs_cor)) ...
         || ((abs_cor < abs_tra) && (abs_tra > abs_sag))),
    rval = 'SMS_TRANSVERSE';
  elseif ( (almost_equal_sc && (abs_sag > abs_tra)) ...
           || (almost_equal_st && (abs_sag < abs_cor)) ...
           || ((abs_sag < abs_cor) && (abs_cor > abs_tra)) ...
           || ((abs_sag > abs_tra) && (abs_sag < abs_cor)) ...
           || ((abs_sag < abs_tra) && (abs_tra < abs_cor)))
    rval = 'SMS_CORONAL';
  elseif ( (almost_equal_tc && (abs_cor < abs_sag)) ...
           || ((abs_sag > abs_cor) && (abs_sag > abs_tra)) ...
           || ((abs_cor > abs_tra) && (abs_cor < abs_sag)) ...
           || ((abs_cor < abs_tra) && (abs_tra < abs_sag)))
    rval = 'SMS_SAGITTAL';
  else
    error('invalid slice orientation encountered');
  end;

  return;


%**************************************************************************%
function bool = isAlmostEqual(val_a, val_b)
%ISALMOSTEQUAL
%
% varargout = isAlmostEqual(varargin)
%
%
% see also ISALMOSTEQUAL_ALT.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2012/nov/17
% $Id: mrir_measdat_vox2ras.m,v 1.2 2014/03/03 16:00:56 jonp Exp $
%**************************************************************************%

  diff = val_a - val_b;

  if ( (diff >= -1.0e-6) && (diff <= 1.0e-6) ),
    bool = 1;
  else,
    bool = 0;
  end;

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_measdat_vox2ras.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
