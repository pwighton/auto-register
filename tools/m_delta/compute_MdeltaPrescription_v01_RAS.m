function compute_MdeltaPrescription_v01_RAS(pos, RO_dir, PE_dir, filestr)
% compute_MdeltaPrescription_RAS(pos, RO_dir, PE_dir, filestr)
% 
% Takes as inputs pos, RO_dir and PE_dir in RAS and converts to LPS to
% pass to compute_MdeltaPrescription_v01

RAS2LPS =  [-1  0  0
             0 -1  0
             0  0  1];

#fprintf('nargin is %d\n', nargin)

if nargin < 3
  error('Must pass at least pos, RO_dir and PE_dir')
elseif nargin < 4
  compute_MdeltaPrescription_v01(RAS2LPS*pos(:), RAS2LPS*RO_dir(:), RAS2LPS*PE_dir(:));
else
  compute_MdeltaPrescription_v01(RAS2LPS*pos(:), RAS2LPS*RO_dir(:), RAS2LPS*PE_dir(:), filestr);
end

return
