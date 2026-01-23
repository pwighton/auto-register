function compute_MdeltaPrescription_fromSS_v01_RAS(pos, SS_dir, filestr)
% compute_MdeltaPrescription_fromSS_v01_RAS(pos, SS_dir, filestr)
% 
% Given a slice-select direction (SS_dir), computes orthogonal RO_dir and PE_dir
% vectors, then calls compute_MdeltaPrescription_v01_RAS
%
% Inputs:
%   pos     - position in RAS coordinates [3x1]
%   SS_dir  - slice-select direction in RAS coordinates [3x1]
%   filestr - (optional) output filename

% Normalize SS_dir
SS_dir = SS_dir(:);
SS_dir = SS_dir / norm(SS_dir);

% Find a vector that's not parallel to SS_dir
% We'll try [1,0,0], [0,1,0], [0,0,1] and pick the one least parallel
candidates = [1 0 0; 0 1 0; 0 0 1]';
dots = abs(candidates' * SS_dir);
[~, idx] = min(dots);
v = candidates(:, idx);

% Compute RO_dir as cross product, then normalize
RO_dir = cross(SS_dir, v);
RO_dir = RO_dir / norm(RO_dir);

% Compute PE_dir as cross product to ensure orthogonality
PE_dir = cross(SS_dir, RO_dir);
PE_dir = PE_dir / norm(PE_dir);

% Verify orthogonality (optional, for debugging)
fprintf('Checking orthogonality:\n');
fprintf('  dot(RO_dir, PE_dir) = %e\n', dot(RO_dir, PE_dir));
fprintf('  dot(RO_dir, SS_dir) = %e\n', dot(RO_dir, SS_dir));
fprintf('  dot(PE_dir, SS_dir) = %e\n', dot(PE_dir, SS_dir));

% Call the main function
if nargin < 3
  compute_MdeltaPrescription_v01_RAS(pos, RO_dir, PE_dir);
else
  compute_MdeltaPrescription_v01_RAS(pos, RO_dir, PE_dir, filestr);
end

return
