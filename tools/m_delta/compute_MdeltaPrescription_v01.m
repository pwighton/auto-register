function compute_MdeltaPrescription_v01(pos, RO_dir, PE_dir, filestr)
% compute_MdeltaPrescription_v01(pos, RO_dir, PE_dir, file_str)
% -- pos: the desired slice center in LPH coords
% -- RO_dir: the desired readout direction in LPH coords
% -- PE_dir: the desired phase-encoding direction in LPH coords
% -- filestr: e.g. 'mat01.txt' (for auto_register.py's -prescrip argument)
%
% NOTE: the initial prescription (which AAhijack will subsequently update)
% must be SAGITTAL (for now, at least)
%
% (mukundb, 2025/10/06)

% reshape the inputs into column vectors
pos = pos(:);
RO_dir = RO_dir(:);
PE_dir = PE_dir(:);

% normalize the direction vectors
RO_dir = RO_dir / norm(RO_dir);
PE_dir = PE_dir / norm(PE_dir);

% check that RO and PE are orthogonal
if abs(dot(PE_dir, RO_dir)) > sqrt(eps)
  error('RO and PE dirs need to be orthogonal to each other!')
end

% calculate the SS direction
SS_dir = cross(PE_dir, RO_dir);

% build the 3x3 "desired prescription" matrix
Mprescription = [PE_dir RO_dir SS_dir];

% build the 3x3 "initial prescription" matrix (NOTE: the below is only
% correct if the initial prescription is chosen to be the Siemens default
% SAGITTAL: at isocenter with PE: A >> P and 0 degrees in-plane rotation)
Minitial = [  0  0 -1
             +1  0  0
              0 -1  0];
%            PE RO SS

% now we can compute MdeltaPrescription (the CHANGE in prescription)
MdeltaPrescription = Mprescription / Minitial;
% i.e., MdeltaPrescription = Mprescription * inv(Minitial);
% i.e., MdeltaPrescription * Minitial = Mprescription;

% make it a 4x4 transformation matrix
MdeltaPrescription(:,4) = pos;       % 4th col: add translation vector
MdeltaPrescription(4,:) = [0 0 0 1]; % 4th row: for homogenous coords

% write out MdeltaPrescription either to stdout or to a file
if nargin < 4
  fid = 1; % stdout
else
  fid = fopen(filestr, 'wt');
end

fprintf(fid, '%f %f %f %f\n', MdeltaPrescription(1,:));
fprintf(fid, '%f %f %f %f\n', MdeltaPrescription(2,:));
fprintf(fid, '%f %f %f %f\n', MdeltaPrescription(3,:));
fprintf(fid, '%f %f %f %f\n', MdeltaPrescription(4,:));

if nargin < 4
  fprintf('NOTE: above matrix NOT written to file\n');
else
  fclose(fid);
  fprintf('matrix written to file: %s\n', filestr);
end

return
