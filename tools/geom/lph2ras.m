function ras = lph2ras(lph)
% convert an affine transform from LHP to RAS coordinates

    flip = eye(4);
    flip(1, 1) = -1;
    flip(2, 2) = -1;

    ras = flip * lph * flip;
