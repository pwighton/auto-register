function lhp = ras2lhp(ras)
% convert an affine transform from RAS to LHP coordinates

    flip = eye(4);
    flip(1, 1) = -1;
    flip(2, 2) = -1;

    lhp = flip * ras * flip;
