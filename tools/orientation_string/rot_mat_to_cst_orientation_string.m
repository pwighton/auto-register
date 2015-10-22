function [or1, or2, deg1, or3, deg2] = rot_mat_to_cst_orientation_string(R)

    alpha1 = atan2(-R(1, 3), R(1, 1));

    R1 = get_plane_rotation('cor', 'sag', alpha1);

    R2 = R * inv(R1);

    alpha2 = atan2(R2(2, 3), R2(2, 2));

    or1 = 'cor';
    or2 = 'sag';
    or3 = 'tra';

    deg1 = rad2deg(alpha1);
    deg2 = rad2deg(alpha2);

end
