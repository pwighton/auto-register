function [deg1, deg2] = rot_mat_to_tsc_orientation_string(R)

    alpha1 = atan2(-R(1, 3), R(1, 1));
    alpha2 = atan2(R(3, 2), R(2, 2));

    deg1 = rad2deg(alpha1);
    deg2 = rad2deg(alpha2);

end
