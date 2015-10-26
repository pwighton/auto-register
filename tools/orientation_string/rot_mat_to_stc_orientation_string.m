function [deg1, deg2] = rot_mat_to_stc_orientation_string(R)

    alpha2 = atan2(R(1, 1), R(2, 1));

    sin_a = -R(1, 2) / cos(alpha2);

    alpha1 = atan2(sin_a, -R(3, 2));

    deg1 = rad2deg(alpha1);
    deg2 = rad2deg(alpha2);

end
