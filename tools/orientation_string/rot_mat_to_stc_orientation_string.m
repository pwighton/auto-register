function [deg1, deg2, deg3] = rot_mat_to_stc_orientation_string(R)

    gamma = atan2(R(2, 2), R(2, 1));
    alpha1 = atan2(-R(3, 3), R(1, 3));
    alpha2 = atan2(-R(2, 3), sqrt(R(1,3)^2 + R(3, 3)^2));

    deg1 = rad2deg(alpha1);
    deg2 = rad2deg(alpha2);
    deg3 = rad2deg(gamma);

end
