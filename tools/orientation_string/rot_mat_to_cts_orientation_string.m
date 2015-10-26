function [deg1, deg2] = rot_mat_to_cts_orientation_string(R)

    alpha1 = atan2(-R(3, 3), -R(3, 2));
    alpha2 = atan2(R(2, 1), R(1, 1));

    deg1 = rad2deg(alpha1);
    deg2 = rad2deg(alpha2);

end
