function [deg1, deg2] = rot_mat_to_cts_orientation_string(R)

    sin_b = -R(1, 3);

    cos_a_cos_b = R(2, 3);

    sin_a_cos_b = -R(3, 3);

    alpha2 = -asin(R(1, 3));
    alpha1 = -acos(R(2, 3) / cos(alpha2));

    deg1 = rad2deg(alpha1);
    deg2 = rad2deg(alpha2);

end
