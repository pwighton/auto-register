function [deg1, deg2] = rot_mat_to_cst_orientation_string(R)

    alpha2 = -asin(R(3, 3));
    alpha1 = -acos(R(2, 3) / cos(alpha2));

    deg1 = rad2deg(alpha1);
    deg2 = rad2deg(alpha2);

end
