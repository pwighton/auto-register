function R = orientation_string_to_rot_mat(or1, or2, deg1, or3, deg2)

  R1 = get_plane_rotation(or1, or2, deg1, 1)

  R2 = get_plane_rotation(or1, or3, deg2, 1);

  R = R2 * R1;
end
