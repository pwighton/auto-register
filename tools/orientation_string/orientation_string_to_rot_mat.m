function R = orientation_string_to_rot_mat(or1, or2, deg1, or3, deg2)

  var_info = whos('deg1');
  if (strcmp(var_info.class, 'sym'))
      is_degrees = 0;
  else
      is_degrees = 1;
  end

  R1 = get_plane_rotation(or1, or2, deg1, is_degrees);
  R2 = get_plane_rotation(or1, or3, deg2, is_degrees);
  R = R2 * R1;
end
