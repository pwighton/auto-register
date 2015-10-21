function R = get_plane_rotation(or1, or2, alpha, is_degrees)

  or1 = lower(or1);
  or2 = lower(or2);

  if (nargin > 3 & is_degrees)
      alpha = deg2rad(alpha);
  end

  if (strcmp(or1, 'cor'))
      if (strcmp(or2, 'sag'))
          R = get_z_rot(alpha);
      elseif (strcmp(or2, 'tra'))
          R = get_x_rot(-alpha);
      end
  elseif (strcmp(or1, 'sag'))
      R = get_z_rot(pi/2);
      if (strcmp(or2, 'cor'))
          R = R * get_z_rot(-alpha);
      elseif (strcmp(or2, 'tra'))
          R = R * get_x_rot(alpha);
      else
          R = [];
      end
      R = get_z_rot(-pi/2) * R;
  elseif (strcmp(or1, 'tra'))
      R = get_y_rot(pi/2);
      if (strcmp(or2, 'cor'))
          R = R * get_y_rot(alpha);
      elseif (strcmp(or2, 'sag'))
          R = R * get_x_rot(alpha);
      else
          R = [];
      end
      R = get_y_rot(-pi/2) * R;
  end
end
