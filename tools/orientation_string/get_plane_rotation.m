function R = get_plane_rotation(or1, or2, alpha, is_degrees)

  or1 = lower(or1);
  or2 = lower(or2);

  if (nargin > 3 & is_degrees)
      alpha = deg2rad(alpha);
  end

  rstr = sprintf('R = r%s%s(alpha);', or1, or2);
  eval(rstr);

end
