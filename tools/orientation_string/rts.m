function R = rts(angle, is_degrees)

  if (nargin > 1 && is_degrees)
      angle = deg2rad(angle);
  end

  R = [ cos(angle) 0 -sin(angle);
                 0 1           0;
        sin(angle) 0  cos(angle)];
end
