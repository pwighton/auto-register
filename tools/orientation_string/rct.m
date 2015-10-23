function R = rct(angle, is_degrees)

  if (nargin > 1 && is_degrees)
      angle = deg2rad(angle);
  end

  R = [1           0           0;
       0 -sin(angle)  cos(angle);
       0 -cos(angle) -sin(angle)];
end
