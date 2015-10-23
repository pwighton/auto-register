function R = rtc(angle, is_degrees)

  if (nargin > 1 && is_degrees)
      angle = deg2rad(angle);
  end

  R = [1          0           0;
       0 cos(angle) -sin(angle);
       0 sin(angle)  cos(angle)];
end
