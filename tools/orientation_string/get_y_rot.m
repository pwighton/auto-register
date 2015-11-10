function R = get_y_rot(alpha)
  R = [ cos(alpha) 0  -sin(alpha);
                 0 1           0;
        sin(alpha) 0  cos(alpha)];
end
