function R = get_x_rot(alpha)
  R = [1          0           0;
       0 cos(alpha) -sin(alpha);
       0 sin(alpha)  cos(alpha)];
 end
