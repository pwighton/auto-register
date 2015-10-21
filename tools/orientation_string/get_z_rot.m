function R = get_z_rot(alpha)
  R = [cos(alpha) -sin(alpha) 0;
       sin(alpha)  cos(alpha) 0;
                0           0 1];
end
