function [orientation_str, Tx, Ty, Tz] = affine_to_orientation_string(T)

    orientation_str = rot_mat_to_orientation_string(T);
    Tx = T(1, 4);
    Ty = T(2, 4);
    Tz = T(3, 4);

end
