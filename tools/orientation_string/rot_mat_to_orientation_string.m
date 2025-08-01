function orientation_str = rot_mat_to_orientation_string(R)

    orientation_str = '';

    orders = perms({'t', 's', 'c'});
    for order=1:length(orders)
        eval_str = sprintf('[a, b, g] = rot_mat_to_%s_orientation_string(R);', ...
                           strcat(orders{order, :}));
        eval(eval_str);

        if abs(g) < 45 && abs(a) <= 45 && abs(a) > abs(b)
            if length(orientation_str) > 0
                fprintf('ambiguous orientation_string, already had %s\n', ...
                        orientation_str);
            end

            orientation_str = sprintf('%s > %s (%0.2f) > %s (%0.2f)', ...
                                      orders{order, 1}, ...
                                      orders{order, 2}, ...
                                      a, ...
                                      orders{order, 3}, ...
                                      b);
        end
    end
end
