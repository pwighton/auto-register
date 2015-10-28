function T = read_lta(filename)

    T = eye(4);
    cur_row = 1;

    f = fopen(filename);
    while (~feof(f))
        fline = fgets(f);
        if (~(isstrprop(fline(1), 'alphanum') || ...
              fline(1) == '-') || ...
            isletter(fline(1)) || ...
            strcmp(fline(1:5), '1 4 4'))

            continue
        end

        T(cur_row, :) = sscanf(fline, '%f %f %f %f')';
        cur_row = cur_row + 1;
    end


end
