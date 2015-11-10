function M = header_to_m(ref_mgh_file, mov_mgh_file)

  [v, h_ref] = load_mgh(ref_mgh_file);
  [v, h_mov] = load_mgh(mov_mgh_file);

  M = h_ref * inv(h_mov);
