function segments = get_segments(mask, min_len)
    mask = mask(:)';
    d = diff([0 mask 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    lengths = ends - starts + 1;
    valid = lengths >= min_len;
    segments = [starts(valid); ends(valid)]';
end