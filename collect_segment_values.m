function values = collect_segment_values(signal, segments)
    values = [];
    for i = 1:size(segments,1)
        idx = segments(i,1):segments(i,2);
        values = [values; signal(idx)];
    end
end