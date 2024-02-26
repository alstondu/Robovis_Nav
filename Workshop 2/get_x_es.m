function [r_es, v_es] = get_x_es()
    Pseudo_ranges = readtable("Workshop2_Pseudo_ranges.csv");
    num_s = table2array(Pseudo_ranges(1,2:end));
    time = table2array(Pseudo_ranges(2:end,1));
    r_es = zeros(length(num_s),3, length(time));
    v_es = zeros(length(num_s),3, length(time));
    for i = 1: length(time)
        for j = 1:length(num_s)
            [r_es(j, :, i),v_es(j, :, i)] = Satellite_position_and_velocity(time(i),num_s(j));
        end
    end
end
