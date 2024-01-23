function r_ea = get_r_ea()
    L_b = deg2rad(-33.821075);
    lambda_b = deg2rad(151.188496);
    h_b = 120;
    v_eb_n = 0;
    [r_ea,v_ea] = pv_NED_to_ECEF(L_b,lambda_b,h_b,v_eb_n);
end

