function delta_rho = delta_rho_function(x, x_d, E0, E1, a, t, w)
    G0_terms = G0_function(x, x_d, a, t, E0, w) .* G0_function(x_d, x_d, a, t, E0, w);
    E1_minus_G0 = E1 - G0_function(x_d, x_d, a, t, E0, w);
    delta_rho = -1/pi * imag(G0_terms ./ E1_minus_G0);
end