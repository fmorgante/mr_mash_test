out <- simulate_data(n, p, p_causal, r, intercepts=rep(1, r),
                     pve, Sigma_cor_offdiag, Sigma_scale,
                     Gamma_cor_offdiag, Gamma_scale,
                     V_cor_offdiag, V_offdiag_scale, prop_testset)