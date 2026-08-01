function [mac] = var_MACdd(psi_m,phi_m,cov_psi_csqrt_ref,cov_phi_csqrt,N)

        n_dof = length(phi_m);        
        I_d = eye(n_dof);

        % Build the Hessian for MAC

        ipr = psi_m'*psi_m;

        m_xx_psi_psi = real(psi_m)*real(psi_m)';
        m_yy_psi_psi = imag(psi_m)*imag(psi_m)';

        m_xy_psi_psi = real(psi_m)*imag(psi_m)';
        m_yx_psi_psi = m_xy_psi_psi';

        k_limit = psi_m'*phi_m/(ipr) ;
        
        kd_limit = 2/(abs(k_limit)^2*ipr);
        h1_limit = [real(k_limit)*I_d -imag(k_limit)*I_d;imag(k_limit)*I_d real(k_limit)*I_d];
        
        h2_limit = [m_xx_psi_psi+m_yy_psi_psi-ipr*I_d m_xy_psi_psi-m_yx_psi_psi ; m_yx_psi_psi-m_xy_psi_psi m_xx_psi_psi+m_yy_psi_psi-ipr*I_d];
        p_limit = [-eye(2*n_dof) zeros(2*n_dof);zeros(2*n_dof) h1_limit.'];
        
        H_mac_limit = -kd_limit*p_limit*[h2_limit h2_limit;h2_limit h2_limit]*p_limit/2; 
        
        stack_sqrtcov_psi_phi = [cov_psi_csqrt_ref;cov_phi_csqrt];
        
        a_sig = H_mac_limit*N*(stack_sqrtcov_psi_phi*stack_sqrtcov_psi_phi');
        c_2 = trace((a_sig)^2);
        c_3 = trace((a_sig)^3);
        s_1 = c_3/c_2^(1.5);

        mu_q = trace(a_sig);
        sig_q = sqrt(2*c_2);

        a = 1/s_1;
        delta_h = 0; %
        l_h = a^2; % 

        a = sqrt(l_h + 2*delta_h);

        mu_chi2 = l_h + delta_h;
        std_chi2 = sqrt(2)*a;
        
        alpha_ = sig_q./std_chi2;
        beta_ = (mu_q - (mu_chi2.*sig_q)./std_chi2); 
        
        qq = ncx2inv(1-eps,l_h,delta_h) - ncx2inv(0.01,l_h,delta_h);
        tmac = 1 - beta_/N - alpha_/N.*qq;
        
        if tmac<0
            tmac= 1;
        end
        
        meval = (abs(psi_m'*phi_m)).^2 ./ (psi_m'*psi_m * diag(phi_m'*phi_m))';
        
        mac.meval = meval;
        mac.tmac = tmac;
        mac.bmac = meval>=tmac;

    
end




