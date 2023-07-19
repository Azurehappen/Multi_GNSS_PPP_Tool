function [prior_state, prior_cov] = ekfPredict(p, state, cov, dt)

Qc = p.ekf_para.q_clkDrift * [dt^3 / 3, dt^2 / 2;
                   dt^2 / 2, dt];
Qp_diag = p.ekf_para.q_pos * ones(1,3);
Qv_diag = p.ekf_para.q_vel * ones(1,3);
Qa_diag = p.ekf_para.q_acc * ones(1,3);
num_user_states = p.modeToNumUserStates(p.state_mode);
phi = diag([ones(1, num_user_states+2), ...
     exp(p.ekf_para.u_isb * dt) * ones(1, length(state) - num_user_states-2)]);
if p.state_mode == p.pva_mode
    % x = [x,y,z,vx,vy,vz,ax,ay,az,clk,clk_drift,isb..]
    phi(1:9,1:9) = [eye(3,3),dt*eye(3,3), dt^2/2 * eye(3,3);
                zeros(3,3),eye(3,3),dt*eye(3,3);
                zeros(3,3),zeros(3,3),eye(3,3)];
    % Qd_p = Q*dt+(F*Q+Q*F')*dt^2/2 + (...)*dt^3/3, see Aided Nav. Book
    % eqn. 4.119
    Q_p = [diag(Qp_diag)*dt + diag(Qv_diag)*dt^3/3, diag(Qv_diag)*dt^2/2+diag(Qa_diag)*dt^4/8, diag(Qa_diag)*dt^3/6;
           diag(Qv_diag)*dt^2/2+diag(Qa_diag)*dt^4/8, diag(Qa_diag)*dt^3/3 + diag(Qv_diag)*dt, diag(Qa_diag)*dt^2/2;
           diag(Qa_diag)*dt^3/6, diag(Qa_diag)*dt^2/2, diag(Qa_diag)*dt];
    Q_pc = [Q_p, zeros(num_user_states, 2);
            zeros(2, num_user_states), Qc];
elseif p.state_mode == p.pos_mode
     % x = [x,y,z,clk,clk_drift,isb..]
    Q_pc = [diag(Qp_diag), zeros(num_user_states, 2);
            zeros(2, num_user_states), Qc];
end
phi(num_user_states+1, num_user_states+2) = dt; % clk(k) = clk(k-1) + clk_drift * dt

Lbias = length(state) - num_user_states - 2;
if Lbias ~= 0
    Qd = [Q_pc, zeros(num_user_states+2, Lbias);
          zeros(Lbias, num_user_states+2), p.ekf_para.q_isb * eye(Lbias, Lbias)];
else
    Qd = Q_pc;
end

prior_state = phi * state;
prior_cov = phi * cov * phi' + Qd;