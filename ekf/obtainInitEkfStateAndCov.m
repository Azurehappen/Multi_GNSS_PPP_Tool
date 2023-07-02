function [state, cov] = obtainInitEkfStateAndCov(p, estState)

% position, clock bias (m), clock drift (m/s)
state = [estState.pos;estState.clock_bias;0];
cov_diag = [30^2 * ones(1,3), 40^2, 200^2];
if ~isnan(estState.isb_dict(p.glo.sys_num))
    state = [state; estState.isb_dict(p.glo.sys_num)];
    cov_diag = [cov_diag, 5^2];
end
if ~isnan(estState.isb_dict(p.gal.sys_num))
    state = [state; estState.isb_dict(p.gal.sys_num)];
    cov_diag = [cov_diag, 5^2];
end
if ~isnan(estState.isb_dict(p.bds.sys_num))
    state = [state; estState.isb_dict(p.bds.sys_num)];
    cov_diag = [cov_diag, 5^2];
end

cov = diag(cov_diag);
