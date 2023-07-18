function log = save_result(p,cpt,log,i,estState,res,grdpos,epoch_t)

log.epoch_t = [log.epoch_t, epoch_t];
log.pos_ecef(:,i) = estState.pos;
%------------------------%
% [pos_llh,~,~]=ecef2llh_iter(estState.pos);
% R_e2g=ll2R(pos_llh); % rotation matrix from ecef 2 geodetic frame
% err_pos = grdpos - estState.pos;
% ned_err=R_e2g*err_pos;
lla_gt = ecef2lla(grdpos', 'WGS84');
wgs84 = wgs84Ellipsoid('meter');
pos = estState.pos;
[xNorth,yEast,zDown] = ecef2ned(pos(1),pos(2),pos(3),lla_gt(1),lla_gt(2),lla_gt(3),wgs84);
ned_err = [xNorth;yEast;zDown];
log.ned_err(:,i) = ned_err;
%------------------------%
log.ned_err_norm(i) = norm(ned_err);
log.hor_err(i) = norm(ned_err(1:2));
log.err(i) = norm(grdpos - estState.pos);
log.rover_clk(i) = estState.clock_bias;
log.sv_num_GPS(i) = cpt.num_sv(1);log.sv_num_GLO(i) = cpt.num_sv(2);
log.sv_num_GAL(i) = cpt.num_sv(3);log.sv_num_BDS(i) = cpt.num_sv(4);
ind_mark = cpt.svprn_mark ~= 0;
log.res(ind_mark,i) = res;
log.elev(ind_mark,i) = cpt.elev;
start = 1; endi = log.num_obs_gps;
log.res_GPS(:,i) = log.res(start:endi,i);
log.elev_GPS(:,i) = log.elev(start:endi,i);
start = start + log.num_obs_gps; endi = endi + log.num_obs_glo;
log.res_GLO(:,i) = log.res(start:endi,i);
log.elev_GLO(:,i) = log.elev(start:endi,i);
start = start + log.num_obs_glo; endi = endi + log.num_obs_gal;
log.res_GAL(:,i) = log.res(start:endi,i);
log.elev_GAL(:,i) = log.elev(start:endi,i);
start = start + log.num_obs_gal; endi = endi + log.num_obs_bds;
log.res_BDS(:,i) = log.res(start:endi,i);
log.elev_BDS(:,i) = log.elev(start:endi,i);
log.isb_glo(i) = estState.isb_dict(p.glo.sys_num);
log.isb_gal(i) = estState.isb_dict(p.gal.sys_num);
log.isb_bds(i) = estState.isb_dict(p.bds.sys_num);
log.state_cov(:,i) = diag(p.state_cov)';
R_e2g=ll2R(lla_gt');
ned_cov = R_e2g * p.state_cov(1:3, 1:3) * R_e2g';
log.ned_cov(:,i) = diag(ned_cov);
if p.est_mode == p.raps_est
    log.state_info(:,i) = diag(p.raps_J)';
    log.raps_spec_xyz(:,i) = p.raps_spec_xyz;
    J_ned = R_e2g' * p.raps_J(1:3, 1:3) * R_e2g;
    log.pos_info_ned(:,i) = diag(J_ned);
end
end