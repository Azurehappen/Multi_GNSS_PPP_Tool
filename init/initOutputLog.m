function log = initOutputLog(p, obs)

N = length(obs.tr_sow); % The number of positioning points
% Initialize output
log.epoch_t = datetime.empty;
log.gpst = obs.tr_sow-obs.tr_sow(1);
log.err = NaN(1,N); % The position (Norm) error between estimated pos and true pos
log.hor_err = NaN(1,N); % The horizontal position (Norm) error between estimated pos and true pos
log.ned_err_norm = NaN(1,N); % NED frame norm error
log.ned_err = NaN(3,N); % NED frame error
log.pos_ecef = NaN(3,N); % Estimated position in ECEF
log.rover_clk = NaN(1,N); % Receiver clock bias (meter)
log.sv_num_GPS = NaN(1,N); % The amount of GPS satellite be used
log.sv_num_GLO = NaN(1,N); % The amount of GLO satellit  used
log.sv_num_GAL = NaN(1,N); % The amount of GAL satellite be used
log.sv_num_BDS = NaN(1,N); % The amount of BDS satellite be used
log.isb_glo = NaN(1,N);
log.isb_gal = NaN(1,N);
log.isb_bds = NaN(1,N);

numOfState = p.modeToNumUserStates(p.state_mode) + 2; % user_states,clk,clk_drift
log.state_cov = NaN(numOfState+p.enableGLO+p.enableGAL+p.enableBDS, N);
log.ned_cov = NaN(3, N);
if p.est_mode == p.raps_est
    log.state_info = log.state_cov;
    log.raps_spec_xyz = NaN(3, N);
    log.pos_info_ned = NaN(3, N);
end
if ~isempty(obs.gps)
    log.num_obs_gps = size(obs.gps(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_gps = 0;
end
if ~isempty(obs.glo)
    log.num_obs_glo = size(obs.glo(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_glo = 0;
end
if ~isempty(obs.gal)
    log.num_obs_gal = size(obs.gal(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_gal = 0;
end
if ~isempty(obs.bds)
    log.num_obs_bds = size(obs.bds(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_bds = 0;
end
log.res_GPS = NaN(log.num_obs_gps,N); % The residual at the end
log.res_GLO = NaN(log.num_obs_glo,N); % The residual at the end
log.res_GAL = NaN(log.num_obs_gal,N); % The residual at the end
log.res_BDS = NaN(log.num_obs_bds,N); % The residual at the end
log.elev_GPS = NaN(log.num_obs_gps,N); % The elevation of satellites
log.elev_GLO = NaN(log.num_obs_glo,N);
log.elev_GAL = NaN(log.num_obs_gal,N);
log.elev_BDS = NaN(log.num_obs_bds,N);
log.res = [log.res_GPS;log.res_GAL;log.res_GLO;log.res_BDS];
log.elev = [log.elev_GPS;log.elev_GAL;log.elev_GLO;log.elev_BDS];

end