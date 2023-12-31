function  log = compute_gnss_ecef(p,eph,obs)
% This function is to implement GNSS positioning with
% standard mode (without Iono, Trop, Es correction)
% or PPP mode.
% Output: log is a data struct that store the results.
%----------------------------------%
N = length(obs.tr_sow); % The number of positioning points
log = initOutputLog(p, obs);
% Mark the sat prn that be computed
gpslog.svprn_mark = zeros(log.num_obs_gps,1);glolog.svprn_mark = zeros(log.num_obs_glo,1);
gallog.svprn_mark = zeros(log.num_obs_gal,1);bdslog.svprn_mark = zeros(log.num_obs_bds,1);
% Record the prn of each system for whose satllite been used
gpslog.prn_record = zeros(log.num_obs_gps,1);glolog.prn_record = zeros(log.num_obs_glo,1);
gallog.prn_record = zeros(log.num_obs_gal,1);bdslog.prn_record = zeros(log.num_obs_bds,1);
% Satellite position in ECEF frame
gpslog.s_pos_ecef = [];glolog.s_pos_ecef = [];gallog.s_pos_ecef = [];bdslog.s_pos_ecef = [];
if p.post_mode == 1
    % precise satellite position in ECEF frame
    gpslog.s_pos_prc = [];glolog.s_pos_prc = [];gallog.s_pos_prc = [];bdslog.s_pos_prc = [];
end
% Satellite velocity in ECEF frame
gpslog.s_v_ecef = [];glolog.s_v_ecef = [];gallog.s_v_ecef = [];bdslog.s_v_ecef = [];
% Signal propagation time
gpslog.tp = [];glolog.tp = [];gallog.tp = [];bdslog.tp = [];
% corrected pseudorange
gpslog.corr_range = [];glolog.corr_range = [];gallog.corr_range = [];bdslog.corr_range = [];
% The number of satellite be computed
gpslog.num_sv = 0;glolog.num_sv = 0;gallog.num_sv = 0;bdslog.num_sv = 0;
for i = 1:p.inval:N
    if ~isnan(p.Grdpos.t(1))
        index = round(p.Grdpos.t) == round(obs.tr_sow(i));
        grdpos = p.Grdpos.pos(:,index);
        if isempty(grdpos)
            continue;
        end
    else
        grdpos = p.Grdpos.pos;
    end
    if mod(i,4000)==0
        i
    end
    p.i = i; % To debug
    if ~isempty(obs.gps)&& p.freq==1 && p.enableGPS == 1
        % GPS satellite position computation, Single frenquency receiver mode
        gpslog = satpost_corrpsedR_singlefreq(p,eph,obs,i,log.num_obs_gps,'gps');
    end
    if ~isempty(obs.glo)&& p.freq==1 && p.enableGLO == 1
        % GLO satellite position computation, Single frenquency receiver mode
        glolog = satpost_corrpsedR_singlefreq(p,eph,obs,i,log.num_obs_glo,'glo');
    end
    if ~isempty(obs.gal)&& p.freq==1 && p.enableGAL == 1
        % GAL satellite position computation, Single frenquency receiver mode
        gallog = satpost_corrpsedR_singlefreq(p,eph,obs,i,log.num_obs_gal,'gal');
    end
    if ~isempty(obs.bds)&& p.freq==1 && p.enableBDS == 1
        % BDS satellite position computation, Single frenquency receiver mode
        bdslog = satpost_corrpsedR_singlefreq(p,eph,obs,i,log.num_obs_bds,'bds');
    end
    cpt.prn_record = [gpslog.prn_record;glolog.prn_record;gallog.prn_record;bdslog.prn_record];
    cpt.svprn_mark = [gpslog.svprn_mark;glolog.svprn_mark;gallog.svprn_mark;bdslog.svprn_mark];
    cpt.corr_range = [gpslog.corr_range;glolog.corr_range;gallog.corr_range;bdslog.corr_range];
    ind = find(cpt.prn_record==0);
    cpt.prn_record(ind) = [];
    cpt.svprn_mark(ind) = [];
    ind = find(cpt.corr_range==0);
    cpt.corr_range(ind) = [];
    cpt.s_pos_ecef = [gpslog.s_pos_ecef,glolog.s_pos_ecef,gallog.s_pos_ecef,bdslog.s_pos_ecef];
    cpt.s_pos_ecef(:,ind) = [];
    if p.post_mode == 1
        cpt.s_pos_prc = [gpslog.s_pos_prc,glolog.s_pos_prc,gallog.s_pos_prc,bdslog.s_pos_prc];
        cpt.s_pos_prc(:,ind) = [];
    end
    cpt.s_v_ecef = [gpslog.s_v_ecef,glolog.s_v_ecef,gallog.s_v_ecef,bdslog.s_v_ecef];
    cpt.s_v_ecef(:,ind) = [];
    cpt.tp = [gpslog.tp;glolog.tp;gallog.tp;bdslog.tp];
    cpt.tp(ind) = [];
    cpt.num_sv = [gpslog.num_sv,glolog.num_sv,gallog.num_sv,bdslog.num_sv];
    if sum(cpt.num_sv) < p.min_sv + sum(cpt.num_sv~=0) - 1 ...
            || (p.enableGPS  == 1 && gpslog.num_sv == 0)
        continue;
    end
    % elevation & azimuth
    cpt.elev = NaN(length(cpt.corr_range),1); cpt.az = NaN(length(cpt.corr_range),1);
    % trop delay and iono delay
    cpt.trop_delay = NaN(length(cpt.corr_range),1); cpt.iono_delay = NaN(length(cpt.corr_range),1);
    if isempty(log.epoch_t) || (seconds(obs.datetime(i) - log.epoch_t(end)) > 1.5 && p.post_mode ~= p.mode_sps)
        [estState,res] = userpos(p,cpt);
        if p.post_mode == 1
            p.mk = 1;
        end
        if p.post_mode == p.mode_dgnss && ~isempty(p.eph_b) && ~isempty(p.obs_b)
            [cpt,n] = diff_corr_compute(p,cpt,obs.tr_posix(i));
            if ~isempty(n) && sum(cpt.num_sv)>=p.min_sv
                cpt = cpt_clear(cpt); % Clear the data where don't have diff correction
                cpt.corr_range = cpt.corr_range - cpt.diff_corr;
                [estState,res] = userpos(p,cpt);
            end
        end
        if ~isempty(estState.pos)
            % Save the initial state but not the first epoch.
            [p.state0, p.state_cov] = obtainInitEkfStateAndCov(p, estState);
            log.epoch_t = [log.epoch_t, obs.datetime(i)];
            %log = save_result(p,cpt,log,i,estState,res,grdpos,obs.datetime(i));
        end
        dt = 0;
    elseif p.post_mode ~= p.mode_sps
        dt = seconds(obs.datetime(i) - log.epoch_t(end));
        % EKF predict
        [p.state0, p.state_cov] = ekfPredict(p, p.state0, p.state_cov, dt);
    end
    % Rotate the sat pos to common reference frame
    [estState,~] = userpos(p,cpt);
    cpt = earth_rotation_corr(p,cpt,estState.clock_bias/p.c);
    % Check elevation
    cpt = elevaz_check(p,cpt,p.state0(1:3));
    if sum(cpt.num_sv) <= p.min_sv + sum(cpt.num_sv~=0) - 1
        continue;
    end
    switch p.post_mode
        case p.mode_sps % Standard GNSS
            % if p.elev_mark ==0
            % log = save_result(p,cpt,log,i,re_pos,clock_bias,res);
            % else
            % Compute the final position
            % [re_pos,clock_bias,res] = userpos_Rcorr(p,cpt);
            tdoy = doy(obs.tr_prime(1:3,i)); % Day of year
            [rt.week, rt.dow, rt.sow] = date2gnsst(obs.tr_prime(:,i)');
            cpt.IoFac = zeros(length(cpt.corr_range),1);
            cpt = trop_iono_compute(p,eph,cpt,obs,p.state0(1:3),tdoy,[],rt);
            cpt.corr_range = cpt.corr_range - cpt.trop_delay - cpt.iono_delay;
            [estState,res] = userpos(p,cpt);
            log = save_result(p,cpt,log,i,estState,res,grdpos,obs.datetime(i));
            % end
        case p.mode_ppp % PPP
            tdoy = doy(obs.tr_prime(1:3,i)); % Day of year
            [rt.week, rt.dow, rt.sow] = date2gnsst(obs.tr_prime(:,i)');
            %%-------------%%
            cpt.IoFac = zeros(length(cpt.corr_range),1);
            cpt = correctBeidouCodeError(p,cpt);
            cpt = trop_iono_compute(p,eph,cpt,obs,p.state0(1:3),tdoy,rt,obs.tr_posix(i));
            % Using the correction to the measurements
            if ~isempty(find(cpt.iono_delay~=0, 1))
                cpt.corr_range = cpt.corr_range - cpt.trop_delay - cpt.iono_delay;
                %%-------------%%
                % Compute the final position
                if p.double_diff == 0
                    %[estState,res] = userpos(p,cpt);
                    [p,estState,res] = stateUpdate(p,cpt,dt);
                else
                    [re_pos,clock_bias,res] = userpos_2diff(p,cpt);
                end
                if ~isempty(estState.pos)
                    log = save_result(p,cpt,log,i,estState,res,grdpos,obs.datetime(i));
                end
            end
        case p.mode_dgnss % DGNSS
            if ~isempty(p.eph_b) && ~isempty(p.obs_b)
                [cpt,n] = diff_corr_compute(p,cpt,obs.tr_posix(i));
                if ~isempty(n) && sum(cpt.num_sv)>=p.min_sv
                    cpt = cpt_clear(cpt); % Clear the data where don't have diff correction
                    cpt.corr_range = cpt.corr_range - cpt.diff_corr;
                    [p,estState,res] = stateUpdate(p,cpt,dt);
                    %[estState,res] = userpos(p,cpt);
                    if ~isempty(estState.pos)
                        log = save_result(p,cpt,log,i,estState,res,grdpos,obs.datetime(i));
                    end
                end
            else
                warning('No differential source given')
            end
        otherwise
            warning('Unsupport positioning option');
    end
end






end