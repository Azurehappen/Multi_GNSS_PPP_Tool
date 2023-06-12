function [satlog] = satpost_corrpsedR_singlefreq(p,eph,obs,ind,len_prn,sys_type)
% This function is to compute the satellite positions and correct 
% the psedoranges by sat clock bias
% Input: 
%       p --parameters
%       eph --ephemeris data
%       obs --observables
%       ind --index of observables data
%       sys_type --The system that be computed, 'gps', 'gal', 'glo' and 'bds'
% Outpu:
%       satlog.svprn_mark -- Mark the sat prn that be computed
%       satlog.s_pos_ecef -- Satellite position in ECEF frame
%       satlog.corr_range -- corrected pseudorange
%----------------------------%
% Initialize
satlog.svprn_mark = zeros(len_prn,1);
satlog.prn_record = zeros(len_prn,1);
satlog.s_pos_ecef = zeros(3,len_prn);
satlog.s_pos_prc = zeros(3,len_prn);
satlog.s_v_ecef = zeros(3,len_prn);
satlog.corr_range = zeros(len_prn,1);
satlog.tp = zeros(len_prn,1);
obs_tr_gps = obs.tr_posix(ind);
obs_tr = obs_tr_gps;
if p.post_mode == 1 && ~isempty(p.orbit_dict) && ~isempty(p.clock_dict)
    [~, orbit_single] = closestKeyAndValue(p.orbit_dict, obs_tr_gps, 0, 60);
    if isempty(orbit_single)
        satlog.num_sv = 0;
        return;
    end
end

switch(sys_type)
    case 'gps'
        obs_range = obs.gps(1).data.P(:,ind);
        Strength = obs.gps(1).data.S(:,ind);
        message_duration = p.gps.message_duration;
        eph_info = eph.gps;
    case 'glo'
        obs_range = obs.glo(1).data.P(:,ind);
        Strength = obs.glo(1).data.S(:,ind);
        message_duration = p.glo.message_duration;
        obs_tr = obs_tr_gps - p.glo.lps_gps; % Correct time diff from GPS time to GLO time
        eph_info = eph.glo;
    case 'gal'
        obs_range = obs.gal(1).data.P(:,ind);
        Strength = obs.gal(1).data.S(:,ind);
        message_duration = p.gal.message_duration;
        obs_tr = obs_tr_gps - p.gal.lps_gps; % Correct time diff from GPS time to GAL time
        eph_info = eph.gal;
    case 'bds'
        obs_range = obs.bds(1).data.P(:,ind);
        Strength = obs.bds(1).data.S(:,ind);
        message_duration = p.bds.message_duration;
        obs_tr = obs_tr_gps - p.bds.lps_gps; % Correct time diff from GPS time to BDS time
        eph_info = eph.bds;
end

for prn = 1 :len_prn
    if p.post_mode == 1 && (~isfield(orbit_single, sys_type)...
            || ~isKey(orbit_single.(sys_type), prn))
        continue;
    end
    orbit_corr = [];
    clock_corr_m = 0;
    if p.post_mode == 1
        orbit_corr = orbit_single.(sys_type)(prn);
        clock_corr_m = obtainSatClockCorr(p.clock_dict, obs_tr_gps, 0, 60, prn, ...
            orbit_corr.iod, sys_type);
        if isnan(clock_corr_m)
            continue;
        end
    end

    if (obs_range(prn)~=0)&&(~isnan(obs_range(prn)))&& prn<=size(eph_info.a_f0,1)
        tp_prime = obs_range(prn)/p.c;
        t_sv = obs_tr-tp_prime;
        tidx = ephtidx(eph_info.t_oc{prn},t_sv,eph_info.SV_health(prn,:),message_duration);
        % Check the signal strength and sv health (health message 0 means ok)
        if ~isempty(tidx) && Strength(prn)>=p.sig_strg
            switch(sys_type)
                % compute ephemeris to satellite position and clock bias
                case 'gps'
                    [sat, dt_sv] = eph2pos(p,eph_info,prn,tidx,t_sv,'gps',orbit_corr);
                    if ~isnan(sat.pos_ecef(1))
                        satlog.svprn_mark(prn) = p.gps.sys_num;satlog.prn_record(prn) = prn;
                    end
                case 'glo'
                    [sat, dt_sv] = geph2pos(p,eph_info,prn,tidx,t_sv,'glo',orbit_corr);
                    if ~isnan(sat.pos_ecef(1))
                        satlog.svprn_mark(prn) = p.glo.sys_num;satlog.prn_record(prn) = prn;
                    end
                 case 'gal'
                    [sat, dt_sv] = eph2pos(p,eph_info,prn,tidx,t_sv,'gal',orbit_corr);
                    if ~isnan(sat.pos_ecef(1))
                        satlog.svprn_mark(prn) = p.gal.sys_num;satlog.prn_record(prn) = prn;
                    end
                case 'bds'
                    [sat, dt_sv] = eph2pos(p,eph_info,prn,tidx,t_sv,'bds',orbit_corr); 
                    if ~isnan(sat.pos_ecef(1))
                        satlog.svprn_mark(prn) = p.bds.sys_num;satlog.prn_record(prn) = prn;
                    end
            end
            if ~isnan(sat.pos_ecef(1))
                satlog.tp(prn) = tp_prime+dt_sv;
                satlog.s_pos_ecef(:,prn) = sat.pos_ecef;
                satlog.s_v_ecef(:,prn) = sat.v_ecef;
                if p.post_mode == 1
                    satlog.s_pos_prc(:,prn) = sat.pos_prc;
                end
                % ts = tr - rho/c - dt_sv - dt_corr + code_bias;
                satlog.corr_range(prn) = obs_range(prn)+p.c*dt_sv+clock_corr_m;
                
            end
        end
    end
end
satlog.num_sv = sum(satlog.svprn_mark~=0);
end