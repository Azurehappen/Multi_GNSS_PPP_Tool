function [cpt,n] = diff_corr_compute(p,cpt,tr_posix)
% Compute differential correction
% Input:
%       p       -- parameters which include base station's eph and obs
%       cpt     -- Satellite computation results
%       obs_tr  -- The time that rover received the data
diff_corr = NaN(length(cpt.corr_range),1);
ind_sv = find(cpt.svprn_mark~=0);
n = find(p.obs_b.tr_posix>(tr_posix-180) & p.obs_b.tr_posix<tr_posix);
if isempty(n)
    return;
end

delta = abs(p.obs_b.tr_sow(n)-tr_posix);
tridx = find(delta == min(delta));
tridx = n(tridx(end));
obstr_gps = p.obs_b.tr_posix(tridx);
obs_tr = obstr_gps;
for i = 1:length(cpt.corr_range)
    prn = cpt.prn_record(ind_sv(i));
    sys_type = cpt.svprn_mark(ind_sv(i));
    switch sys_type
        case p.gps.sys_num
            pseR = p.obs_b.gps(1).data.P(prn,tridx);
            Strength = p.obs_b.gps(1).data.S(prn,tridx);
            message_duration = p.gps.message_duration;
            eph_info = p.eph_b.gps;
        case p.glo.sys_num
            pseR = p.obs_b.glo(1).data.P(prn,tridx);
            Strength = p.obs_b.glo(1).data.S(prn,tridx);
            message_duration = p.glo.message_duration;
            eph_info = p.eph_b.glo;
            obs_tr = obstr_gps - p.glo.lps_gps;% Correct time diff from GPS time to GLO time
        case p.gal.sys_num
            pseR = p.obs_b.gal(1).data.P(prn,tridx);
            Strength = p.obs_b.gal(1).data.S(prn,tridx);
            message_duration = p.gal.message_duration;
            eph_info = p.eph_b.gal;
            obs_tr = obstr_gps - p.gal.lps_gps;% Correct time diff from GPS time to GAL time
        case p.bds.sys_num
            pseR = p.obs_b.bds(1).data.P(prn,tridx);
            Strength = p.obs_b.bds(1).data.S(prn,tridx);
            message_duration = p.bds.message_duration;
            eph_info = p.eph_b.bds;
            obs_tr = obstr_gps - p.bds.lps_gps;% Correct time diff from GPS time to BDS time
    end
    if pseR == 0 || isnan(pseR) || prn>size(eph_info.SV_health,1)
        cpt.num_sv(sys_type) = cpt.num_sv(sys_type)-1;
        cpt.svprn_mark(ind_sv(i)) = 0;
        cpt.prn_record(ind_sv(i)) = 0;
        continue;
    end
    tp_prime = pseR/p.c;
    t_sv = obs_tr - tp_prime;
    %---------------------------------% Find the time index in eph data
    teph = ephtidx(eph_info.t_oc{prn},t_sv,eph_info.SV_health(prn,:),message_duration);
    %---------------------------------%
    % Check the signal strength and sv health (health message 0 means ok)
    if isempty(teph) || Strength < p.sig_strg
        cpt.num_sv(sys_type) = cpt.num_sv(sys_type)-1;
        cpt.svprn_mark(ind_sv(i)) = 0;
        cpt.prn_record(ind_sv(i)) = 0;
        continue;
    end
    switch(sys_type)
        % compute ephemeris to satellite position and clock bias
        case p.gps.sys_num
            [sat, dt_sv] = eph2pos(p,eph_info,prn,teph,t_sv,'gps',[]);
        case p.glo.sys_num
            [sat, dt_sv] = geph2pos(p,eph_info,prn,teph,t_sv);
        case p.gal.sys_num
            [sat, dt_sv] = eph2pos(p,eph_info,prn,teph,t_sv,'gal',[]);
        case p.bds.sys_num
            [sat, dt_sv] = eph2pos(p,eph_info,prn,teph,t_sv,'bds',[]);
    end
    if isnan(sat.pos_ecef(1))
        continue;
    end
    range = norm(p.P_base-sat.pos_ecef)+sagnac(p,sat.pos_ecef,p.P_base);
    diff_corr(i) = pseR + p.c*dt_sv - range; % The correctoin here include receiver clock bias.
end
% Get the field names
fields = fieldnames(cpt);
% Find indices to delete
del_ind = find(isnan(diff_corr));
% Loop over each field
for i = 1:length(fields)
    % Get the field name
    field = fields{i};
    if strcmp(field,'num_sv')
        continue;
    end
    % Check the size of the field
    if size(cpt.(field), 1) == 3
        % If the field is pos, vel...
        cpt.(field)(:, del_ind) = [];
    else
        cpt.(field)(del_ind, :) = [];
    end
end
% Update diff_corr
diff_corr(del_ind) = [];
cpt.diff_corr = diff_corr;
end