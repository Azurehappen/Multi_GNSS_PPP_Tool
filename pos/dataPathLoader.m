function files = dataPathLoader(data_num)

switch(data_num)
    case 1
        files.eph = 'data/PPP/BRDC_20231440.rnx';
        files.obs = 'data/PPP/ppp_raps_test_1.obs';
        files.ssr = 'data/PPP/SSRA00CNE01440.23C';
        files.vtec = 'data/PPP/SSRA00CNE01440.23C';
        files.code_bias = 'data/PPP/CAS0MGXRAP_20231420_OSB.BIA';
        files.ustec_data = 'data/PPP/ustec_data/';
        files.Grdpos.pos = [-2430697.636;-4704189.056;3544329.070];
        files.Grdpos.t = NaN;
        files.preload = 'data/PPP/preload.mat';
    case 2
        files.eph = 'data/DGNSS_Moving/COM3_200712.nav';
        files.obs = 'data/DGNSS_Moving/COM3_200712.obs';
        files.ssr = [];
        files.vtec = [];
        files.code_bias = [];
        files.ustec_data = [];
        files.data_base = 'data/DGNSS_Moving/rbst0711base.obs';
        files.base_pos = [-2430697.636;-4704189.056;3544329.070];
        M  = readmatrix('data/DGNSS_Moving/groundtruth.csv');
        files.Grdpos.pos = M(:,3:5)';
        files.Grdpos.t = M(:,2);
        files.preload = 'data/DGNSS_Moving/preload.mat';
end

end 