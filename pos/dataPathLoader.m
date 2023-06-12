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
end

end 