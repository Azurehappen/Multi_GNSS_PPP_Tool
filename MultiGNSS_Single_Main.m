% This code implement single frequency code-based GNSS (GPS, GLO, GAL, BDS)
% PPP solution using real-time PPP correction.
% PPP corrections:
%   SSR satellite orbit and clock correction (WHU stream)
%   Troposphere correction: IGGtrop or UNB3M
%   Ionosphere correction: SSR VTEC (CNES stream)
%   Satellite code bias (OSR): GIPP product
% clear all
% close all
%--------------------------------%
addpath('data')
addpath('parser')
addpath('time_compute')
addpath('eph')
addpath('pos')
addpath('corr')
addpath('init')
addpath('estimation')
%--------------------------------%
% Pick the Data Number
data_num = 2;
files = dataPathLoader(data_num);
%--------------------------------%

% Initialize parameters
if exist(files.preload,'file')==2 % Check if the data already been parsed
    load(files.preload);
else
    [p, eph, obs] = readDataFiles(files);
    p.post_mode  = p.mode_dgnss; %%%% 0=Standard GNSS, 1 = PPP, 2= DGNSS
    %--------------------------------%
    [p, obs] = loadDataAndCorr(p, files, eph, obs);
    save(files.preload, 'p', 'eph', 'obs');
end
%-------------%
p = initialParameters(p, files, eph);
p.run_mode = 0;
p.post_mode  = p.mode_dgnss; %%%% sps=Standard GNSS, ppp = PPP, dgnss = DGNSS
p.IGS_enable = 1;
p.VRS_mode = 0;
p.double_diff = 0;
p.elev_mark  = 10*pi/180;
p.enableGPS  = 1; % Enable GPS: 1 means enable, 0 means close
p.enableGLO  = 0; % Enable GLO: 1 means enable, 0 means close
p.enableGAL  = 1; % Enable GAL: 1 means enable, 0 means close
p.enableBDS  = 1; % Enable BDS: 1 means enable, 0 means close
p.inval = 1; % Computation time interval 
p.tec_tmax = 15;
p.tec_tmin = 0;
p.L2enable = 0;
p.enable_vtec = false;
p.est_mode = p.raps_est;
p.state_mode = p.pva_mode;

if p.state_mode == p.pos_mode
    p.ekf_para.q_pos = 30^2;
end

output = compute_gnss_ecef(p,eph,obs);

%%
figure
scatter(p.t,output.err,'.')
xtickformat('yyyy-MM-dd HH:mm:ss')
title('ECEF positioning error')
xlabel('Local time')
ylabel('Error, unit: meter');grid on

figure
scatter(p.t,output.hor_err,'.')
xtickformat('yyyy-MM-dd HH:mm:ss')
title('Horizontal positioning error')
xlabel('Local time')
ylabel('Error, unit: meter');grid on

total = output.sv_num_GPS + output.sv_num_GAL + output.sv_num_BDS;
figure
scatter(p.t,output.sv_num_GPS,'.')
hold on
scatter(p.t,output.sv_num_GAL,'.')
hold on
scatter(p.t,output.sv_num_BDS,'.')
hold on
scatter(p.t,total,'.')
title('total satellites been used')
legend('GPS','GAL','BDS','Total')
xlabel('Receiver time using GPS second')
ylabel('Distance, unit: meter');grid on
% % 
figure
scatter(p.t,output.rover_clk/p.c,'.')
title('Local bias')
xlabel('Receiver time using GPS second');
ylabel('Clock bias, seconds');grid on

figure
scatter(p.t,output.isb_gal,'.')
title('ISB GAL')
xlabel('Receiver time using GPS second');
ylabel('GPS-GAL ISB, meter');grid on

figure
scatter(p.t,output.isb_bds,'.')
title('ISB BDS')
xlabel('Receiver time using GPS second');
ylabel('GPS-BDS ISB, meter');grid on

figure
subplot(311)
scatter(p.t,sqrt(output.ned_cov(1,:)),'.')
ylabel('Cov N, meter');grid on
subplot(312)
scatter(p.t,sqrt(output.ned_cov(2,:)),'.')
ylabel('Cov E, meter');grid on
subplot(313)
scatter(p.t,sqrt(output.ned_cov(3,:)),'.')
ylabel('Cov D, meter');grid on
title('NED Pos Covariance')
xlabel('Receiver time using GPS second');

figure
subplot(311)
scatter(p.t,sqrt(output.state_cov(1,:)),'.')
ylabel('Cov x, meter');grid on
subplot(312)
scatter(p.t,sqrt(output.state_cov(2,:)),'.')
ylabel('Cov y, meter');grid on
subplot(313)
scatter(p.t,sqrt(output.state_cov(3,:)),'.')
ylabel('Cov z, meter');grid on
title('ECEF Pos Covariance')
xlabel('Receiver time using GPS second');

if p.est_mode == p.raps_est
    figure
    subplot(311)
    scatter(p.t,sqrt(output.pos_info_ned(1,:)),'.')
    hold on
    yline(sqrt(1/p.raps.hor_cov_spec))
    ylabel('Infor N, meter');grid on
    subplot(312)
    scatter(p.t,sqrt(output.pos_info_ned(2,:)),'.')
    hold on
    yline(sqrt(1/p.raps.hor_cov_spec))
    ylabel('Infor E, meter');grid on
    subplot(313)
    scatter(p.t,sqrt(output.pos_info_ned(3,:)),'.')
    hold on
    yline(sqrt(1/p.raps.ver_cov_spec))
    ylabel('Infor D, meter');grid on
    title('NED Pos Information')
    xlabel('Receiver time using GPS second');

    figure
    subplot(311)
    scatter(p.t,sqrt(output.state_info(1,:)),'.')
    hold on
    plot(p.t, sqrt(output.raps_spec_xyz(1,:)))
    ylabel('Infor X, meter');grid on
    subplot(312)
    scatter(p.t,sqrt(output.state_info(2,:)),'.')
    hold on
    plot(p.t, sqrt(output.raps_spec_xyz(2,:)))
    ylabel('Infor Y, meter');grid on
    subplot(313)
    scatter(p.t,sqrt(output.state_info(3,:)),'.')
    hold on
    plot(p.t, sqrt(output.raps_spec_xyz(3,:)))
    ylabel('Infor Z, meter');grid on
    title('ECEF Pos Information')
    xlabel('Receiver time using GPS second');
end

percentage = (sum(output.hor_err < 1.0) / length(output.hor_err)) * 100
percentage = (sum(output.hor_err < 1.5) / length(output.hor_err)) * 100
percentage = (sum(output.err < 3.0) / length(output.hor_err)) * 100

plotEstPosOnMap(output.pos_ecef, output.hor_err)
% figure
% subplot(311)
% scatter(p.t,output.ned_err(1,:),'.')
% title('North Error in NED');grid on;
% subplot(312)
% scatter(p.t,output.ned_err(2,:),'.')
% title('East Error in NED');grid on;
% subplot(313)
% scatter(p.t,output.ned_err(3,:),'.')
% title('Down Error in NED')
% xlabel('Receiver time using GPS second');grid on;
% figure
% for i=1:32
%     scatter(output.gpst,output.res_GPS(i,:),'.')
%     hold on
% end
% hold off
% grid on
% title('GPS residual')
% xlabel('Receiver time using GPS second');
% ylabel('Residual, unit: meter');
% figure
% scatter(output.gpst,output.sv_num_GPS,'.')
% title('Number of GPS satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on
% 
% figure
% scatter(output.gpst,output.sv_num_GLO,'.')
% title('Number of GLO satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on
% 
% figure
% scatter(output.gpst,output.sv_num_GAL,'.')
% title('Number of GAL satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on
% 
% figure
% scatter(output.gpst,output.sv_num_BDS,'.')
% title('Number of BDS satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on