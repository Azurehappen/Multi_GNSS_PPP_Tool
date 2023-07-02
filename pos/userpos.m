function [estState,res] = userpos(p,cpt)
% This is solver for computing user position, receiver clock
% bias, satellite system offset.
% Measurement selection applied
% Input: 
%       s_pos_ecef: 3-by-N Satellite position in ECEF frame.
%       x0 : 3-by-1 initial interative coordinates in ECEF frame.
%       y: m-by-1 Corrected pseudorange.
%
% Output:
%       
%       
%       
%       

%-------------------%
% Initialize
estState.isb_dict = dictionary;
estState.isb_dict(p.glo.sys_num) = NaN;
estState.isb_dict(p.gal.sys_num) = NaN;
estState.isb_dict(p.bds.sys_num) = NaN;

x0 = p.state0(1:4);
[H_isb,x_isb] = formIsbStatesAndH(cpt.num_sv);
xk = [x0;x_isb];
%------------------%
[estState.pos,estState.clock_bias,isb_est,res] = LSsolver(p,xk,H_isb,cpt);
j = 1;
for i = 2:length(cpt.num_sv)
    if cpt.num_sv(i) ~= 0
        estState.isb_dict(i) = isb_est(j);
        j = j+1;
    end
end
% switch p.select
%     case 0
%         [pos,clock_bias,res] = LSsolver(p,xk,H_offset,s_pos_ecef,y);
%     case 1
%         [pos,clock_bias,res,cost] = TDsolver(p,xk,H_offset,s_pos_ecef,y);
%     case 2
%         [pos,clock_bias,res,cost] = LSSsolver(p,xk,H_offset,s_pos_ecef,y);
%     case 3
%         [pos,clock_bias,res,cost] = MSsolver(p,xk,H_offset,s_pos_ecef,y);
%     case 4
%         [pos,clock_bias,res,cost] = LTSsolver(p,xk,H_offset,s_pos_ecef,y);
% end

end
