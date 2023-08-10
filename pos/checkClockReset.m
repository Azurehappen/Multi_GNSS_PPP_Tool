function [x, cov, flag] = checkClockReset(p, x0, cov0, clk_ind, res, cpt)
x = x0;
cov = cov0;
flag = false;
if sum(abs(res) > 5000) > 0.8*length(res)
    [estState0,~] = userpos(p,cpt);
    x(clk_ind) = estState0.clock_bias;
    cov(:, clk_ind:clk_ind+1) = zeros(length(x0),2);
    cov(clk_ind:clk_ind+1,:) = zeros(2,length(x0));
    cov(clk_ind, clk_ind) = 100^2;
    cov(clk_ind+1, clk_ind+1) = 100^2;
    flag = true;
end
end