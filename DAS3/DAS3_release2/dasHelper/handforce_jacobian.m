% calculating hand force jacobians
function [df_dhandF] = handforce_jacobian(x,xdot,u,M,exF,handF)
    df_dhandF = zeros(298,3);
    [f, ~, ~, ~, ~] = das3h(x, xdot, u, M, exF, handF);
    h = 1e-7;
    for i = 1:3
        tmp = handF(i);
        handF(i) = handF(i) + h;
        [f_new, ~, ~, ~, ~] = das3h(x, xdot, u, M, exF, handF);
        df_dhandF(:,i) = (f_new - f)/h;
        handF(i) = tmp;
    end
end
