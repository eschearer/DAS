% new function for arm force
function [f, dfdx, dfdxdot, dfdu, FGH] = das3h(x, xdot, u, M, exF, handF)
    [f, dfdx, dfdxdot, dfdu, FGH] = das3('Dynamics',x , xdot, u, M, exF, handF);
end