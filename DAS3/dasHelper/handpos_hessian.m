% calculating hand position hessians
function [d2Phand_dx2, dPhand_dx] = handpos_hessian(x)
    [dPhand_dx, ~, ~] = handpos_jacobian(x);
    d2Phand_dx2 = zeros(3,11,11);
    h = 1e-7;
    for i = 1:11
        tmp = x(i);
        x(i) = x(i) + h;
        dPhand_dx_new = handpos_jacobian(x);
        d2Phand_dx2(:,:,i) = (dPhand_dx_new - dPhand_dx)/h;
        x(i) = tmp;
    end
end