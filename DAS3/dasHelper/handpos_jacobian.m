% calculating hand position jacobians
function [dPhand_dx, Phand, Phand_new] = handpos_jacobian(x)
    dPhand_dx = zeros(3,11);
    stick = das3('Visualization', x);
    Phand = stick(13,1:3)';
    h = 1e-7;
    for i = 1:11
        tmp = x(i);
        x(i) = x(i) + h;
        stick = das3('Visualization',x);
        Phand_new = stick(13,1:3)';
        dPhand_dx(:,i) = (Phand_new - Phand)/h;
        x(i) = tmp;
    end        
end