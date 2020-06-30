function [F,dFdx] = DasEulerFunc(x)
global xp  h  u HandGoal Khand Bhand Kang Bang handF model;


% Adjust forces on end of wrist using a PID controller
[dPhand_dx, Phand] = pos_jacobian(x,model);

Vhand=dPhand_dx*x(12:22);

handF=-Khand*(Phand-HandGoal)-Bhand;

[F, dfdx, dfdxdot, dfdu, FGH] = das3('Dynamics',x,(x-xp)/h, u, zeros(5,1), zeros(2,1), handF);
if nargout>1
    % Update extra terms in the Jacobian
%     d2Phand_dx2 = das_hessian(x);
%     [dPhand_dx,~] = pos_jacobian(x);
%     [df_dhandF] = handforce_jacobian(x,xdot,step_u,M,exF,handF);
%     dfdx(:,1:11) = dfdx(:,1:11) - df_dhandF*(Khand*dPhand_dx - Bhand * mult(d2Phand_dx2,x(12:22)));
%     dfdx(:,12:22) = dfdx(:,12:22) - df_dhandF*Bhand*dPhand_dx;
    
    dFdx=dfdx+dfdxdot/h;
end
