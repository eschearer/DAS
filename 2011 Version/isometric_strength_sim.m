function [state,moments,fval,exitflag,output] = isometric_strength_sim(idof, sign_m, posture, params, x_prev)
% Finds the maximum moment for a given posture of the arm model,
% ensuring glenohumeral and scapular stability
%
% das3mex used here has one extra input: the five moments
% and one extra output: the thoraco-humeral Y-Z-Y angles
%
% Inputs:
% (1) which moment idof (1-5) we want to maximize
% (2) sign of the moment (1: maximum positive, -1: maximum negative)
% (3) the posture we want (3 thorax-humerus angles, 2 elbow angles)
%
% Outputs:
% (1) state
% (2) five moments (for the 3 thorax-humerus angles, 2 elbow angles)
%
% Dimitra Blana, March 2012

% input checks
if idof>5
    disp('idof should be between 1 and 5');
    state=[];
    moments=[];
    fval=0;
    exitflag=0;
    return;
end

if sign_m ~=1 && sign_m ~=-1
    disp('The sign of the moment has to be either 1 (positive) or -1 (negative)');
    state=[];
    moments=[];
    fval=0;
    exitflag=0;
    return;
end

if length(posture)~=5
    disp('posture should have five elements: 3 thorax-humerus angles and 2 elbow angles.');
    state=[];
    moments=[];
    fval=0;
    exitflag=0;
    return;
end

% initialize model
das3mex();

% model sizes
ndof = 11;
nmus = 138;
nstates = 2*ndof + 2*nmus;

% indices to the state variables within the state vector
iq = 1:ndof;
iqdot = max(iq) + (1:ndof);
iLce = max(iqdot) + (1:nmus);
iAct = max(iLce) + (1:nmus);

nvar = nstates+5;					
ncon = 2*ndof+nmus+6;			    
                                        
% range of motion limits
xlimdeg=  [-65		-20;
            5		30;
            0  		83; 
            33 		69;
            -18 	7;
            -17		18;
            -90	    90; 
            -30    	90;
            -90   	90;    
            5      	140;
            5     	160];
        
xlims =  xlimdeg*pi/180;

if nargin>4
    % initial guess: previous solution
    x0 = x_prev;
else
    % initial guess for state: equilibrium
    % initial guess for moments: zero
     xeq = load('hang_equilibrium');
     x0 = [xeq.state; 0; 0; 0; 0; 0];
end

LB = -Inf*ones(nvar,1);
UB = Inf*ones(nvar,1);
% muscle active states between 0 and 1
LB(iAct) = 0;
UB(iAct) = 1;

% Lce
LB(iLce) = -5;
UB(iLce) = 5;
% angles
LB(iq) = xlims(:,1);
UB(iq) = xlims(:,2);
% elbow angles = required posture
LB(10:11) = posture(4:5);
UB(10:11) = posture(4:5);

% velocities are zero, since position is static
LB(iqdot) = 0;
UB(iqdot) = 0;

% the moment we are trying to maximize/minimize should have the right sign!
if sign_m>0
    LB(nstates+idof)=0;
else
    UB(nstates+idof)=0;
end  

% for glenohumeral constraint:
Rgt = glenoid_scap;

gceq_nnz = 1;
xrand = LB + (UB-LB).*rand(size(LB));		% a random vector of unknowns
J = conjac(xrand);
gceq_nnz = nnz(J);
Jpattern = double(J~=0);

if strcmpi(params.solver, 'fmincon')
    % find x using fmincon
    % to check derivatives, set 'DerivativeCheck' to 'on'
    options = optimset('Display','iter','Algorithm',params.alg,...
        'GradConstr','on','GradObj','on','DerivativeCheck','off',...
        'MaxFunEvals',params.maxfulevals,'TolX',params.TolX,...
        'TolCon',params.TolCon,'TolFun',params.TolFun);
    [x,fval,exitflag,output] = fmincon(@myfun,x0,[],[],[],[],LB,UB,@mycon,options);
%     [x,fval,exitflag,output] = fmincon(@(x)0,x0,[],[],[],[],LB,UB,@mycon,options);
else
    
    % find x using IPOPT
    % to check derivatives, set checkderivatives to 1
    checkderivatives=0;
    if (checkderivatives)
        hh = 1e-7;
        x = x0;
        obf = objfun(x);
        gradf = objgrad(x);
        hess = objhess(x);
        cf = confun(x);
        cjac = conjac(x);
        cjac_num = zeros(ncon,nvar);
        grad_num = zeros(nvar,1);
        hess_num = zeros(nvar,nvar);
        for i=1:nvar
            fprintf('checking derivatives for unknown %4d of %4d\n',i,nvar);
            xisave = x(i);
            x(i) = x(i) + hh;
            helpv = confun(x);
            cjac_num(:,i) = (confun(x) - cf)/hh;
            grad_num(i) =   (objfun(x) - obf)/hh;
            hess_num(:,i) = (objgrad(x) - gradf)/hh;
            x(i) = xisave;
        end

        % find the max difference in constraint jacobian and objective gradient
        [maxerr,irow] = max(abs(cjac-cjac_num));
        diffcjac = abs(cjac-cjac_num);
        [maxerr,icol] = max(maxerr);
        fprintf('Max.error in constraint jacobian: %8.5f at %d %d\n', maxerr, irow(icol), icol);
        [maxerr,irow] = max(abs(grad_num-gradf));
        fprintf('Max.error in objective gradient:  %8.5f at %d\n', maxerr, irow);
        [maxerr,irow] = max(abs(hess-hess_num));
        [maxerr,icol] = max(maxerr);
        fprintf('Max.error in objective Hessian:   %8.5f at %d %d\n', maxerr, irow(icol), icol);
        keyboard
    end

    fval=0;
    output=0;
    % tell Matlab where IPOPT is
    addpath('tools');

    funcs.objective = @objfun;
    funcs.gradient  = @objgrad;
    funcs.constraints = @confun;
    funcs.jacobian    = @conjac;
    funcs.jacobianstructure = @conjacstructure;
    options.lb = LB;
    options.ub = UB;
    %             equality constraints     -1<GH force<0         TS and AI<0         
    options.cl = [zeros(2*ndof+nmus+3,1);        -1;             -Inf; -Inf];
    options.cu = [zeros(2*ndof+nmus+3,1);         0;                0;    0];
    options.IPOPT.max_iter = 1000; %MaxIterations
    options.IPOPT.hessian_approximation = 'limited-memory';
    options.IPOPT.mu_strategy = 'adaptive';		% worked better than 'monotone'
    options.IPOPT.bound_frac = 0.001;			% worked better than 0.01 or 0.0001
    options.IPOPT.bound_push = options.IPOPT.bound_frac;
    options.IPOPT.tol = 1e-4; %OptimalityTolerance
    options.IPOPT.limited_memory_max_history = 12;	% 6 is default, 12 converges better, but may cause "insufficient memory" error when N is large
        
    [x, info] = ipopt(x0,funcs,options);
    exitflag = info.status;
    fprintf('IPOPT info status: %d\n', exitflag);
end

% function outputs:
state = x(1:nstates);
moments = x(nstates+1:end);

    % objective function (and gradient) for fmincon
    function [f g] = myfun(x)
        % we want to maximize moment(idof) = x(nstates+idof)
        f = -sign_m*(x(nstates+idof));

        if nargout > 1
            g = objgrad(x);
        end
    end

    % objective function for IPOPT
    function f = objfun(x)
        f = myfun(x);
    end

    % gradient of objective function for IPOPT
    function g = objgrad(x)
        g = zeros(nvar,1);
        g(nstates+idof)=-sign_m;
        
    end

    % hessian of objective function
    function h = objhess(x)
        % hessian of objective function
        hdiag = zeros(nvar,1);        
        
        % store hessian as a diagonal sparse matrix
        h = spdiags(hdiag,0,nvar,nvar);
    end

    % constraints (and gradients) for fmincon
    function [c,ceq,gc,gceq] = mycon(x)
        % Non-linear constraints: c(x) <= 0, ceq(x) = 0
        
        % das3mex has one extra input: the five moments (x(nstates+1:end))
        % and one extra output: the thoraco-humeral Y-Z-Y cardan angles (thorhum)
        [f, ~, ~, ~, FGH, ~, thorhum] = das3mex(x(1:nstates),zeros(nstates,1),zeros(nmus,1),x(nstates+1:end));
        
        % calculate GH constraint
        Fgh0 = Rgt*FGH;  % take glenoid orientation into account
        if norm(Fgh0), Fgh0 = Fgh0/norm(Fgh0); end
        % decompose into polar angles
        thetar = asin(-Fgh0(2));
        
        if ~(sqrt(Fgh0(1)^2+Fgh0(3)^2)), phir = 0.0;
        else phir=asin(Fgh0(3)/sqrt(Fgh0(1)^2+Fgh0(3)^2));
        end
        FGHcontact = (thetar/0.7744)^2 + (phir/0.6728)^2 - 1; % <=0
        
        % solve thorax ellipsoid surface equation for TS and AI, to find
        % out whether they are inside, on, or outside the thorax
        Fscap = das3mex('Scapulacontact', x(1:nstates))';  %<=0
        
        % Inequalities:
        % FGHcontact within elliptical constraint (between -1 and 0)
        % TS on or inside the ellipsoid surface
        % AI on or inside the ellipsoid surface
        c = [FGHcontact Fscap];
        
        % Equalities:
        % f*(x,0,0,M) = 0
        % with f* being only the elements of f without activation dynamics
        
        % in f:
        % the first NDOF rows are: qdot-dq/dt = 0
        % the next NDOF rows are the equations of motion from Autolev
        % the next NMUS rows are the muscle contraction dynamics
        % the final NMUS rows are the muscle activation dynamics: da/dt - (u-a)(c1*u + c2) = 0
        
        fstar = f(1:2*ndof+nmus)';
        
        % thoraco-humeral angles = required posture
        fthorhum = thorhum' - posture(1:3);
        
        ceq = [fstar fthorhum];
        
        % Gradients
        if nargout > 2
            [~,gc,gceq] = congrad(x,FGHcontact,Fscap,fthorhum);
        end
    end

    % constraints for IPOPT
    function allcon = confun(x)
        [c,ceq] = mycon(x);
        allcon = [ceq c];      
    end

    function J = conjacstructure()
    % returns structure of constraint Jacobian matrix
        J = Jpattern;
    end

    % Jacobian and gradients of constraints
    function [J,gc,gceq] = congrad(x,FGHcontact,Fscap,fthorhum)
        
        [f, dfdx, ~, ~,FGH, ~, thorhum] = das3mex(x(1:nstates),zeros(nstates,1),zeros(nmus,1),x(nstates+1:end));
        
        if nargin<4
            % calculate GH, scapula and posture constraints
            Fgh0 = Rgt*FGH;  
            if norm(Fgh0), Fgh0 = Fgh0/norm(Fgh0); end
            thetar = asin(-Fgh0(2));        
            if ~(sqrt(Fgh0(1)^2+Fgh0(3)^2)), phir = 0.0;
            else phir=asin(Fgh0(3)/sqrt(Fgh0(1)^2+Fgh0(3)^2));
            end
            FGHcontact = (thetar/0.7744)^2 + (phir/0.6728)^2 - 1; % <=0

            Fscap = das3mex('Scapulacontact', x(1:nstates))';  %<=0
            fthorhum = thorhum' - posture(1:3);
        end
        
        h = 1e-7;
        dFGHdx = zeros(nvar,1);
        dFscapdx = zeros(nvar,2);
        dthorhum = zeros(nvar,3);

        %	dfdx	(298 x 298 sparse) 	Jacobian of f with respect to x
        dfdxstar = dfdx(1:2*ndof+nmus,:);
        % Gradient matrix: transpose of the form of Jacobians
        % (one column per constraint)
        gceq = spalloc(nvar, 2*ndof+nmus+3, gceq_nnz);
        gceq(1:nstates,1:2*ndof+nmus) = dfdxstar';

        for ix=1:nvar
            saved = x(ix);
            x(ix) = x(ix)+h;
            [fnew, ~, ~, ~, FGH, ~, thorhum] = das3mex(x(1:nstates),zeros(nstates,1),zeros(nmus,1),x(nstates+1:end));
            Fgh0 = Rgt*FGH;
            if norm(Fgh0), Fgh0 = Fgh0/norm(Fgh0); end
            thetar = asin(-Fgh0(2));
            if ~(sqrt(Fgh0(1)^2+Fgh0(3)^2)), phir = 0.0;
            else phir=asin(Fgh0(3)/sqrt(Fgh0(1)^2+Fgh0(3)^2));
            end
            FGHcontact_new = (thetar/0.7744)^2 + (phir/0.6728)^2 - 1; % <=0
            dFGHdx(ix) = (FGHcontact_new-FGHcontact)/h;

            Fscap_new = das3mex('Scapulacontact', x(1:nstates))';
            dFscapdx(ix,:) = (Fscap_new-Fscap)/h;

            fthorhum_new = thorhum' - posture(1:3);
            dthorhum(ix,:) = (fthorhum_new-fthorhum)/h;
            if ix>nstates
                gceq(ix,1:2*ndof+nmus) = (fnew(1:2*ndof+nmus)-f(1:2*ndof+nmus))/h; 
            end
            x(ix) = saved;
        end
        gc = [dFGHdx dFscapdx];
        gceq(:,2*ndof+nmus+1:2*ndof+nmus+3) = dthorhum;  
        
        J = [gceq gc]';
    end

    % Jacobian of constraints for IPOPT
    function J = conjac(x)
        [J,~,~] = congrad(x);
    end
end