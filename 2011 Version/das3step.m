function [xnew, FGH] = das3step(x, u, h, M, exF)
	% Ton van den Bogert
	% (c) 2010-2011 Case Western Reserve University

	% For dynamic simulation.  Advances the model by one time step.

	% Inputs
	%	x		(298 x 1) System state
	%	u		(138 x 1) Muscle excitations, assumed constant during the time step
	%	h		(scalar) Time step size
    % Optional inputs
    % 	M		(5 x 1)		Moments applied to the thorax-humerus YZY axes and the elbow flexion and supination axes
	%   exF     (2 x 1)     Vertical force of amplitude exF(2) applied to the ulna at a distance of exF(1) from the elbow 
	%						(to simulate a mobile arm support)

	% Outputs
	%	xnew	(298 x 1) The system state at the end of the time step
	%	FGH		(3 x 1) The glenohumeral reaction force, acting on the scapula

	% Method: First order Rosenbrock formula on implicit differential equation
	
	% Since this algorithm is just a few lines of Matlab code, we may code it inline inside a
	% Simulink block to speed things up, but I am not familiar enough with Simulink to give
	% specific recommendations.
	
	global das3step_xdot das3step_u		% global so we can preserve this value between function calls
	
	% If this is the very first time this function is called, initialize derivatives and
	% muscle excitations to zero.
	% Consequently, startup will work best when initial state x is a passive equilibrium
	% state which satisfies the implicit differential equation when xdot and u are both zero.
	if ~(exist('das3step_xdot')==1) || isempty(das3step_xdot)
		das3step_xdot = zeros(size(x));
	end
	if ~(exist('das3step_u')==1) || isempty(das3step_u)
		das3step_u = zeros(size(u));
	end

	% Evaluate dynamics in current x and xdot
	if nargin>4
        [f, dfdx, dfdxdot, dfdu, FGH] = das3mex( x , das3step_xdot, das3step_u, M, exF);	
    elseif nargin>3
        [f, dfdx, dfdxdot, dfdu, FGH] = das3mex( x , das3step_xdot, das3step_u, M);
    else
        [f, dfdx, dfdxdot, dfdu, FGH] = das3mex( x , das3step_xdot, das3step_u);
    end
    
	% Solve the change in x from the 1st order Rosenbrock formula
	du = u - das3step_u;
    
    easy = 1;       % easy=0 uses partitioning, used in the (possibly faster) das3step C code
                    % Here, in Matlab, it seems to make no difference in performance
    if (easy)
        dx = (dfdx + dfdxdot/h)\(dfdxdot*das3step_xdot - f - dfdu*du);
        % the line above is from the IUTAM paper, but we can take some advantage
        % of the structure of the matrix dfdx+dfdxdot/h and solve a smaller linear system
    else
        A = dfdx+dfdxdot/h; 
        b = dfdxdot*das3step_xdot - f - dfdu*du;
        % we need to solve A*dx = b, using the fact that A has this
        % structure:  [     A1(160x160) A2(160x138)    ; 
        %                zeros(138x160) diag(D(138x1))  ]
        % first solve dx(161:end) from the lower half:
        dx = zeros(size(x));
        dx(161:end) = b(161:end)./ diag(A(161:end,161:end));
        % now solve dx(1:160) from the upper half:
        dx(1:160)   = A(1:160,1:160) \ (b(1:160) - A(1:160,161:end)*dx(161:end));
    end

    % update x
	xnew = x + dx;
	
	% update of global variables for the next simulation step
	das3step_xdot = dx/h;
	das3step_u = u;

end
