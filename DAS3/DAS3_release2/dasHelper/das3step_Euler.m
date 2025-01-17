% das3step.m function
function [xnew] = das3step_Euler(x, u, h, xdot, step_u, M, exF, handF, K, B)
	% Ton van den Bogert
	% (c) 2010-2011 Case Western Reserve University
    %
	% For dynamic simulation.  Advances the model by one time step.
    %
	% Inputs
	%	x		(298 x 1) System state
	%	u		(138 x 1) Muscle excitations, assumed constant during the time step
	%	h		(scalar)  Time step size
    %   xdot    (298 x 1) System state derivatives (initially zero)
	%	step_u	(138 x 1) Muscle excitations from previous step (initially zero)
    % Optional inputs
    % 	M		(5 x 1)	  Moments applied to the thorax-humerus YZY axes and the elbow flexion and supination axes
	%   exF     (2 x 1)   Vertical force of amplitude exF(2) applied to the ulna at a distance of exF(1) from the elbow 
	%						(to simulate a mobile arm support)
    %   handF   (3 x 1)   Force at the CoM of the hand
    %
	% Outputs
	%	xnew	(298 x 1) The system state at the end of the time step
	%	FGH		(3 x 1)   The glenohumeral reaction force, acting on the scapula
    %
	% Method: First order Rosenbrock formula on implicit differential equation	
	
	% Evaluate dynamics in current x and xdot
    global xp
    
    x0=xp;
    xnew=fsolve(@DasEulerFunc,x0);

end