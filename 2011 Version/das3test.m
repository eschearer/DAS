function [out1] = das3test(command, arg1)

% This program runs various tests on the das3 model

	clear global	% to ensure that das3step is properly initialized
	tic;			% set the timer

	if (nargin < 1)
		error('das3test needs command string as first input');
	end
	
	% Some model related variables
	ndof = 11;
	nmus = 138;
	nstates = 2*ndof + 2*nmus;
	
	% define indices to the state variables within the state vector x
	iq = 1:ndof;
	iqdot = max(iq) + (1:ndof);
	iLce = max(iqdot) + (1:nmus);
	iAct = max(iLce) + (1:nmus);
	
	% define DOF names (ensure their order remains consistent with q[] array defined in das3.al)
	dofnames = { 'SCy' 'SCz' 'SCx' 'ACy' 'ACz' 'ACx' 'GHy' 'GHz', 'GHyy' 'ELx', 'PSy'};
	if (size(dofnames,2) ~= ndof)
		error('dofnames has incorrect size');
	end

	% run make.m to make sure we are using latest version of MEX function
	clear mex
	make
	
	% Initialize the model
	das3mex();
		
	% construct a state that should be close to static equilibrium
	xneutral= zeros(nstates,1);
	% joint angles such that the arm hangs down    
    xneutral(1:11) = [-33.486 20.145 32.914 45.571 0.458 -12.062 43.864 21.042 -34.659 5.000 5.000]*pi/180;
	xneutral(iLce) = 2.0;				% all Lce very long so we don't get much passive force
	
	% these are the range of motion limits in das3.bio (useful to know for some of the tests)
	xlimdeg=  [-65		-19;
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
 
    % ======================= do the stick figure test
	if strcmp(command, 'stick')
		disp('Stick figure test...');
		figure(1);clf;
		das3stick(xneutral);
	end
	
	% ======================= do the speed test
	if strcmp(command, 'speed')
		disp('Speed test (takes about 5 seconds)...');
		Neval = 1000;
        
        tic
		for i=1:Neval
			x = rand(nstates,1);
			xdot = rand(nstates,1);
			u = rand(nmus,1);
			M = rand(5,1);
            exF = rand(2,1);	
			[f, dfdx, dfdxdot, dfdu, Fgh, Scap, qTH] = das3mex(x,xdot,u,M,exF);
		end
		fprintf('Computation time for dynamics with Jacobians:    %6.2f ms\n',1000*toc/Neval);
        
        tic
        for i=1:Neval
			x = rand(nstates,1);
			xdot = rand(nstates,1);
			u = rand(nmus,1);
			M = rand(5,1);
            exF = rand(2,1);	
			f = das3mex(x,xdot,u,M,exF);
		end
		fprintf('Computation time for dynamics without Jacobians: %6.2f ms\n',1000*toc/Neval);	
        
		tic
		for i=1:Neval
			x = rand(nstates,1);
			das3mex('Jointmoments',x);
		end
		fprintf('Computation time without Autolev code:           %6.2f ms\n',1000*toc/Neval);		
	end
	
	% ======================= do the moment arms test
	if strcmp(command, 'moment arms')
		disp('Moment arm test...');
		x = zeros(nstates,1);
		x(1:11) = (xlims(:,1)+xlims(:,2))/2;		% angles at midpoint of range of motion
		[MA] = das3mex('Momentarms', x);
		figure(1);clf;
		spy(MA');
		axis('fill');
		xlabel('muscle element');
		ylabel('DOF number');
		
		% for each DOF, report the largest positive and negative moment arms
		fprintf('DOF   Largest positive moment arm          Largest negative moment arm\n')
		fprintf('----  -----------------------------        -----------------------------------\n');
		MA = full(MA);
		for i = 1:ndof
			[dmin, imin] = min(MA(:,i));
			[dmax, imax] = max(MA(:,i));
			name_min = das3mex('Musclename',imin);
			name_max = das3mex('Musclename',imax);
			fprintf('%-4s  %6.1f mm (%-16s)         %6.1f mm (%-16s)\n',dofnames{i},1000*dmax,name_max,1000*dmin,name_min);
		end
	end
	
	% ======================= find the passive equilibrium
	if strcmp(command, 'equilibrium')	
		disp('Equilibrium test...');
		
		% no muscle excitations or speed
		u = zeros(nmus,1);
		xdot = zeros(nstates,1);
		
		% initial guess for equilibrium state
		x = xneutral;      % this initial guess works for sure, but rand(298,1) also works sometimes
		
		% solve the equilibrium with fsolve
		tic;
		figure(1);clf;
        options = optimoptions('fsolve','SpecifyObjectiveGradient', true, 'Display','none', ...
            'FunctionTolerance', 1e-12, 'OptimalityTolerance', 1e-12, ...
            'StepTolerance', 1e-12);
        neval = 0;
        [x,fval,exitflag,output] = fsolve(@resfun, x, options);
		if (exitflag < 1)
			disp('Warning: equilibrium not found');
            disp('type dbcont to continue anyway');
			disp('type dbquit if you do not want this result saved.');
			keyboard
        else
            fprintf('Equilibrium was found in %d iterations\n', output.iterations);
		end

		disp(['Norm of f was: ',num2str(norm(fval),'%10.4e')]);
		disp('This must be close to zero, otherwise it is not an equilibrium');
		das3stick(x);
		disp('Equilibrium posture is shown in Figure 1');
		
		% save the result
		make_osimm('equilibrium',x(1:11));
		save('equilibrium.txt','x','-ascii','-double');
		disp('Equilibrium state was stored on equilibrium.txt');
		disp('Equilibrium posture was stored on equilibrium.sto');

		disp('');
		disp('Hit ENTER to examine the equilibrium state');
		disp('');
		pause
		disp('joint angles (deg):')
		disp(num2str(180/pi*x(1:ndof)','%8.2f'));
		disp('angular velocities (rad/s):')
		disp(num2str(x(ndof+1:2*ndof)','%8.4f'));
		disp('Lce/Lceopt:')
		Lcerel = reshape([x(2*ndof+(1:nmus))' 0 0],10,14)';
		disp(num2str(Lcerel,'%8.4f'));
		disp('Activation states:')
		Act = reshape([x(2*ndof+nmus+(1:nmus))' 0 0],10,14)';
		disp(num2str(Act,'%8.4f'));
		
		disp('');
		disp('Hit ENTER to see other variables');
		disp('');
		pause	

		[f, dfdx, dfdxdot, dfdu, Fgh, Fscap] = das3mex(x, zeros(nstates,1), zeros(nmus,1));
		disp('Force in GH joint:')
		disp(Fgh');
		disp('Contact force in TS:');
		disp(Fscap(:,1)');
		disp('Contact force in AI:');
		disp(Fscap(:,2)');

		[moments] = das3mex('Jointmoments',x);
		disp('Joint moments:');
		disp(moments');

		[forces] = das3mex('Muscleforces',x);
		disp('Muscle forces:');
		forces = reshape([forces' 0 0],10,14)';
		disp(num2str(forces,'%8.4f'));
		
		disp('');
		disp('Hit ENTER to see eigenvalue analysis');
		disp('');
		pause	
		eigenvalues = eig(full(-inv(dfdxdot)*dfdx));
		figure(2);clf;semilogx(eigenvalues,'o');title('eigenvalues of equilibrium state');
		xlabel('Real');
		ylabel('Imaginary');
		disp('Eigenvalues with largest real part (s^-1):');
		disp('(must all be negative for stability)');
		realparts = sort(real(eigenvalues),'descend');
		disp(realparts(1:5)');
		disp('');
		disp('Smallest time constants (s):');
		timeconst = sort(1./abs(eigenvalues));
		disp(num2str(timeconst(1:5)','%12.4e'));
	end
	
	function [f,J] = resfun(x)
	% returns f and df/dx for equilibrium solver
		neval = neval + 1;
		[f, J, dfdxdot, dfdu] = das3mex(x,xdot,u);
		
		% produce output once every so often
		if (toc >= 0.1)
			tic;
			fprintf('Equilibrium eval # %4d: Norm of f: %12.5e\n', neval, norm(f));
% 			das3stick(x);
% 			drawnow;
		end
	end
		
	% ======================= do the derivatives test
	if strcmp(command, 'derivatives')
		disp('Checking dynamics derivatives...');
		for i=1:1							% increase this to repeat the test several times
			x = rand(nstates,1);
            x = xneutral;
			xdot = .1*randn(nstates,1);
			u = rand(nmus,1);
			[f,dfdx,dfdxdot,dfdu] = das3mex(x,xdot,u);
            dfdx_num	 	= zeros(nstates,nstates);
            dfdxdot_num		= zeros(nstates,nstates);
			dfdu_num 		= zeros(nstates,nmus);
			for j=1:nstates
				% derivatives with respect to x
				saved = x(j);
				h = 1e-7;
				x(j) = x(j)+h;
				[fnew] = das3mex(x,xdot,u);
				dfdx_num(:,j) = (fnew-f)/h;
				x(j) = saved;
				
				% derivatives with respect to xdot
				saved = xdot(j);
				if j>ndof && j<=2*ndof
					h = 1e10;		% for mass matrix, use large finite difference (it is not dependent on xdot!) to avoid roundoff error
				else
					h = 1e-5;		% for the other elements, this seems to be a good finite difference
				end
				xdot(j) = xdot(j)+h;
				[fnew] = das3mex(x,xdot,u);
				dfdxdot_num(:,j) = (fnew-f)/h;				
				xdot(j) = saved;
			end
			
			% derivatives with respect to u
			for j=1:nmus
				saved = u(j);
				h = 1e-7;
				u(j) = u(j) + h;
				[fnew] = das3mex(x,xdot,u);
				dfdu_num(:,j) = (fnew-f)/h;
				u(j) = saved;
			end

			% check for differences between numerical and analytical result
			fprintf('Checking df/dx...\n'); 		matcompare(dfdx, dfdx_num);
			fprintf('Checking df/dxdot...\n'); 		matcompare(dfdxdot, dfdxdot_num);
			fprintf('Checking df/du...\n');			matcompare(dfdu, dfdu_num);
			disp('If you want to inpect the Jacobians, look at dfdx, dfdxdot, and dfdu.');
			disp('The finite difference approximations are dfdx_num, dfdxdot_num, and dfdu_num.');
			disp('When done, type "dbcont"');
			keyboard
		end

    end		
    
    % ======================= do the isometric muscle tests
	if (strcmp(command, 'isometric') 	|| strcmp(command, 'all') )
		disp('Calculating isometric strength curves...');
		clf;
	   
        Lceopt = das3mex('LCEopt');
        
		% load equilbrium state x
		x = load('equilibrium.txt');
                
        % shoulder elevation moment-angle curves, 30-deg intervals in
        % elevation angle
        isometric_curves(x,Lceopt,dofnames,8,10:10:60,7,-90:30:90);

        % elbow flexion moment-angle curves, at 30-deg intervals in pronation angle
        isometric_curves(x,Lceopt,dofnames,10,10:10:140,11,10:30:160);

        % elbow pronation moment-angle curves, at 30-deg intervals in flexion angle
        isometric_curves(x,Lceopt,dofnames,11,10:10:160,10,10:30:140);

 		fprintf('Hit ENTER to continue.\n');
		pause
	end
	
	% ======================= do the isokinetic muscle tests
	if (strcmp(command, 'isokinetic') 	|| strcmp(command, 'all') )
		disp('Calculating isokinetic strength curves...');
		clf;
	   
        Lceopt = das3mex('LCEopt');
        
		% load equilbrium state x
		x = load('equilibrium.txt');
        
        % muscle moment arms in this position
        momentarms = full(das3mex('Momentarms', x));
        
		% moment-angular velocity curves
        for icurve=1:11
            isokinetic_curves(x,momentarms,Lceopt,dofnames,icurve,-1000:100:1000);
        end
        
        fprintf('Hit ENTER to continue.\n');
		pause
	end

	% ======================= do the explicit simulation test
	if strcmp(command, 'explicit')
        close all
		disp('Running simulation with ODE23...');
		neval = 0;
		nodefun = 0;
		printinterval = 3.0;
		x = load('equilibrium.txt');
		times = [0:0.01:4];
        starttime = tic;
		[tout, xout] = ode23(@odefun, times, x);
        cputime = toc - starttime;
        fprintf('Computation time: %8.3f seconds (%6.1f slower than real time)\n', cputime, cputime/max(times)); 
		fprintf('Number of integration steps:    %d\n', nodefun);
		fprintf('Number of function evaluations: %d\n', neval);
		
		% plots on screen
		figure(1);
		plotangles(tout,xout);
        figure(2);
        plot(tout(1:end-1), diff(tout));
        xlabel('time (s)');
        ylabel('ODE23 time step (s)');

		% export to mat file and opensim motion file
		save explicit_simulation tout xout
		make_osimm('explicit_simulation', xout(:,1:11), tout);

	end
	
	% ====================== simulate arm movement ==============
	% this is done with our own implicit algorithm, which will be used for real time simulation
	if strcmp(command, 'simulate')
		disp('Simulating arm movements...');
        
        % first run the simulation with the C function das3step (in das3mex.c)
        disp('das3steptest uses the das3step C function in das3mex.c');
        das3steptest;
				
        % now run it with the Matlab function das3step.m
		% set simulation parameters, default time step is 3 ms
        disp('Now running simulation with the das3step.m code...');
		t = 0;
		tend = 4.0;
		if (nargin < 2)
			tstep = .003;
		else
			tstep = arg1;
		end
		nsteps = round((tend-t)/tstep);
		
		% reserve space to store results
		tout = tstep*(0:nsteps)';
		xout = zeros(nsteps+1, nstates);
		tout(1) = t;
		
		% load equilbrium state x
		x = load('equilibrium.txt');
		
		% run simulation
		xout(1,:) = x';
		tic
		for i=1:nsteps
			u = stimfun(t);
			x = das3step(x, u, tstep);
			xout(i+1,:) = x';			% store result
			i = i+1;
			t = t + tstep;
		end
		simtime = toc;
		fprintf('CPU time per time step: %8.3f ms\n', 1000*simtime/nsteps);
		fprintf('Simulation speed is %8.3f times faster than real time\n',tend/simtime);
		
		% plots on screen
		figure(1);
		plotangles(tout,xout);

		% export to mat file and to opensim motion file
		save simulation tout xout
		make_osimm('simulation', xout(:,1:11), tout);
        
        % compare joint angles to the result of das3steptest
        x1 = load('das3steptest.txt');
        fprintf('joint angle difference between C and Matlab version of das3step: %12.5e deg\n', ...
            rms(rms(x1-xout))*180/pi);
		
		% compare joint angles to the explicit simulation test and report RMS difference
		tout_test = tout;
		xout_test = xout;
		if (exist('explicit_simulation.mat')==2)
			load('explicit_simulation');		% this loads tout and xout from the mat file
		else
			disp('Could not load explicit_simulation.mat');
			disp('Run das3test(''explicit'') first, if you want to quantify simulation errors.');
			return
		end
		
		% resample our implicit results to the same times as the explicit simulation
		xout_test = interp1(tout_test, xout_test, tout, 'linear', 'extrap');
		
		% compute the RMS error in joint angles, by comparing to explicit result
		angdiff = 180/pi*(xout(:,1:11) - xout_test(:,1:11));
		rmsdiff = sqrt(mean(angdiff.^2));
		disp('RMS errors in each joint angle (deg):');
		for i=1:ndof
			fprintf('%-6s: %8.3f\n', dofnames{i}, rmsdiff(i));
		end
		
		% output the overall RMS error
		out1 = sqrt(mean(mean(angdiff.^2)));
		fprintf('\nOverall RMS error: %8.3f degrees\n', out1);
		
	end

	%==================================================================================
	% The following function solves the state derivatives from state and control,
	% so we can simulate using a Matlab explicit ODE solver.
	% There are more efficient ways to do this, but this is just for model testing, not for efficient simulation.
	function [xdot] = odefun(t,x)
		persistent xdotinitialguess
		nodefun = nodefun+1;
		if ~exist('xdotinitialguess') || isempty(xdotinitialguess)
			xdotinitialguess = zeros(size(x));
		end
		xdot = xdotinitialguess;
		u = stimfun(t);
		if 0 && exist('fsolve.m') > 1
			% use matlab's fsolve if optimization toolbox is installed
			options = optimset(optimset('fsolve'),'Jacobian','on','Display','off');
			[xdot, f] = fsolve(@zero,xdot,options);
            resnorm = norm(f);
		else
			% use a simple Newton-Raphson method
			% this actually turns out to be much faster than fsolve and equally accurate
			resnorm = 1e10;
			while resnorm > 1e-8
				[f,J] = zero(xdot);
				xdot = xdot - (J\f);		% do one full Newton step
				resnorm = norm(f);			% calculate norm of residuals
			end
		end
		% give some output on screen, every once in a while
		if (toc > printinterval)
			tic;
			fprintf('Step %d: Neval = %d -- t = %20.14g -- Norm(f) = %10.3e\n', nodefun, neval, t, resnorm);
			if (printinterval == 0.0) pause; end
		end

		% update the xdot initial guess for a better start of the next step
		xdotinitialguess = xdot;
		nodefun = nodefun + 1;
		
		function [f,J] = zero(xdot)
			[f,~,J,~] = das3mex(x, xdot, u);
			neval = neval+1;
		end
    end	
end
%=====================================================
function matcompare(a,b);
	% compares two matrices and prints element that has greatest difference
	[maxerr,irow] = max(abs(a-b));
	[maxerr,icol] = max(maxerr);
	irow = irow(icol);
	fprintf('Max. difference: %e at %d %d (%e analytical vs. %e numerical)\n', ...
		maxerr, irow, icol, full(a(irow,icol)),full(b(irow,icol)));
end
%==========================================================
function plotangles(t,x)
	subplot(2,2,1)
	plot(t,180/pi*x(:,1:3));
	legend('SC\_roty','SC\_rotz','SC\_rotx');
	xlabel('time (s)');
	ylabel('angle (deg)');
	
	subplot(2,2,2)
	plot(t,180/pi*x(:,4:6));
	legend('AC\_roty','AC\_rotz','AC\_rotx');
	xlabel('time (s)');
	ylabel('angle (deg)');
	
	subplot(2,2,3)
	plot(t,180/pi*x(:,7:9));
	legend('GH\_roty','GH\_rotz','GH\_rotyy');
	ylabel('angle (deg)');
	xlabel('time (s)');
	
	subplot(2,2,4)
	plot(t,180/pi*x(:,10:11));
	legend('EL\_rotx','PS\_roty');
	ylabel('angle (deg)');
	xlabel('time (s)');

end
%==========================================================
function [u] = stimfun(t)
	nmus = 138;
	trampup = 0.2;		% time for the ramp up part of the stim pattern
	tpassive = 2.0;		% time for muscle to become passive

	if t<trampup
			u = t/trampup + zeros(nmus,1);
	elseif t<tpassive
			u = ones(nmus,1);
	else 
		u = zeros(nmus,1);
	end
end
%==========================================================
function isometric_curves(x, Lceopt, dofnames, joint1, range1, joint2, range2)
    % produces isometric moment-angle curves for a joint
    % one for each value in a range of angles in second joint

    % angles in degrees
    angles = x(1:11)'*180/pi;
    angvel = zeros(size(angles));

    pascurves = [];
    poscurves = [];
    negcurves = [];
    legends = {};
    for angle2 = range2
        angles(joint2) = angle2;
        fprintf('isometric simulations for %s at %8.3f deg\n', dofnames{joint2},angle2);
        x(joint2) = angle2*pi/180;
        pasmoments = [];
        posmoments = [];
        negmoments = [];
        for angle1 = range1
            angles(joint1) = angle1;
            x(joint1) = angle1*pi/180;

            % find sign of moment arm in this position
            momentarms = full(das3mex('Momentarms', x));

            mom = maxmoment(joint1, angles, angvel, momentarms, Lceopt, 0);
            pasmoments = [pasmoments mom];
            mom = maxmoment(joint1, angles, angvel, momentarms, Lceopt, 1);
            posmoments = [posmoments mom];
            mom = maxmoment(joint1, angles, angvel, momentarms, Lceopt, -1);
            negmoments = [negmoments mom];
        end
        pascurves = [pascurves  pasmoments'];
        poscurves = [poscurves  posmoments'];
        negcurves = [negcurves  negmoments'];
        legends = [legends ; [char(dofnames(joint2)) ' ' num2str(angle2)] ];
    end

    % find max and min moment to set axis limits
    allcurves = [pascurves poscurves negcurves poscurves-pascurves negcurves-pascurves];
    allcurves = reshape(allcurves,[],1);
    minmom = min(allcurves);
    maxmom = max(allcurves);
    momrange = maxmom-minmom;

    figure;
    % plot total moments on left side of figure
    subplot(2,3,1);
    plot(range1, poscurves,'x-');
    ylim([minmom-0.1*momrange maxmom+0.1*momrange]);
    title([char(dofnames(joint1)) ': positive moment']);
    ylabel('moment (Nm)');
    subplot(2,3,4);
    plot(range1, negcurves,'x-');
    ylim([minmom-0.1*momrange maxmom+0.1*momrange]);
    xlabel('angle (deg)');
    ylabel('moment (Nm)');
    title([char(dofnames(joint1)) ': negative moment']);

    % plot passive moments in middle column of figure
    subplot(2,3,2);
    plot(range1, pascurves,'x-');
    ylim([minmom-0.1*momrange maxmom+0.1*momrange]);
    title([char(dofnames(joint1)) ': passive moment']);
    xlabel('angle (deg)');

    % subtract passive moments and plot in rightmost column of figure
    subplot(2,3,3);
    plot(range1, poscurves-pascurves,'x-');
    ylim([minmom-0.1*momrange maxmom+0.1*momrange]);
    title([char(dofnames(joint1)) ': positive - passive']);
    subplot(2,3,6);
    plot(range1, negcurves-pascurves,'x-');
    ylim([minmom-0.1*momrange maxmom+0.1*momrange]);
    title([char(dofnames(joint1)) ': negative - passive']);
    xlabel('angle (deg)');

    legend(legends);        
end

%=============================================================================================================
function isokinetic_curves(x, momentarms, Lceopt, dofnames, joint, range)
    % produces isokinetic moment-angular velocity curves for a joint
    %
    % x             the state vector
    % momentarms    the 138 muscle moment arms
    % Lceopt        the 138 muscle optimal fibre lengths
    % dofnames      the names of the 11 degrees of freedom
    % joint         for which joint we will calculate curve
    % range         the range of velocities
    
    % angles in degrees
    angles = x(1:11)'*180/pi;
    angvel = zeros(size(angles));
    pasmoments = [];
    posmoments = [];
    negmoments = [];
    for vel = range
        angvel(joint) = vel;
        fprintf('isokinetic simulations for %s at %8.3f deg/s\n', dofnames{joint},vel);
        pasmoments = [pasmoments maxmoment(joint, angles, angvel, momentarms, Lceopt, 0)];
        posmoments = [posmoments maxmoment(joint, angles, angvel, momentarms, Lceopt, 1)];
        negmoments = [negmoments maxmoment(joint, angles, angvel, momentarms, Lceopt, -1)];
    end

    % find max and min moment to set axis limits
    allcurves = [pasmoments posmoments negmoments posmoments-pasmoments negmoments-pasmoments];
    allcurves = reshape(allcurves,[],1);
    minmom = min(allcurves);
    maxmom = max(allcurves);
    momrange = maxmom-minmom;

    figure;
    % plot total moments on left side of figure
    subplot(2,3,1);
    plot(range, posmoments,'x-');
    ylim([minmom-0.1*momrange maxmom+0.1*momrange]);
    title([char(dofnames(joint)) ': positive moment']);
    ylabel('moment (Nm)');
    subplot(2,3,4);
    plot(range, negmoments,'x-');
    ylim([minmom-0.1*momrange maxmom+0.1*momrange]);
    xlabel('angular velocity (deg/s)');
    ylabel('moment (Nm)');
    title([char(dofnames(joint)) ': negative moment']);

    % plot passive moments in middle column of figure
    subplot(4,3,[2 3]);
    plot(range, pasmoments,'x-');
    ylim([minmom-0.1*momrange maxmom+0.1*momrange]);
    title([char(dofnames(joint)) ': passive moment']);
    xlabel('angular velocity (deg/s)');

    % subtract passive moments and plot in rightmost column of figure
    subplot(2,3,3);
    plot(range, posmoments-pasmoments,'x-');
    ylim([minmom-0.1*momrange maxmom+0.1*momrange]);
    title([char(dofnames(joint)) ': positive - passive']);
    subplot(2,3,6);
    plot(range, negmoments-pasmoments,'x-');
    ylim([minmom-0.1*momrange maxmom+0.1*momrange]);
    title([char(dofnames(joint)) ': negative - passive']);
    xlabel('angular velocity (deg/s)');

end

%=============================================================================================================
function mom = maxmoment(joint, angles, angvel, momentarms, Lceopt, sign)
    % simulate maximum moment at one joint, as function of all joint angles and angular velocities
    % joint         for which joint we will calculate moment
    % angles        the eleven joint angles (deg)
    % angvel        the eleven angular velocities (deg/s)
    % momentarms    the 138 muscle moment arms
    % Lceopt        the 138 muscle optimal fibre lengths
    % sign          0: passive, 1: max positive moment, -1: max negative moment

    angles = angles*pi/180;		% convert to radians
    angvel = angvel*pi/180;
    ndof = length(angles);
    nmus = length(Lceopt);

    Act =  sign*momentarms(:,joint) > 0;	% vector that has a 1 for all muscles we want to activate

    % determine lengthening velocities of the muscle-tendon complexes, normalize to Lceopt
    Vmuscle = -(momentarms * angvel') ./ Lceopt;

    % determine the Lce's at which there is contraction equilibrium (dF/dt = 0, or Lcedot = Vmuscle)
    % use Newton's method
    % we want these state derivatives:
    xdot = [zeros(2*ndof,1);Vmuscle;zeros(nmus,1)];	
    u = zeros(nmus,1);				% no stim, we don't care about activation dynamics here
    x = [angles angvel zeros(1,nmus) Act']';

    % unfortunately, fzero (1-dimensional root finding) is the only solver
    % that is robust enough here.  This makes the code slow.
    % Multidimensional solvers, doing all muscles at the same time, fsolve,
    % lmopt, and simple Newton method all failed.
    options = optimset('TolX', 1e-6);
	for imus=1:nmus                 % go through all muscles
        Lce = 1.0;                  % initial guess for this muscle's Lce
        [Lce, Fval, flag] = fzero(@contraction_equilibrium, Lce, options);  % find Lce at contraction equilibrium
        if (flag < 0)
            fprintf('maxmoment: muscle contraction equilibrium not found within max number of iterations.\n');
            keyboard
        end
	end 

% now determine the joint moments at this state of the system
    moments = das3mex('Jointmoments', x);
    mom = moments(joint);
    
    function [F] = contraction_equilibrium(Lce)
        x(2*ndof+imus) = Lce;
        f = das3mex(x, xdot, u);
        F = f(2*ndof+imus);
    end

end


