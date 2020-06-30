function das3gui
% GUI to run real-time simulations with DAS3.
% The model is shown as a simple stick figure.
% User can press buttons to activate muscles.  Each button click will
% increase muscle excitation by 10%.  Each 50 ms, muscle excitations will
% decay by 1%.
%
% 1. Separated the deltoid into three functional parts
% (anterior-middle-posterior) instead of two anatomical ones
% (clavicular-scapular)
% 2. Changed degrees to radians for Opensim motion file
% 3. Corrected timer stop error
% 4. Changed default "data" value from zero to equilibrium position
% Dimitra Blana, February 2012
%
% Added plot of glenohumeral force in glenoid
% Dimitra Blana, March 2012

	global muscles stim x FGH buttons mainfigure data idata steptime das3timer ha hb


	% some constants
	ndof = 11;
	nmus = 138;						% number of muscle elements
	nstates = 2*ndof + 2*nmus;
	steptime = 0.05;				% updating screen 20 times per second is fast enough
	nframes = 200;					% number of frames stored in circular buffer
	    
	% define the 30 muscle groups
	muscles = { ...
		% name								first element	last element
		'trapezius scap. part'                       1       11  ; ...
		'trapezius clav. part'                       12      13  ; ...
		'levator scapulae'                           14      15  ; ...
		'pectoralis minor'                           16      19  ; ...
		'rhomboideus'                                20      24  ; ...
		'serratus anterior'                          25      36  ; ...
		'deltoideus, posterior'                      37      41  ; ...
		'deltoideus, middle'                         42      47  ; ...
		'deltoideus, anterior'                       48      51  ; ...
		'coracobrachialis'                           52      54  ; ...
		'infraspinatus'                              55      60  ; ...
		'teres minor'                                61      63  ; ...
		'teres major'                                64      67  ; ...
		'supraspinatus'                              68      71  ; ...
		'subscapularis'                              72      82  ; ...
		'biceps, caput longum'                       83      83  ; ...
		'biceps, caput breve'                        84      85  ; ...
		'triceps, caput longum'                      86      89  ; ...
		'latissimus dorsi'                           90      95  ; ...
		'pect. major, thor. part'                    96     101  ; ...
		'pect. major, clav. part'                   102     103  ; ...
		'triceps, medial part'                      104     108  ; ...
		'brachialis'                                109     115  ; ...
		'brachioradialis'                           116     118  ; ...
		'pronator teres, hum-rad'                   119     119  ; ...
		'pronator teres, uln-rad'                   120     120  ; ...
		'supinator, uln-rad'                        121     125  ; ...
		'pronator quadratus'                        126     128  ; ...
		'triceps, lateral part'                     129     133  ; ...
		'anconeus'                                  134     138 };
	nmuscles = size(muscles,1);
	
	% initialize the stim levels
	stim = zeros(nmuscles,1);

	% Create the GUI
	
	% screen size
	screen = get(0,'MonitorPositions');
	screenwidth = screen(3);
	screenheight = screen(4);
	
	% window size (pixels)
	width = 600;
	height = 650;
	
	% width and spacing of the pushbuttons (pixels)
	spacing = 5;
	xsize = 120;
	
	% make the figure window
	close all
	mainfigure = figure('Position',[7,screenheight-height-70,width,height], ...
		'Name', 'DAS3 Real Time Simulation', ...
		'NumberTitle', 'off', ...
		'ToolBar', 'figure');
   
	% construct the plot area for the stick figure, the left part of the window
	buttonarea = xsize/width;
	ha = axes('Position',[0,0,1-buttonarea,0.65]);     
	hb = axes('Position',[0,0.7,1-buttonarea,0.2]); 
	
	% construct the push buttons for each muscle, on right side of window
	ysize = height/nmuscles - spacing;
	for i=1:nmuscles
		xpos = width - xsize;
		ypos = height - (ysize+spacing)*i;
		buttons(i) = uicontrol('Style','pushbutton','String',muscles{i,1},...
          'Position',[xpos,ypos,xsize,ysize],...
          'Callback',{@button_Callback});
	end
	
	% construct the quit button
	uicontrol('Style','pushbutton','String','QUIT',...
          'Position',[1,1,60,30],...
          'Callback',{@quit_button});
	
	% construct the reset button
	uicontrol('Style','pushbutton','String','RESET',...
          'Position',[65,1,60,30],...
          'Callback',{@reset_button});
	
	% initialize the das3 model
	das3mex();

	% load equilbrium state x
	load('equilibrium.mat');
    [~, FGH] = das3step(x, zeros(nmus,1), 0.005);

	% create file to store simulation result
	%data = zeros(nframes,ndof);
	data = repmat(x(1:11)',nframes,1);   % equilibrium position
	idata = 1;

	% initialize timer     
	das3timer = timer('TimerFcn',@das3update,'StopFcn',@das3stop,...
        'Period', steptime, 'ExecutionMode', 'fixedDelay', 'BusyMode','queue');
	start(das3timer);

end

% Callback functions

function quit_button(source, eventdata)
	global das3timer
    
    stop(das3timer);    
end

function das3stop(source, eventdata)
	global mainfigure steptime data idata

	delete(mainfigure);
	    
	% write the .sto file for Opensim
	filename = 'das3gui_result';
	numframes = size(data,1);
	fid = fopen([filename '.sto'],'wt');
	fprintf(fid,'%s\n',filename);
	fprintf(fid,'%s%i\n','nRows=',numframes);
	fprintf(fid,'%s\n','nColumns=12');
	fprintf(fid,'%s\n','endheader');
	fprintf(fid,'%s\n','time  SC_y  SC_z SC_x  AC_y  AC_z  AC_x  GH_y  GH_z  GH_yy  EL_x  PS_y');
	for i = idata:(idata+numframes-1)
		j = mod(i,numframes) + 1;
		fprintf(fid,'%f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n',(i-idata)*steptime,data(j,:));
	end
	fclose(fid);
end


function reset_button(source, eventdata)
	global x stim hb FGH
	load('equilibrium.mat');
	stim = zeros(size(stim));
    axes(hb);
    hold off;
    [~, FGH] = das3step(x, zeros(138,1), 0.005);
end

function button_Callback(source,eventdata) 
% respond to push button click

	global stim buttons

	% determine which button was pressed
	num = find(source == buttons);
	
	% increase the stim of all elements of this muscle
	stim(num) = stim(num) + 0.1;
	stim = min(stim,1.0);
	
end

function das3update(obj, event, string_arg)
% display and simulate the das3, whenever timer requests it

	global x FGH stim buttons muscles data idata ha hb

	% draw the model
    axes(ha);
	das3stick(x);
    axes(hb);
    plot_GH(FGH);
	drawnow;
	
	% color the buttons according to muscle's stim level
	bgcolor = [0.831373 0.815686 0.784314];  % the standard grey button collor
	red = [1 0 0];
	for i = 1:size(stim,1)
		button = buttons(i);
		color = bgcolor + stim(i)*(red - bgcolor);
		set(button,'BackgroundColor',color);
	end
	
	% stimulate the muscle elements
	u = zeros(138,1);
	for i=1:size(stim,1)
		elements = muscles{i,2}:muscles{i,3};
		u(elements) = stim(i);
	end
	
	% simulate the model for 50 ms, in steps of 5 ms
	step = 0.005;
	for i=1:10
		[x, FGH] = das3step(x, u, step);
	end
	
	% store in circular buffer, so we can save result on file later
%	data(idata,:) = x(1:11)*180/pi;  % (in degrees) for the Opensim <--
	data(idata,:) = x(1:11);  		 % time and angles (in radians!) for the Opensim file
	idata = idata + 1;
	if ( idata > size(data,1) )
		idata = 1;
	end

	% decay all muscle stimulations by 1%
	stim = stim*0.99;
	
end

function das3stick(x)
% draw a stick figure of the DAS3 model in state(s) x.

	% radii of thorax ellipsoid (must be same as in das3.al)
	Ax  =  0.1470;
	Ay  =  0.2079;
	Az  =  0.0944;	
	
	% define the points used to draw each bone, and their color, and whether to draw filled polygon
	% indices refer to how they come back from the MEX function: IJ, SC, AC, TS, AI, GH, UL, RD, HD
	bonepoints = {
		[2 3],		    [0    0    1   ]	0;			% clavicle: SC and AC
		[3 4 5 3 6],	[0    0.50 0   ]	1;			% scapula: AC, TS, AI, back to AC, and GH
		[6 7],			[1    0    0   ]	0;			% humerus: GH and UL
		[7 9],			[0    0.75 0.75]	0;			% ulna: UL and HD
		[8 9],			[0.75 0    0.75]	0;			% radius: RD and HD
	};
	
	cla

	% draw the thorax ellipsoid
	% (in das3.al, the ellipsoid is centered at the origin)
	[xx,yy,zz] = ellipsoid(0,0,0,Ax,Ay,Az);
	[row,col,v] = find(zz>0);
	newx = xx(row,col);
	newy = yy(row,col);
	newz = zz(row,col);
	mesh(newz,newx,newy);
	
	% get the 3D coordinates of the stick figure points
	d = das3mex('Stick',x);
		
	% draw the stick figure
	hold on
	for j = 1:size(bonepoints,1);
		list = bonepoints{j,1};
		x = d(1,list);
		y = d(2,list);
		z = d(3,list);
		filled = bonepoints{j,3};
		c = bonepoints{j,2};		% color
		if (filled)
			fill3(z,x,y,c);
		end
		% plot necessary also after fill, to remove black edge of polygon
		plot3(z,x,y,'Color',c,'Marker','o','LineWidth',2)
	end
	
	axis('equal');
	axis([-0.5 0.5 -0.5 0.5 -0.5 0.5]);
	axis('off');

end

function plot_GH(FGH)

    Rgt = glenoid_scap;
    Fgh0 = Rgt*FGH;  
    if norm(Fgh0), Fgh0 = Fgh0/norm(Fgh0); end
    % decompose into polar angles
    thetar = asin(-Fgh0(2));
    if ~(sqrt(Fgh0(1)^2+Fgh0(3)^2))
        phir = 0.0;
    else
        phir=asin(Fgh0(3)/sqrt(Fgh0(1)^2+Fgh0(3)^2));
    end

    aphi=tand(38.55);
    atheta=tand(44.37);
    
    Plellips(0,0,aphi,atheta);
    axis equal
    axis off
    box off
    ylabel('Posterior','fontsize',12);
    set(gca,'xtick',[],'ytick',[]);
    hold on;
    scatter(-tan(phir),-tan(thetar),'o'); 

end 


function Plellips(Mx,My,ax,ay)
    % routine to plot an ellipse
    step=(2*ax)/100;
    x=Mx-ax:step:Mx+ax;
      y1= sqrt((ones(1,101)-((x-Mx*ones(1,101))./(ax*ones(1,101))).^2).*((ay*ones(1,101)).^2)) + My*ones(1,101);
      y2=-sqrt((ones(1,101)-((x-Mx*ones(1,101))./(ax*ones(1,101))).^2).*((ay*ones(1,101)).^2)) + My*ones(1,101);

    plot(x,y1,'k-',x,y2,'k-')
end