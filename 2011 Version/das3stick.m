function das3stick(state, makemovie)
% das3stick.m: draw a stick figure of the DAS3 model in state(s) x.
% 
% Inputs:
%	state 		(298 x N matrix) 	Model state (N=1), or a series of states
%	makemovie	(logical)			Flag to indicate if movie file should be made

	% default is no movie
	if (nargin < 2)
		makemovie = 0;
	end

	% radii of thorax ellipsoid (must be same as in das3.al)
	Ax  =  0.1470;
	Ay  =  0.2079;
	Az  =  0.0944;	
	
	% set plotting volume (TODO: Dimitra, I am not sure which region should be plotted, we don't need to see the whole ellipsoid
	xrange = [-0.3 0.5];
	yrange = [-0.2 0.7];
	zrange = [-0.4 0.4];
	
	% define the points used to draw each bone, and their color, and whether to draw filled polygon
	% indices refer to how they come back from the MEX function: IJ, SC, AC, TS, AI, GH, UL, RD, HD
	bonepoints = {
		[2 3],		    [0    0    1   ]	0;			% clavicle: SC and AC
		[3 4 5 3 6],	[0    0.50 0   ]	1;			% scapula: AC, TS, AI, back to AC, and GH
		[6 7],			[1    0    0   ]	0;			% humerus: GH and UL
		[7 9],			[0    0.75 0.75]	0;			% ulna: UL and HD
		[8 9],			[0.75 0    0.75]	0;			% radius: RD and HD
	};
	
	if makemovie
		avi = avifile('movie.avi', 'fps', 15);
	end

	cla
	for i=1:size(state,2)
		
		% draw the thorax ellipsoid
		% (in das3.al, the ellipsoid is centered at the origin)
		[xx,yy,zz] = ellipsoid(0,0,0,Ax,Ay,Az);
        [row,col] = find(zz>0);
        newx = xx(row,col);
        newy = yy(row,col);
        newz = zz(row,col);
        mesh(newz,newx,newy);
		
		% get the 3D coordinates of the stick figure points
		d = das3mex('Stick',state(:,i));
		
		% draw the stick figure
		hold on
		for j = 1:size(bonepoints,1)
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
		% set the axes and viewpoint
		axis('equal');
		axis([xrange yrange zrange]);
		%axis('off')
        xlabel('x');
        ylabel('y');
        zlabel('z');
        view(157, 6);

		if (makemovie)
			if (i==1)
				F = getframe(gca);
				frame = [1 1 size(F.cdata,2) size(F.cdata,1)];
			else
				F = getframe(gca,frame);
			end
			avi = addframe(avi,F);
			cla
		end
	end
	if (makemovie)
		avi = close(avi);
	end
	hold off;
end
