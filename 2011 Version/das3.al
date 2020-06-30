% das3.al
%
% This is an input file for Autolev to generate dynamic equations of motion for the DAS3 model.
%
% Author: Dimitra Blana
% Based on walk2d.al by Ton van den Bogert

AUTOZ ON
OVERWRITE ALL
PAUSE 0

%--------------------------------------------------------------------
%	Thorax reference frame, units
%--------------------------------------------------------------------
Newtonian   Thorax				% reference frame that does not move
UnitSystem  kg,meter,sec

% Origin		Center of thorax ellipsoid
% X				Lateral (= to the right, for right shoulder which is modeled here)
% Y				Superior
% Z				Posterior
% Bone reference frames are described in das3implicit.doc

%----------------------------------------
% Define bodies and their mass properties
%----------------------------------------
Bodies Clavicle, Scapula, Humerus, Ulna, Radius, Hand

% masses (values are set in das3mex.c)
Constants par__Mass_Clavicle, par__Mass_Scapula, par__Mass_Humerus, par__Mass_Ulna, par__Mass_Radius, par__Mass_Hand
Mass	Clavicle 	= par__Mass_Clavicle
Mass	Scapula		= par__Mass_Scapula
Mass	Humerus		= par__Mass_Humerus
Mass	Ulna		= par__Mass_Ulna
Mass	Radius		= par__Mass_Radius
Mass	Hand		= par__Mass_Hand

% coordinates of segment centers of mass (in their own reference frame)
Constants par__ClavCMx, par__ClavCMy, par__ClavCMz 
Constants par__ScapCMx, par__ScapCMy, par__ScapCMz 
Constants par__HumCMx,  par__HumCMy , par__HumCMz  
Constants par__UlnaCMx, par__UlnaCMy, par__UlnaCMz 
Constants par__RadCMx,  par__RadCMy , par__RadCMz  
Constants par__HandCMx, par__HandCMy, par__HandCMz 

% inertia matrices
Constants par__ClavIxx, par__ClavIyy, par__ClavIzz, par__ClavIxy, par__ClavIyz, par__ClavIxz
Constants par__ScapIxx, par__ScapIyy, par__ScapIzz, par__ScapIxy, par__ScapIyz, par__ScapIxz
Constants par__HumIxx,  par__HumIyy,  par__HumIzz,  par__HumIxy,  par__HumIyz,  par__HumIxz
Constants par__UlnaIxx, par__UlnaIyy, par__UlnaIzz, par__UlnaIxy, par__UlnaIyz, par__UlnaIxz
Constants par__RadIxx,  par__RadIyy,  par__RadIzz,  par__RadIxy,  par__RadIyz,  par__RadIxz
Constants par__HandIxx, par__HandIyy, par__HandIzz, par__HandIxy, par__HandIyz, par__HandIxz
Inertia  Clavicle, 	par__ClavIxx, par__ClavIyy, par__ClavIzz, par__ClavIxy, par__ClavIyz, par__ClavIxz
Inertia  Scapula, 	par__ScapIxx, par__ScapIyy, par__ScapIzz, par__ScapIxy, par__ScapIyz, par__ScapIxz
Inertia  Humerus, 	par__HumIxx,  par__HumIyy,  par__HumIzz,  par__HumIxy,  par__HumIyz,  par__HumIxz
Inertia  Ulna, 		par__UlnaIxx, par__UlnaIyy, par__UlnaIzz, par__UlnaIxy, par__UlnaIyz, par__UlnaIxz
Inertia  Radius, 	par__RadIxx,  par__RadIyy,  par__RadIzz,  par__RadIxy,  par__RadIyz,  par__RadIxz
Inertia  Hand, 		par__HandIxx, par__HandIyy, par__HandIzz, par__HandIxy, par__HandIyz, par__HandIxz

%------------------------------------------------
% Define points
%------------------------------------------------
% Bony landmarks from Wu et al., 2004
% IJ:    deepest point of Incisura Jugularis(thorax)
% SC:    most ventral point on the sternoclavicular joint(clavicle)
% AC:    most dorsal point on the acromioclavivular joint(scapula)
% TS:    Trigonum Spinae (scapula)
% AI:    Angulur Inferior (scapula)

% Joint centers, from van der Helm, 1994:
% GH:    center of glenoid (scapula)
% UL:    center of humeral-ulnar joint (ulna)
% RD:    center of ulna-radial joint (radius)
% HD:    center of radio-carpal joint (hand)

Points IJ, SC, AC, GH, UL, RD, HD

% Note: the points TS and AI will be created later by calling the contact.a script

%--------------------------------------------------------------------
% morphological constants (values are set in das3mex.c)
%--------------------------------------------------------------------
Constants  par__IJx, par__IJy, par__IJz				% coordinates of IJ in thorax frame
Constants  par__SCx, par__SCy, par__SCz				% coordinates of SC in thorax frame
Constants  par__ACx, par__ACy, par__ACz				% coordinates of AC in clavicle frame
Constants  par__GHx, par__GHy, par__GHz				% coordinates of GH in scapula frame
Constants  par__ULx, par__ULy, par__ULz				% coordinates of UL in humerus frame
Constants  par__RDx, par__RDy, par__RDz				% coordinates of RD in ulna frame
Constants  par__HDx, par__HDy, par__HDz				% coordinates of HD in radius frame
Constants  par__Ax, par__Ay, par__Az				% radii of the thorax ellipsoid
Constants  par__TSx, par__TSy, par__TSz				% coordinates of TS in scapula frame
Constants  par__AIx, par__AIy, par__AIz				% coordinates of AI in scapula frame
Constants  par__FlexAxisX, par__FlexAxisY, par__FlexAxisZ	% Elbow flexion axis orientation, in humerus frame
Constants  par__PronAxisX, par__PronAxisY, par__PronAxisZ	% Forearm pronation axis orientation, in ulna frame

%---------------------------------------------------------------------
% other constants (values are set in das3mex.c)
%---------------------------------------------------------------------
Constants par__epscontact, par__kcontact					% parameters for scapula-thorax contact model
Constants par__epsconoid, par__conK, par__conL				% parameters for the conoid ligament model
Constants par__conOx, par__conOy, par__conOz                % coordinates of conoid origin in clavicle frame
Constants par__conIx, par__conIy, par__conIz                % coordinates of conoid insertion in scapula frame

%---------------------------------------------------
% Generalized coordinates
%---------------------------------------------------

% degrees of freedom from Wu et al. 2004

% sternoclavicular joint:
% 	q1	SC_roty: clavicle protraction/retraction
% 	q2	SC_rotz: clavicle elevation/depression
% 	q3	SC_rotx: clavicle axial rotation
% acromioclavicular joint:
% 	q4	AC_roty: scapular protraction/retraction
% 	q5	AC_rotz: scapular lateral/medial rotation
% 	q6	AC_rotx: scapular tilt
% glenohumeral joint:
% 	q7	GH_roty: humerus plane of elevation
% 	q8	GH_rotz: humerus elevation angle
% 	q9	GH_rotyy: humerus axial rotation
% humero-ulnar joint:
% 	q10	EL_rotx: elbow flexion/extension
% ulno-radial joint:
% 	q11	PS_roty: forearm pronation/supination

MotionVariables' q{11}''

%--------------------------------------------------------------------------------------------
% Clavicle Segment, attached to thorax via three hinges and two massless intermediate frames
%--------------------------------------------------------------------------------------------
Frames ClavicleA, ClavicleB
Simprot(Thorax,    ClavicleA, 2, q1)
Simprot(ClavicleA, ClavicleB, 3, q2)
Simprot(ClavicleB, Clavicle,  1, q3)
P_ThoraxO_SC> = Vector(Thorax, par__SCx, par__SCy, par__SCz)					% position of SC joint
V_SC_Thorax> = 0>																% velocity of SC joint
P_SC_ClavicleO> = Vector(Clavicle, par__ClavCMx, par__ClavCMy, par__ClavCMz)	% position of Clavicle CM						
V2pts(Thorax, Clavicle, SC, ClavicleO) 											% velocity of Clavicle CM
	
%-----------------------------------
% Scapula Segment
%-----------------------------------
Frames ScapulaA, ScapulaB
Simprot(Clavicle, ScapulaA, 2, q4)
Simprot(ScapulaA, ScapulaB, 3, q5)
Simprot(ScapulaB, Scapula,  1, q6);							% orientation of Scapula
P_SC_AC> = Vector(Clavicle, par__ACx, par__ACy, par__ACz)						% position of AC joint
V2pts(Thorax, Clavicle, SC, AC)													% velocity of AC joint													
P_AC_ScapulaO> = Vector(Scapula, par__ScapCMx, par__ScapCMy, par__ScapCMz)		% position of scapula CM
V2pts(Thorax, Scapula, AC, ScapulaO) 											% velocity of scapula CM
		
%-----------------------------------
% Humerus Segment
%-----------------------------------
Frames HumerusA, HumerusB
Simprot(Scapula,  HumerusA, 2, q7)
Simprot(HumerusA, HumerusB, 3, q8)
Simprot(HumerusB, Humerus,  2, q9)												% orientation of Humerus
P_AC_GH> = Vector(Scapula, par__GHx, par__GHy, par__GHz)						% position of GH joint
V2pts(Thorax, Scapula, AC, GH)													% velocity of GH joint																							
P_GH_HumerusO> = Vector(Humerus, par__HumCMx, par__HumCMy, par__HumCMz)   		% position of humerus CM
V2pts(Thorax, Humerus, GH, HumerusO) 											% velocity of humerus CM
		
%-----------------------------------
% Ulna Segment
%-----------------------------------
FlexionAxis> = Vector(Humerus, par__FlexAxisX, par__FlexAxisY, par__FlexAxisZ)
% FlexionAxis> = Vector(Humerus, 1, 0, 0)											% is not exact, but code will run faster if we do this
Simprot(Humerus, Ulna, FlexionAxis>, q10);										% orientation of Ulna
P_GH_UL> = Vector(Humerus, par__ULx, par__ULy, par__ULz)						% position of UL
V2pts(Thorax, Humerus, GH, UL)													% velocity of UL											
P_UL_UlnaO> = Vector(Ulna, par__UlnaCMx, par__UlnaCMy, par__UlnaCMz)			% position of Ulna CM
V2pts(THorax, Ulna, UL, UlnaO) 													% velocity of Ulna CM				
		
%-----------------------------------
% Radius Segment
%-----------------------------------
PronationAxis> = Vector(Ulna, par__PronAxisX, par__PronAxisY, par__PronAxisZ)
% PronationAxis> = vector(Ulna, 0, 1, 0)										% is not exact, but code will run faster if we do this
Simprot(Ulna, Radius, PronationAxis>, q11);										% orientation of Radius
P_UL_RD> = Vector(Ulna, par__RDx, par__RDy, par__RDz)							% position of RD
V2pts(Thorax, Ulna, UL, RD)														% velocity of RD
P_RD_RadiusO> = Vector(Radius, par__RadCMx, par__RadCMy, par__RadCMz)			% position of Radius CM
V2pts(Thorax, Radius, RD, RadiusO) 												% orientation of Radius CM
		
%-----------------------------------
% Hand Segment
%-----------------------------------
Radius_Hand = [1, 0, 0; 0, 1, 0; 0, 0, 1]						% no rotation between Radius and Hand frames
P_RD_HD> = Vector(Radius, par__HDx, par__HDy, par__HDz)						% position of HD
V2pts(Thorax, Radius, RD, HD)									% velocity of HD
P_HD_HandO> = Vector(Hand, par__HandCMx, par__HandCMy, par__HandCMz)			% position of Hand CM
V2pts(Thorax, Hand, HD, HandO) 									% velocity of Hand CM

%--------------------------------------------------------------------
% Apply gravity
%--------------------------------------------------------------------
Constants   g=9.81				% Local gravitational acceleration
Gravity( -g*Thorax2> )			% gravity acts along -y axis of Thorax

%--------------------------------------------------------------------
% Apply force on the forearm (ulna) to simulate a mobile arm support
%--------------------------------------------------------------------
Points ForF
Variables  distF, ampF
P_UL_ForF> = Vector(Ulna,0,-distF,0)								
V2pts(Thorax, Ulna, UL, ForF)
Force_ForF> = ampF*Thorax2> 

%--------------------------------------------------------------------
% Calculate and apply conoid force
%--------------------------------------------------------------------
Points conO, conI							% Conoid origin (on the clavicle) and insertion (on the scapula)
P_SC_conO> = Vector(Clavicle,par__conOx,par__conOy,par__conOz)					
P_AC_conI> = Vector(Scapula,par__conIx,par__conIy,par__conIz)				
V2pts(Thorax, Clavicle, SC, conO)
V2pts(Thorax, Scapula, AC, conI)		
LOI = Mag( P_conO_conI> ) 												% Distance between conoid origin and insertion
Stretch = LOI - par__conL												% Stretch of conoid "spring"
StretchPositive := 0.5*(Stretch + sqrt(Stretch^2 + par__epsconoid^2))	% Eliminate negative stretches
Uvec> = P_conO_conI> / LOI 												% Unit vector from conoid origin to insertion
Force(conO/conI, -par__conK*StretchPositive*Uvec> ) 				    % Spring force
 
%----------------------------------------------------------------------------
% Moments (exactly one for each of the 11 rotational DOF of the model)
%----------------------------------------------------------------------------

Variables MSCx, MSCy, MSCz					% moments on clavicle
Variables MACx, MACy, MACz					% moments on scapula
Variables MGHy, MGHz, MGHyy					% moments on shoulder
Variables MELx								% moment on elbow flexion
Variables MPSy								% moment on forearm pronation

%--------------------------------------------------------------------
% Apply a joint moment to each hinge
%--------------------------------------------------------------------
Torque(Thorax/ClavicleA,    MSCy*Thorax2>)
Torque(ClavicleA/ClavicleB, MSCz*ClavicleA3>)
Torque(ClavicleB/Clavicle,  MSCx*ClavicleB1>)

Torque(Clavicle/ScapulaA,   MACy*Clavicle2>)
Torque(ScapulaA/ScapulaB,   MACz*ScapulaA3>)
Torque(ScapulaB/Scapula,    MACx*ScapulaB1>)

% problem with singularities at GH when q8=0? -- not if we limit range of motion (in das3mex.c) to make gimbal lock unlikely
Torque(Scapula/HumerusA,    MGHy*Scapula2>)
Torque(HumerusA/HumerusB,   MGHz*HumerusA3>)
Torque(HumerusB/Humerus,    MGHyy*HumerusB2>)

Torque(Humerus/Ulna, 	MELx*FlexionAxis>)
Torque(Ulna/Radius, 	MPSy*PronationAxis>)

%--------------------------------------------------------------------
% Determine the YZY rotation angles from thorax to humerus
% Apply the external torques on the axes of this joint coordinate system
%--------------------------------------------------------------------

% This is the symbolic form of an YZY rotation matrix:
%  frames a,b,c,d
%  variables q1,q2,q3
%  simprot(a,b,2,q1);
%  simprot(b,c,3,q2);
%  simprot(c,d,2,q3);
%  a_d
% Result[1,1] = COS(q1)*COS(q2)*COS(q3) - SIN(q1)*SIN(q3)
% Result[1,2] = -SIN(q2)*COS(q1)
% Result[1,3] = SIN(q1)*COS(q3) + SIN(q3)*COS(q1)*COS(q2)
% Result[2,1] = SIN(q2)*COS(q3)
% Result[2,2] = COS(q2)
% Result[2,3] = SIN(q2)*SIN(q3)
% Result[3,1] = -SIN(q3)*COS(q1) - SIN(q1)*COS(q2)*COS(q3)
% Result[3,2] = SIN(q1)*SIN(q2)
% Result[3,3] = COS(q1)*COS(q3) - SIN(q1)*SIN(q3)*COS(q2)

% so we can solve the YZY angles as follows (if we are not in gimbal lock where sin(q2) is zero)
qTH = [														  &
	sym_atan2(Thorax_Humerus[3,2] , (-Thorax_Humerus[1,2])) 	; &		% plane of elevation, between -PI and PI
	acos(Thorax_Humerus[2,2])								; &		% elevation angle, between 0 and PI
	sym_atan2(Thorax_Humerus[2,3] , Thorax_Humerus[2,1])	  &		% humerus internal rotation, between -PI and PI
	]	
encode qTH

% calculate the floating axis (where we need to apply the Z moment)
floatingaxis> = cross(Thorax2>,Humerus2>)
floatingaxis> = floatingaxis>/mag(floatingaxis>)		% normalize to unit vector

% apply an external moment on humerus that is equivalent to the JCS moments that we wanted

variables MTHy, MTHz, MTHyy								% the three JCS moments between thorax and humerus
variables MHx, MHy, MHz									% the XYZ component of the (equivalent) external moment
MH> = vector(Thorax, MHx, MHy, MHz)

% equations that say that JCS moments are projections of the external moment:
equation[1] = dot(MH>, Thorax2>) 		- MTHy
equation[2] = dot(MH>, floatingaxis>) 	- MTHz
equation[3] = dot(MH>, Humerus2>) 		- MTHyy

% now we can solve the external moment (will fail in gimbal lock!):
solve(equation, MHx, MHy, MHz)

% and we apply this external moment to the humerus
Torque(Thorax/Humerus, MH>)		

%-----------------------------------------------------------------------------------------
% Scapulothoracic contact forces on TS and AI, uses the macro in the file contact.a
%-----------------------------------------------------------------------------------------
contact(TS);
contact(AI);

% Also encode those to we have them in the MEX function
F_SCAP = [FxTS , FyTS, FzTS; FxAI , FyAI, FzAI];
encode F_SCAP

%---------------------------------------------------------------------------------------
% Produce coordinates of points for a stick figure

% This is only used for testing in Matlab. 
% The actual 3D visualization is implemented elsewhere and only needs the joint angles q
%---------------------------------------------------------------------------------------

% the first row are the coordinates of the IJ point, which has no mechanical function in the model
% the second row are the coordinates of the SC joint, which are defined in das3mex.c
Stick = [	0 							,  0.1486 						,  -0.0591; &
			par__SCx 					,  par__SCy 				  	,  par__SCz ; &
			dot(P_ThoraxO_AC>,Thorax1>) ,  dot(P_ThoraxO_AC>,Thorax2>)	,  dot(P_ThoraxO_AC>,Thorax3>); &			
			PxTS 						,  PyTS 						,  PzTS; &			
			PxAI 						,  PyAI 						,  PzAI; &					
			dot(P_ThoraxO_GH>,Thorax1>) ,  dot(P_ThoraxO_GH>,Thorax2>)	,  dot(P_ThoraxO_GH>,Thorax3>); &
			dot(P_ThoraxO_UL>,Thorax1>) ,  dot(P_ThoraxO_UL>,Thorax2>)	,  dot(P_ThoraxO_UL>,Thorax3>); &
			dot(P_ThoraxO_RD>,Thorax1>) ,  dot(P_ThoraxO_RD>,Thorax2>)	,  dot(P_ThoraxO_RD>,Thorax3>); &
			dot(P_ThoraxO_HD>,Thorax1>) ,  dot(P_ThoraxO_HD>,Thorax2>)	,  dot(P_ThoraxO_HD>,Thorax3>)]

encode Stick

% Tthe coordinates of TS and AI were already calculated in the contact model,
% so we reuse those here rather than compute them again.

%--------------------------------------------------------------------
% Generate equations of motion
%--------------------------------------------------------------------

Zero = Fr() + FrStar()
% Kane()							% without Kane(), C code is gives same results and size of code is almost identical

%----------------------------------------------------------
% GH net reaction force acting on the Scapula
%----------------------------------------------------------

% compute net reaction force from weight and accelerations of arm segments
% this will be the actual reaction force if the accelerations satisfy the equations of motion (Zero = 0)
FGH> = 		par__Mass_Humerus * (-A_HumerusO_Thorax> - g*Thorax2>)	&
		+	par__Mass_Ulna    * (-A_UlnaO_Thorax>    - g*Thorax2>)	&
		+	par__Mass_Radius  * (-A_RadiusO_Thorax>  - g*Thorax2>)	&
		+	par__Mass_Hand    * (-A_HandO_Thorax>    - g*Thorax2>)	

% express the result as XYZ force components in Scapula reference frame
F_GH = [	dot(FGH>, Scapula1>) ; dot(FGH>, Scapula2>) ; dot(FGH>, Scapula3>) ]

% we need to report F_GH back to the simulation program to help detect shoulder dislocation
% Note that this is the net intersegmental force, not accounting for muscles.  Muscle contributiosn
% are done in das3mex.c
Encode F_GH						

%-----------------------------------------------------------------------------------------------------
% For implicit dynamics: generate symbolic expressions for Zero and Jacobians d(Zero)/d(q,qd,qdd,moments)
%-----------------------------------------------------------------------------------------------------

% matrices needed for Jacobians of ZZ with respect to q, qdot, qdotdot, and moments
q   = [q1  , q2  , q3  , q4  , q5  , q6  , q7  , q8  , q9  , q10  , q11  ]
qd  = [q1' , q2' , q3' , q4' , q5' , q6' , q7' , q8' , q9' , q10' , q11' ]
qdd = [q1'', q2'', q3'', q4'', q5'', q6'', q7'', q8'', q9'', q10'', q11'']
mom = [MSCy, MSCz, MSCx, MACy, MACz, MACx, MGHy, MGHz, MGHyy, MELx, MPSy];

% generate expressions for implicit equation of motion and its Jacobians
dz_dq = ZEE(D(Zero,q))
dz_dqd = ZEE(D(Zero,qd))
dz_dqdd = ZEE(D(Zero,qdd))
dz_dmom = ZEE(D(Zero,mom))
Encode Zero, dz_dq, dz_dqd, dz_dqdd, dz_dmom

% write all Encoded expressions to C file
Code Algebraic() das3_al_raw.c

EXIT

