% create rotation matrix r, used in calculating the position, jacobians, and hessians
function [ r ] = r3( theta,axis,d )
% rotate [theta] (radians) about
% 'x' 'y' or 'z' [axis], with
% order derivative [d], 0 or not given for none (pure rotation)
if nargin<3;
    d=0;
end
r=eye(3);
x=[cos(theta),-sin(theta);sin(theta),cos(theta)];
dx=[-sin(theta),-cos(theta);cos(theta),-sin(theta)];
ddx=[-cos(theta),sin(theta);-sin(theta),-cos(theta)];
if axis == 'x';
    if d==1
        r(2:3,2:3)=dx;r(1,1)=0;
    elseif d==2
        r(2:3,2:3)=ddx;r(1,1)=0;
    else
        r(2:3,2:3)=x;
    end
elseif axis == 'y';
    if d==1
        r(1,1)=dx(2,2);r(3,3)=dx(1,1);r(3,1)=dx(1,2);r(1,3)=dx(2,1);r(2,2)=0;
    elseif d==2
        r(1,1)=ddx(2,2);r(3,3)=ddx(1,1);r(3,1)=ddx(1,2);r(1,3)=ddx(2,1);r(2,2)=0;
    else
        r(1,1)=x(2,2);r(3,3)=x(1,1);r(3,1)=x(1,2);r(1,3)=x(2,1);
    end
elseif axis == 'z';
    if d==1
        r(1:2,1:2)=dx;r(3,3)=0;
    elseif d==2
        r(1:2,1:2)=ddx;r(3,3)=0;
    else
        r(1:2,1:2)=x;
    end;
end;
end