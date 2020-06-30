function make_osimm(filename,angles,time)
% Creates motion file for OpenSim
%
% Inputs:
% filename: the name of the motion file, without the extension
% angles: an 11xn or nx11 matrix of angles in radians, n the number of
% frames
% time (optional): a 1xn or nx1 vector of time values. If this is not
% provided, the timestep is assumed to be 0.01s.
%
% Corrected: angles now in radians instead of degrees
% Dimitra Blana, February 2012

[nrows,ncolumns]=size(angles);
if ncolumns~=11
    angles=angles';
    nrows = ncolumns;
end

if nargin <3
    time = 0.01:0.01:0.01*nrows;
end

if size(time,2)~=1, time = time'; end
if size(time,1)~=nrows
    errordlg('The time vector does not have the same length as the angle data.','Dimension error');
    return;
end
% data = [time angles*180/pi];  % in degrees for the Opensim model <--
data = [time angles];  

% create motion file
% the header of the motion file is:
%
% <motion name>
% nRows=x
% nColumns=y
% endheader
% time SC_y SC_z SC_x AC_y AC_z AC_x GH_y1 GH_z GH_y2 EL_x PS_y
%
fid = fopen([filename '.sto'],'wt');
fprintf(fid,'%s\n',filename);
fprintf(fid,'%s%i\n','nRows=',nrows);
fprintf(fid,'%s\n','nColumns=12');
fprintf(fid,'%s\n','endheader');
fprintf(fid,'%s\n','time  SC_y  SC_z SC_x  AC_y  AC_z  AC_x  GH_y  GH_z  GH_yy  EL_x  PS_y');
fprintf(fid,'%f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n',data');
fclose(fid);

