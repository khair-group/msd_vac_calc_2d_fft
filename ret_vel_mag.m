function [vel_mag] = ret_vel_mag(xvel,yvel)
%Retruns the magnitude of the two-dimensional velocity vector
len=size(xvel);
size_list=len(1);
xpos=zeros(size_list,1);
ypos=zeros(size_list,1);
vel_mag=zeros(size_list,1);

for i=1:(size_list)
    vel_mag(i)=sqrt(xvel(i).^2+yvel(i).^2);
end

