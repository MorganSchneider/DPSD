function vx = get_velocity_axis(va, M)
% GET_VELOCITY_AXIS Get the unambiguous velocity range
%   GET_VELOCITY_AXIS()

vx = ((-ceil(M/2)+1):1:floor(M/2)) * 2 * va / M;
% From [+va to -va) for plotting
vx = fliplr(vx);

end