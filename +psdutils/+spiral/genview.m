function R = genview(shotn, nshots, spiraln, nspirals, ...
                rtype, rlist)
% shotn = shot index
% nshots = total number of shots
% spiraln = spiral index
% nspirals = total number of spirals
% rtype = 3D transformation type (sos = stack of spirals,
%   r1ax = single axis rotations, r2ax = dual-axis rotations,
%   fromlist = from list of parameters)
% rlist = list of rotations/translations if rtype is 'fromlist'

    PHI = (3 - sqrt(5)) / 2; % 1D golden ratio
    phi1 = 0.4656; % 2D golden ratio 1
    phi2 = 0.6823; % 2D golden ratio 2
    
    switch lower(rtype)
        case 'sos'
            rx = 0;
            ry = 0;
            rz = 2*pi * shotn / PHI;
            dz = (-1)^(spiraln)/nspirals * 2*floor((spiraln + 1) / 2);
            dz = dz + (-1)^(shotn)/(nshots*nspirals) * 2*floor((shotn + 1) / 2);
        case 'r1ax'
            rx = 2*pi*(shotn*nspirals + spiraln)*PHI;
            ry = 0;
            rz = 0;
            dz = 1;
        case 'r2ax'
            rx = 2*acos(mod((shotn*nspirals + spiraln)*phi1,1)); % polar angle
            ry = 0;
            rz = 2*pi*mod((shotn*nspirals + spiraln)*phi2,1); % azimuthal angle
            dz = 1;
        case 'fromlist'
            m = shotn*nspirals + spiraln;
            rx = rlist(m,1);
            ry = rlist(m,2);
            rz = rlist(m,3);
            dz = rlist(m,4);
    end
    
    Tdz = eye(3);
    Tdz(3,3) = dz;
    R = rmat3D('zyx',[rz,ry,rx])*Tdz;

end

function R = rmat3D(axis,angle)
% axis = order of axises to perform rotation
% angle = euler angles (rad) to perform rotation
%
% R = rotation matrix
%

    % Recurse for multiple angles
    if length(axis) > 1
        R = rmat3D(axis(1:end-1),angle(1:end-1))* ...
            rmat3D(axis(end),angle(end));
    else
        switch lower(axis)
            case 'x' % X rotation matrix
                R = [
                    1,              0,              0;
                    0,              cos(angle),     -sin(angle);
                    0,              sin(angle),     cos(angle)
                    ];
            case 'y' % Y rotation matrix
                R = [
                    cos(angle),     0,              sin(angle);
                    0,              1,              0;
                    -sin(angle),    0,              cos(angle)
                    ];
            case 'z' % Z rotation matrix
                R = [
                    cos(angle),     -sin(angle),    0;
                    sin(angle),     cos(angle),     0;
                    0,              0,              1
                    ];
            otherwise
                error('invalid axis %s (argument 1)', axis);
        end
    end
end