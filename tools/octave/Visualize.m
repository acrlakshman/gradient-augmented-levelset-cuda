% Visualize fields.
%
% Fields that will be read are:
%   phi
%   psi
%   velocity
%
%   Data formats: x y z <scalar> (or) x y z <cmpt-1> <cmpt-2> ... <cmpt-N>

function Visualize(path, dim)

    phi_str = strcat(path, '/phi');
    psi_str = strcat(path, '/psi');
    velocity_str = strcat(path, '/velocity');

    if exist(velocity_str, 'file') == 2
        ShowVector(velocity_str, dim);
    end

end

function ShowVector(path, dim)

    fid = fopen(path);
    ts = textscan(fid, '%f\t%f\t%f\t%f\t%f\t%f');
    fclose(fid);

    x = ts{1};
    y = ts{2};
    z = ts{3};

    u = ts{4};
    v = ts{5};
    w = ts{6};

    if (dim == 2)
        quiver(x, y, u, v);
    end

end
