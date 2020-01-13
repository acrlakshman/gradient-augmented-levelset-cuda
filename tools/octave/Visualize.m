% Visualize fields.% % Fields that will be read are : % phi % psi % velocity % %
                                                      Data formats
    : x y z<scalar>(or) x y z<cmpt - 1><cmpt - 2>...<cmpt - N>

      function Visualize(path, dim, nx, ny, nz)

          num_dirs = length(ls(path));

fprintf('num_dirs = %d\n', num_dirs);
    for
        d = 1 : num_dirs

                    phi_str = strcat(path, '/', num2str(d), '/phi');
    psi_str = strcat(path, '/', num2str(d), '/psi');
    velocity_str = strcat(path, '/', num2str(d), '/velocity');

    if
        exist(phi_str, 'file') == 2 ShowScalar(phi_str, dim, nx, ny, nz);
    end

        if exist(velocity_str, 'file') == 2 ShowVector(velocity_str, dim);
    end

        fprintf('d = %d\n', d);

    end

end

function ShowScalar(path, dim, nx, ny, nz)

    fid = fopen(path);
    ts = textscan(fid, '%f\t%f\t%f\t%f');
    fclose(fid);

    x = reshape(ts{ 1 }, nx, ny);
    y = reshape(ts{ 2 }, nx, ny);
    z = ts{ 3 };

    a = reshape(ts{ 4 }, nx, ny);

    if (dim == 2)
        figure(1) contour(x, y, a, [0 0]);
    axis square;
    % surfc(x, y, a);
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
        figure(2);
            quiver(x, y, u, v);
            axis square;
            end

                end
