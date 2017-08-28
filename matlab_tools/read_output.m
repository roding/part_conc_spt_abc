function [  diagnostic_diffusion_coefficient_ratio, ...
            t, ...
            msd, ...
            msd_x, ...
            msd_y, ...
            msd_z, ...
            execution_time] = read_output(file_path)

file_string = fileread(file_path);

diagnostic_diffusion_coefficient_ratio = read_key(file_string, 'diagnostic_diffusion_coefficient_ratio', 'scalar');
t = read_key(file_string, 'time', 'array');
msd = read_key(file_string, 'mean_square_displacement', 'array');
msd_x = read_key(file_string, 'mean_square_displacement_x', 'array');
msd_y = read_key(file_string, 'mean_square_displacement_y', 'array');
msd_z = read_key(file_string, 'mean_square_displacement_z', 'array');
execution_time = read_key(file_string, 'execution_time', 'scalar');

end

