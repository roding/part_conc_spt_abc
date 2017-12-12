function [  model, ...
            number_of_components, ...
            number_of_abc_samples, ...
            m, ...
            c, ...
            az, ...
            dist, ...
            w, ...
            weighting_scheme, ...
            gamma, ...
            t_exec, ...
            number_of_iterations, ...
            number_of_simulations] = read_output(file_path)

    file_string = fileread(file_path);

    model = read_key(file_string, 'model', 'string');
	number_of_components = read_key(file_string, 'number_of_components', 'scalar');
	number_of_abc_samples = read_key(file_string, 'number_of_abc_samples', 'scalar');
	m = reshape(read_key(file_string, 'm', 'array'), [number_of_components, number_of_abc_samples]);
	c = reshape(read_key(file_string, 'c', 'array'), [number_of_components, number_of_abc_samples]);
    
    az = read_key(file_string, 'az', 'array');
    if numel(az) == number_of_components * number_of_abc_samples
        az = reshape(az, [number_of_components, number_of_abc_samples]);
    else
        az = az(:);
        az = repmat(az, [3 1]);
        az = reshape(az, [number_of_components, number_of_abc_samples]);
    end
	
    dist = read_key(file_string, 'dist', 'array');
	w = read_key(file_string, 'w', 'array');
	weighting_scheme = read_key(file_string, 'weighting_scheme', 'string');
	gamma = read_key(file_string, 'gamma', 'scalar');
	t_exec = read_key(file_string, 't_exec', 'scalar');
	number_of_iterations = read_key(file_string, 'number_of_iterations', 'scalar');
	number_of_simulations = read_key(file_string, 'number_of_simulations', 'scalar');
    
end