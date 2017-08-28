function [  distribution_class, ...
            number_of_components, ...
            number_of_abc_samples, ...
            m, ...
            s, ...
            c, ...
            az, ...
            dist, ...
            w, ...
            epsilon] = read_output(file_path)

    file_string = fileread(file_path);

    distribution_class = read_key(file_string, 'distribution_class', 'string');
    number_of_components = read_key(file_string, 'number_of_components', 'scalar');
    number_of_abc_samples = read_key(file_string, 'number_of_abc_samples', 'scalar');
    
    m = read_key(file_string, 'm', 'array');
    m = reshape(m, [number_of_components, number_of_abc_samples]);
    s = read_key(file_string, 's', 'array');
    s = reshape(s, [number_of_components, number_of_abc_samples]);
    c = read_key(file_string, 'c', 'array');
    c = reshape(c, [number_of_components, number_of_abc_samples]);    
    
    az = read_key(file_string, 'az', 'array');
    dist = read_key(file_string, 'dist', 'array');
    w = read_key(file_string, 'w', 'array');
    epsilon = read_key(file_string, 'epsilon', 'scalar');

end