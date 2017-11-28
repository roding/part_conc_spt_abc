function [  ax, ...
            ay, ...
            number_of_frames, ...
            deltat, ...
            K, ...
            DE] = read_data(file_path)

    file_string = fileread(file_path);

    ax = read_key(file_string, 'ax', 'scalar');
    ay = read_key(file_string, 'ax', 'scalar');
    number_of_frames = read_key(file_string, 'number_of_frames', 'array');
    deltat = read_key(file_string, 'deltat', 'scalar');
    K = read_key(file_string, 'K', 'array');
    DE = read_key(file_string, 'DE', 'array');
    
end