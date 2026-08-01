function [omega_n_Hz, eigenvectors] = compute_omega_n(A)
    % Compute the eigenvalues and eigenvectors of the state-space matrix A
    [eigenvectors, eigenvalues] = eig(A);
    
    % Extract the eigenvalues from the diagonal matrix
    eigenvalues = diag(eigenvalues);
    
    % Compute the natural frequencies (in rad/s)
    natural_frequencies_rad_per_s = abs(eigenvalues);
    
    % Convert the natural frequencies to Hertz (Hz)
    omega_n_Hz = natural_frequencies_rad_per_s / (2 * pi);
    
    % Display the natural frequencies in Hz
    disp('The natural frequencies (Hz) of the system are:');
    disp(omega_n_Hz);
    
    % Display the eigenvectors (mode shapes)
    disp('The eigenvectors (mode shapes) of the system are:');
    disp(eigenvectors);
end
