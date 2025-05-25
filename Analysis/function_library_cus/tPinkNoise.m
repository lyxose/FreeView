function truncatedNoise = tPinkNoise(N, seed, contrast, show)
% generate a 2d 1/f noise texture in gray level, truncated at 2*sd.
% -----Parameter-----
% N: Size of the image (N x N)
    
    if nargin<3
        contrast=1;
    end
    
    if nargin<4
        show=false;
    end
    currentSeed = rng;  % save current rng seed state
    rng(seed);

    % Create a grid of frequencies
    [u, v] = meshgrid(-N/2:N/2-1, -N/2:N/2-1);  % Frequency grid
    freq = sqrt(u.^2 + v.^2);  % Radial frequency
    
    % Avoid division by zero for the DC component
    freq(freq == 0) = 1;
    
    % Generate random phases

    randomPhase = exp(1i * 2 * pi * rand(N, N));  % Complex random phases
    
    % Scale the amplitude by 1/f
    amplitude = 1 ./ freq;  % 1/f scaling
    amplitude = amplitude / max(amplitude(:));  % Normalize amplitude
    
    
    % Create the frequency-domain representation of the noise
    frequencyDomainNoise = amplitude .* randomPhase;
    
    % Transform back to the spatial domain
    spatialNoise = real(ifft2(ifftshift(frequencyDomainNoise)));
    
    % Calculate standard noise
    stdNoise = std(spatialNoise(:));
    
    % Truncat by 2*sd
    truncatedNoise = spatialNoise;
    truncatedNoise(abs(truncatedNoise) > 2*stdNoise) = 2*sign(truncatedNoise(abs(truncatedNoise) > 2*stdNoise))*stdNoise;
    
    % Normalize the spatial noise to range [0, 1]
    truncatedNoise = (truncatedNoise - min(truncatedNoise(:))) / (max(truncatedNoise(:)) - min(truncatedNoise(:)));
    
    % apply contrast
    truncatedNoise = contrast * (truncatedNoise-0.5) + 0.5;
    
    % Display the 1/f noise image
    if show
        figure;
        imshow(truncatedNoise, []);
        title('1/f Noise Image');
    end
    rng(currentSeed);  % reset rng as before

end
