function res = imnoiseSNR(A, snr, noiseStrength)
%
% Applies white Gaussian noise with dB intensity given by 'snr'.
% Alternatively, the noise variance can be given as input in
% 'noiseStrength'.
% I is considered grey-level (dimensionality one)
%
% Examples
%       A2 = imnoiseSNR(A, 2);
%        will bury A in 2db noise.
%       A9 = imnoiseSNR(A, [], 9);
%        will bury A in noise of variance = 9 intensity units (std = 3).
%
% G.Sfikas 21 Jul 2007
% Update 7 Feb 2008
% Update 18 Jan 2009 - peut traiter des images RGB
% Update 1 Feb 2009 - peut traiter matrices egalement "double" ou "uint8".
temp2 = double(A(:)');

signalStrength = inv(numel(temp2)) * (temp2*temp2');
if nargin < 3
    noiseStrength = signalStrength / (10^(snr/10));
else
    snr = (log(signalStrength)-log(noiseStrength))*(10/log(10));
end
temp2 = temp2 + sqrt(noiseStrength) * randn(1, numel(temp2));
fprintf('Added white additive Gaussian noise of %f dBel intensity (st.dev. = %f)\n', snr, sqrt(noiseStrength)); 

if isfloat(A) == 0
    res = uint8(reshape(temp2, size(A)));
else
    res = reshape(temp2, size(A));
end
return;