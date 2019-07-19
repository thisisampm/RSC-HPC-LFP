function energies = waveform_energy(samples)
%calculates the energy of every spike
%
%samples is output by nauralyns and is a 32 waveform features by 4 tetrodes
%by n spike events 3d matrix

%for each spike
%square each value in the waveform
%sum the squares
%divide by the number of waveform values

energies = squeeze(sqrt(sum(samples.^2))./size(samples,1))';