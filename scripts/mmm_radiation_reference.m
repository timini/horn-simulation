% Authored adapter for the separately installed, pinned MMM Toolbox.
% No upstream acoustic implementation is copied into this repository.
addpath('/reference/src');
fprintf('Octave %s\n', version());
for name = {'MMM_init', 'MMM_calculateMatrices', 'MMM_ASradiatedPressure', 'MMM_ASbaffledradzmatrix'}
  assert(strcmp(which(name{1}), ['/reference/src/' name{1} '.m']));
end
p = jsondecode(fileread('/study/protocol.json'));
load('/study/MMM_besselzeros.mat', 'bz');
frequencies = p.frequencies_hz(:)';
coords = dlmread('/study/coordinates_1000.csv', ',');
maximum_modes = 64;
wavenumbers = 2*pi*frequencies/p.air.c;
% Direct integration, no precomputed interpolation file or HF approximation.
radiation = MMM_ASbaffledradzmatrix(wavenumbers, p.air.rho, p.air.c, ...
                                  pi*p.geometry.mouth_radius^2, maximum_modes, bz, false, false);
for case_index = 1:numel(p.cases)
  item = p.cases(case_index);
  coords = dlmread(sprintf('/study/coordinates_%d.csv', item.sections), ',');
  output = zeros(numel(frequencies), 7);
  for frequency_index = 1:numel(frequencies)
    data = MMM_init(frequencies(frequency_index), item.modes, coords, 'axi', p.air.rho, p.air.c);
    data.Zrad = radiation(1:item.modes, 1:item.modes, frequency_index);
    data = MMM_calculateMatrices(data, false);
    data.nIntegrationPoints = 201;
    data = MMM_ASradiatedPressure(data, [0, p.observer_distance_m], false);
    pressure_201 = data.pRad;
    data.nIntegrationPoints = 401;
    data = MMM_ASradiatedPressure(data, [0, p.observer_distance_m], false);
    output(frequency_index,:) = [frequencies(frequency_index), ...
      real(data.Z00*data.St), imag(data.Z00*data.St), ...
      real(pressure_201), imag(pressure_201), real(data.pRad), imag(data.pRad)];
    assert(all(isfinite(output(frequency_index,:))));
    fprintf('Completed %s %.6f Hz\n', item.id, frequencies(frequency_index));
    fflush(stdout);
  end
  dlmwrite(sprintf('/study/%s.csv', item.id), output, 'precision', '%.17g');
end
