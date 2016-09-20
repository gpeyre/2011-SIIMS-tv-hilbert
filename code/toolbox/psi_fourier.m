function res = psi_fourier(I)
  res = fft2(I)/sqrt(numel(I));

  