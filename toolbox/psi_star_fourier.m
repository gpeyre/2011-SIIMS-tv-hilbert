function res = psi_star_fourier(I)
  res = ifft2(I)*sqrt(numel(I));