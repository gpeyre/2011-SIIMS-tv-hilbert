function Chi = non_convex_weight(Orient,Amplit,Weight, q)
% construit les points chi pour la norme non convexe 
% Orient: frequences instantanées (estime_orientation)
% Amplit: amplitude (estimate_amplitude)
% Weight: poids associés au fourier local
% q: taille des fenetres
Chi = ones(q,q,size(Orient,1),size(Orient,2));

for i=1:size(Orient,1)
    for j=1:size(Orient,2)
        x_or = Orient(i,j,1);
        y_or = Orient(i,j,2);
        Chi(:,:,i,j) = zeros(q,q);

        if (x_or~=0 || y_or ~=0)
            freq1 = x_or;
            freq2 = -y_or;


            t = 1:q;
            [Y,X] = meshgrid(t,t);
            Z = 2*pi/q*(freq1*X + freq2*Y);
            Sinu = cos(Z)/q;
            Sinu = Weight(:,:,i,j).*Sinu;

            Chi(:,:,i,j) = Amplit(i,j)*abs(fftshift(fft2(Sinu)));
        end
    end
end


