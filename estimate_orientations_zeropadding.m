function Orientations = estimate_orientations_zeropadding(v,mask,q,Dx,tolerance)
% calcule frequence instantanée
% après zoom par zéro padding pour obtenir frequence non-entiere
PM = compute_patch(mask, q, Dx);
Pv = compute_patch(v, q, Dx);
if (~exist('tolerance'))
    tolerance=0;
end

tau = 2;

Orientations = zeros(size(Pv,3),size(Pv,4),2);

for p1=1:size(Pv,3)
    for p2=1:size(Pv,4)
        S = Pv(:,:,p1,p2);
        multi = 1;
        S2 = zeros(multi*q,multi*q);
        S2((multi-1)*q/2+1:(multi+1)*q/2,(multi-1)*q/2+1:(multi+1)*q/2) = S;

        Tmp = abs(fftshift(fft2(S2)));
        Tmp(end/2+1-tau*multi:end/2+1+tau*multi,end/2+1-tau*multi:end/2+1+tau*multi) = 0;
        [f1,f2,maxi] = instantaneous_frequency(Tmp);
        f1 = f1/multi; f2 = f2/multi;

        Orientations(p1,p2,1) = 0;
        Orientations(p1,p2,2) = 0;
         if (sum(sum(PM(:,:,p1,p2)))>.6*q*q) %dans le cas de l'inpainting sous-fenetre plus de 60% dans masque
             continue;
         end
        f_mean = mean(Tmp(:));
        if (maxi>tolerance*f_mean)
            Orientations(p1,p2,1) = f1;
            Orientations(p1,p2,2) = -f2;
        end
    end
end