function Orientations = estimate_orientations(F,mask,q,Dx,tolerance)
% calcule frequence instantanée
% F: transformée fourier local
% mask: masque eventuel (si pas inpainting mask = 0)
% q: taille fenetres fourier local
% Dx: espacement fenetres
% tolerance: presence/asbence texture

P = compute_patch(mask, q, Dx);
if (~exist('tolerance'))
    tolerance=0;
end

F = abs(F);        
k = 1;
F(end/2+1-k:end/2+1+k,end/2+1-k:end/2+1+k,:,:) = 0;

Orientations = zeros(size(F,3),size(F,4),2);

for p1=1:size(F,3)
    for p2=1:size(F,4)
        Tmp = F(:,:,p1,p2);
        [f1,f2,maxi] = instantaneous_frequency(Tmp);
        Orientations(p1,p2,1) = 0;
        Orientations(p1,p2,2) = 0;
         if (sum(sum(P(:,:,p1,p2)))>.4*q*q) %dans le cas de l'inpainting sous-fenetre plus de 40% dans masque (en fait c'est con, on perd de l'info sur les bords a cause de la fenetre)
             continue;
         end
        f_mean = mean(Tmp(:));
        if (maxi>tolerance*f_mean) %texture signifivative présente ?
            Orientations(p1,p2,1) = f1;
            Orientations(p1,p2,2) = -f2;
        end
    end
end


% %plus une passe test pour affiner
% 
% for p1=1:size(Pv,3)
%     for p2=1:size(Pv,4)
%         x_or = Orientations(p1,p2,1);
%         y_or = Orientations(p1,p2,2);
%         if (x_or == 0 && y_or == 0)
%             continue;
%         end
%         S = Pv(:,:,p1,p2);
%         Tmp = abs(fftshift(fft2(S)));
% 
%         t = 1:q;
%         [Y,X] = meshgrid(t,t);
% 
%         freq1 = x_or;
%         freq2 = -y_or;
%         Sinu = cos(2*pi/q*(freq1*X + freq2*Y));
%         Tmp2 = abs(fftshift(fft2(Sinu)));
%         tmp = abs(Tmp.*Tmp2);
%         tmp2 = norm(Tmp2,'fro');
%         Chi_ij = sum(tmp(:))/(tmp2*tmp2)*Tmp2;
%         tmp3=norm(Tmp-Chi_ij,'fro');
% 
%         dx = 1/8; dy = 1/8;
%         mieux=true;
%         change = false;
%         while(mieux)
%             mieux = false;epsilon_x = 0;epsilon_y=0;
%             for eps_x=-1:1
%                 for eps_y=-1:1
%                     if (eps_x==0 && eps_y==0)
%                         continue;
%                     end
%                     x_or_new = x_or+eps_x*dx;
%                     y_or_new = y_or+eps_y*dy;
%                     freq1 = x_or_new;
%                     freq2 = -y_or_new;
%                     Sinu = cos(2*pi/q*(freq1*X + freq2*Y))/q;
%                     Tmp2 = abs(fftshift(fft2(Sinu)));
%                     tmp = abs(Tmp.*Tmp2);
%                     tmp2 = norm(Tmp2,'fro');
%                     Chi_ij = sum(tmp(:))/(tmp2*tmp2)*Tmp2;
%                     tmp4=norm(Tmp-Chi_ij,'fro');
%                     if (tmp4<tmp3)
%                         tmp3 = tmp4;
%                         epsilon_x = eps_x;
%                         epsilon_y = eps_y;
%                         mieux = true;
%                          %disp(['mieux ex=' num2str(eps_x) ' ey=' num2str(eps_y) ' p1=' num2str(p1) ' p2=' num2str(p2)] );
%                     end
%                 end
%             end
%             x_or = x_or + epsilon_x*dx;
%             y_or = y_or + epsilon_y*dy;
%         end
%         Orientations(p1,p2,1) = x_or;
%         Orientations(p1,p2,2) = y_or;
%     end
% end