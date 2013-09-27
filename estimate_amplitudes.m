function A = estimate_amplitudes(F_estim, Chi)
% "essaie" d'estimer l'amplitude des frequences présentes dans Chi

A = zeros(size(Chi,3),size(Chi,4));
for i=1:size(Chi,3)
    for j=1:size(Chi,4)
       if (norm(Chi(:,:,i,j))~=0)
            tmp = abs(F_estim(:,:,i,j)).*Chi(:,:,i,j);
            tmp2 = norm(Chi(:,:,i,j),'fro');
            A(i,j) = sum(tmp(:))/(tmp2*tmp2);
       % A(i,j) = norm(abs(F_estim(:,:,i,j)),'fro')/tmp2;
       % A(i,j) = max(max(abs(F_estim(:,:,i,j))))/ max(max(Chi(:,:,i,j)));
       end
    end
end