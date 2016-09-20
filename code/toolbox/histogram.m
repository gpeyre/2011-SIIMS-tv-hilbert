function histogram(I,valeurs)

if (nargin<2)
	mini = min(I(:));
    maxi = max(I(:));
    valeurs = mini:(maxi-mini)/100:maxi;
end

I = I(:);

Histo= zeros(size(valeurs,2));
for i=1:(size(valeurs,2)-1)
    Histo(i) = sum((I>valeurs(i)).*(I<valeurs(i+1)));
end

plot(valeurs,Histo);