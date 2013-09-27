function Text=image_texture(n,freq,a)
% a largeur parabole

t = 1:n;
[X,Y] = meshgrid(t,t);

Text = cos(pi*freq*(((X-n/2).*(X-n/2))/a+(Y-4*n/5)));