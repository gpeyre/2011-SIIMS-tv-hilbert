function display_orientations(Ori,Dx,n,leng,Image)

if (nargin==5)
    imageplot(Image);
end
if (nargin<4)
    leng=1;
end

for i=1:size(Ori,1)
    for j=1:size(Ori,2)
        x = Ori(i,j,1);
        y = Ori(i,j,2);
        if (x~=0 || y ~=0)
            x_c = 1+(j-1)*Dx;
            y_c = 1+(i-1)*Dx;
            X = [x_c x_c+leng*x];
            Y = [y_c y_c+leng*y];
            line(X, Y,'color','k','LineWidth',2);
        end
    end
end
