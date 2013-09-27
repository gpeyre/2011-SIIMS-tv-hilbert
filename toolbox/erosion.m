function E=erosion(M,p)
E = M;
if (p>0)
    if (sum(M(:))==0)
        return;
    end
    for i=2:size(M,1)-1
        for j=2:size(M,2)-1
            if (M(i,j)==0)
                E(i-1,j)=0;
                E(i+1,j)=0;
                E(i,j-1)=0;
                E(i,j+1)=0;
            end
        end
    end
    if (p>1)
        E = erosion(E,p-1);
    end

end