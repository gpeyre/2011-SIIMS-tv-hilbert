function l1 = L1_curv(M,options)
% norme L1 de la decomposition curvelet
l1 = 0;
x = perform_curvelet_transform(M,options);
for i=1:size(x,1)
    for j=1:size(x,2)
        tmp = x{i,j};
        for k=1:size(tmp,1)
            for l=1:size(tmp,2)
                l1 = l1 + sum(sum(abs(tmp{k,l})));
            end
        end
    end
end