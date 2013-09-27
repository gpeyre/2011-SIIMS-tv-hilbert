function dhv = test_diff_h(v,a,b)
%dhv = a./(pi*(1+a*a*(v+b).*(v+b))); %creneaux
v = clamp(v);
 if (a<1)
     epsilon = .1;
 else 
     epsilon=0;
 end


dhv = a*max(e,abs(v-b)).^(a-1);
