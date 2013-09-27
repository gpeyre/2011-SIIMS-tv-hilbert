function dhv = test_diff_h2(v,a,b)
v = clamp(v,-1,1);
dhv = a*(1-abs(v-b)).^(a-1);