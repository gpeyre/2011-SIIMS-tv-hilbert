function hv = test_h2(v,a,b)
v = clamp(v,-1,1);
hv=sign(v-b).*(1-(1-abs(v-b)).^a);