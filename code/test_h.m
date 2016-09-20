function hv = test_h(v,a,b)
hv=sign(v-b).*abs(v-b).^a;