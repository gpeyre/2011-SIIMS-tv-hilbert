function n = non_convex_norm(F, Chi)

F = abs(F);
tmp = F - Chi;

n = sum(tmp(:).*tmp(:));