% just some code snippets to help me understand why the long-term BPR
% looks so weird in comparison to other stations

[ttemp,ia,ib]=intersect(round(A(1,:),6),round(tF,6));

tinv=ttemp-ttemp(1); tinv=tinv/tinv(end);
tinv=tinv; Ainv=A(2,ia)'; Finv=F{11}(ib);
G=[ones(size(tinv)), tinv, Finv];
m=inv(G'*G)*G'*Ainv;
Amod=G*m;