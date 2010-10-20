n=1000;
s=5;
t=4;
A=gallery('lehmer',n);
r=rand(n,1);
%[T,V]=lanczos_s(A,r,s,t);
[T,V]=ca_lanczos(A,r,s,t);
%[T,V]=ca_lanczos_second(A,r,s,t);
%[T,V]=ca_lanczos_overlap(A,r,s,t);
%for i=1:10
%    [T,V]=ca_lanczos(A,r,s,t);
%    r = V(:,s*t+1);
%end
[T2,V2]=lanczos(A,r,s*t);
%[T3,V3]=sStepLanczos(A,r,s,t);
%E=eig(T(1:s,1:s));
%Es=sort(E,'descend');
EA=eig(A);
EAs=sort(EA,'descend');
ET=eig(T);
ETs=sort(real(ET),'descend');
ET2=eig(T2);
ET2s=sort(ET2,'descend');
%ET3=eig(T3);
%ET3s=sort(ET3,'descend');
%Es(1:s)
EAs(1:s*t)
ETs(1:size(T,1))
ET2s(1:size(T2,1))
%ET3s(1:s*t)