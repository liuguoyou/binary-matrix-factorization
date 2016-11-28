function [u,v]=sssvd(X,lambda)
  %
  % coordinate descent minimization for
  %
  % ||X-uv||_1 + lu||u||_1 + lv||v||_1
  %
  [m,n]=size(X);
  lu=sqrt(m)*lambda;
  lv=sqrt(n)*lambda;
  u=rand(m,1);
  %u=u*(1/norm(u));
  v=rand(n,1);
  %v=v*(1/norm(v));
  f=norm(X-u*v',1)
  for J=1:10
    uOld = u;
    vOld = v;
    fOld = f;
      % update u
      % update v
      for i=1:m
	  Xi=[X(i,:) 0];
	  u(i)=mincoso(Xi,[v; lu]');
      end 
      for j=1:n
	  Xj=[X(:,j)' 0];
	  v(j)=mincoso(Xj,[u; lv]');
      end 
    nu=norm(u);
    ndu=norm(u-uOld);
    nv=norm(v);
    ndv=norm(v-vOld);
    f=norm(X-u*v',1)+lu*sum(abs(u))+lv*sum(abs(v));
    df=fOld-f;
    fprintf('i=%d\tf=%6f\tnu=%6.0f\tnv=%6.0f\tdu=%6.0f\tdv=%6.0f\tdf=%6.0f\n',i,f,nu,nv,ndu,ndv,df);

  end
end
