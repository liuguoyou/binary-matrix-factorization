%
% given vectors x and v of the same size,
% find u which minimizes \sum_i |x_i - u*v_i|
   % this requires the following steps:
   % set J={ i:v_i != 0 }, sort J
% set w_i = |v_i|, z_i = x_i/v_i, i \in J
% choose u = z_{j'} where j' is the smallest index in J
% such that \sum_{j <= j'}w_j >= \sum_{j > j'}w_j 
function u=mincoso(x,v)
	 J=find(v);
	 x0=x(J);
	 v0=v(J);
	 z0=x0./v0;
	 [z0,J]=sort(z0); 
	 v0=v0(J); % synchronize indexes of z0 and v0
	 w0=abs(v0);
	 sleft=0;
	 sright=sum(w0);
	 u=z0(end);
	 for i=1:length(J)
	     sleft=sleft+w0(i);
	     sright=sright-w0(i);
%	     fprintf('i=%d\tw=%f\tz=%f\tsl=%f\tsr=%f\n',i,w0(i),z0(i),sleft,sright);
	     if sleft >= sright
		u=z0(i);
		break
	     end
	 end
end
	 
	 
	 
	 
