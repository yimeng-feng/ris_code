function y = array_response(a1,a2,N,dx,lambda,p0)
y=zeros(N,1);
for m= 0:p0-1
    for n= 0:N/p0-1
        y(m*N/p0+n+1) = exp( 1i*(2*pi*dx/lambda)* ( m*sin(a1)*sin(a2) + n*cos(a2) ) );
    end
end 
end