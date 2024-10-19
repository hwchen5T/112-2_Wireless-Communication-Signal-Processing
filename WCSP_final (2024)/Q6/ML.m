function [ symbol_ML ] = ML(Da_Str,y,H,Q,R,NPW)
    x_hat = (1/sqrt(2))*[-1+1i,-1-1i,1+1i,1-1i];
    K=16;
    d=0;
    kvec=zeros(length(x_hat),1);
    z=Q'*y;
    for ii = Da_Str:-1:1
        [~,N]=size(kvec);
        kvec_matrix=repmat(x_hat,N,1);
        kvec_matrix=reshape(kvec_matrix,[N*length(x_hat),1]).';
        kvec=[kvec,kvec,kvec,kvec]+[zeros(ii-1,N*length(x_hat));
        kvec_matrix;
        zeros(Da_Str-ii,N*length(x_hat))]; 
        d=[d,d,d,d];
        [d,index] = sort(d+abs(z(ii)-R(ii,:)*kvec).^2);             
        if (ii==4)
           kvec=kvec(:,index(1:length(x_hat)));                     
           d=d(:,1:length(x_hat));                                  
        elseif (ii==1)                                              
           symbol_ML = kvec(:,index(1));                       
           d=d(:,1);
        else
           kvec=kvec(:,index(1:K));                                 
           d=d(:,1:K);                                              
        end
    end

end