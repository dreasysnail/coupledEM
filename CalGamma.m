function [alpha,beta,sigma] = CalGamma(boundnum,dist,response)
    m1 = sum(boundnum.*dist,2);
    m2 = sum(boundnum.^2.*dist,2);
    nsize = size(response,1);
    alpha = (sum(response)*sum(m2)-sum(m1)*response'*m1)/(nsize*sum(m2)-sum(m1)^2);
    beta = (sum(response)-nsize*alpha)/sum(m1);
    sigma =  sqrt(sum(sum((repmat(response,1,size(boundnum,2)) - beta.*boundnum-alpha).^2.*dist))/nsize);  
end