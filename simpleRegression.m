function [alpha,beta,sigma] = simpleRegression(y,x)
    assert(length(y)==length(x));
    tempCov = cov(y,x);
    beta = tempCov(1,2)/tempCov(2,2);
    alpha = mean(y)-beta*mean(x);
    sigma = sqrt(sum(((y-alpha-beta*x).^2))/length(x));
    disp([alpha,beta,sigma]);
end
