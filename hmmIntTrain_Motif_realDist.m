function [guessT,guessE,alpha,beta,sigma,logliks] = hmmIntTrain_Motif_realDist(seqs,response,guessT,guessE,alpha,beta,sigma,changeT)

%   See also  HMMGENERATE, HMMDECODE, HMMESTIMATE, HMMVITERBI.

%   Reference: Biological Sequence Analysis, Durbin, Eddy, Krogh, and
%   Mitchison, Cambridge University Press, 1998.

%   Copyright 1993-2011 The MathWorks, Inc.
assert(all(guessT>=0)&&all(guessT<=1));

tol = 5e-4;
trtol = tol;
etol = tol;
maxiter = 500;
topK = 40;


verbose = true;

baumwelch = true;

smoothing = false;

avg = true;

numStates = size(guessE,1);

len = (numStates-1)/2;

response = reshape(response,[],1);
% 
% intEnhancer = 20;

% number of rows of e must be same as number of states

[checkE, numEmissions] = size(guessE);
if checkE ~= numStates
    error(message('stats:hmmtrain:InputSizeMismatch'));
end
if (numStates ==0 || numEmissions == 0)
    guessT = [];
    guessE = [];
    return
end


%null distribution (no change T)
% MuMap = [];
% for t = 0.001:0.002:1
%     [pStates,~, ~, ~, ~] = hmmdecode(seqs(1,:),wrap_T(t,len),ones(size(guessE))*0.25);
%     MuMap = [MuMap; sum(1-pStates(1,:)),pStates(:)'];
% end
    
    
    





[numSeqs, seqLength] = size(seqs);


% initialize the counters
if ~changeT
    assert(size(guessT,1)==1&&size(guessT,2)==1);
end

guessT = [guessT/2;guessT/2];

T = zeros(size(guessT));
NT = zeros(size(guessT));
E = zeros(numStates,numEmissions);

pseudoE = E;
pseudoNT = NT;
pseudoT = T;
expBound = zeros(numSeqs,1);
%tempvar = zeros(numSeqs,1);

converged = false;
loglik = 1; % loglik is the log likelihood of all sequences given the TR and E
logliks = zeros(1,maxiter);
guessT_history = [];
for iteration = 1:maxiter
    oldLL = loglik;
    loglik = 0;
    oldGuessE = guessE;
    oldGuessT = guessT;
%    oldAlpha = alpha;
%    oldBeta = beta;
%    oldsigma = sigma;
    muConY = zeros(numSeqs,1);
    %
    distAll = zeros(numSeqs,topK);
    boundAll = zeros(numSeqs,topK);
    %sdConY = zeros(numSeqs,1);
    for count = 1:numSeqs
        %disp(count);

        seq = seqs(count,:);
        %seqLength = length(seq);
        res = response(count);
        
        if baumwelch   % Baum-Welch training
            % get the scaled forward and backward probabilities
            [pStates,logPseq, fs, bs, s] = hmmdecode_changeT_ds(seq,guessT,guessE,changeT);
            expBound(count) = seqLength-sum(pStates(1,:)); 
            pStates2 = Pstate2nd_changeT_ds(numStates,seq,fs, bs,s,guessT,guessE,changeT);
            
            % p(#|x,y)
            %[condVitBound,condVitLogP,boundNum,logP,condVitBound2,condVitLogP2] = hmmviterbi_topK_cond_table(seq,wrap_T(guessT,len),guessE,topK);
            [condVitBound,condVitLogP,boundNum,logP,condVitBound2,condVitLogP2] = hmmviterbi_topK_cond_table_c(seq,wrap_T(guessT,len),guessE,topK);

            % p (#|x)= p(#,x)/p(x)
            pp = exp(logP-logPseq);
            if sum(pp)<0.80
                fprintf('Warning! Viterbi paths contain only %d prob of %dth seq\n',sum(pp),count);
            end
                
%             newK = topK;
%             while sum(pp)<0.95&&newK<topK*4
%                 newK = 2*newK;
%                 fprintf('Warning! Dynamically changing K to %d\n',newK);
%                 [condVitBound,condVitLogP,currentState,logP] = hmmviterbi_topK_cond_table(seq,wrap_T(guessT,len),guessE,newK);
%                 pp = exp(logP-logPseq);
%             end
            %simulate the truth by topK
            %pp = pp./sum(pp);
            %quick and dirty fix
            pp = pp./sum(pp);
            % p(#|x,y)\propto p(#|x)p(y|#)
            pp = normpdf(res-alpha-beta.*boundNum,0,sigma).*pp';
            distAll(count,:) = pp./sum(pp);
            boundAll(count,:) = boundNum;
            
            
            % need to be optimized
            %l,j,zl,zj
            Pmat = hmmCondDecode(numStates,seq,pStates,pStates2);
            % mu = expected bound number E(mu|z_l)=sum(p(z_i|x,z_l))
            % we assume first state is unbound.
            % p(#) = norm (mu(#),sd(#))
           
         
            mu = seqLength-squeeze(sum(Pmat(:,:,:,1),2))';
            % #_i~N(E#_i,var#_i) var#_i=tempvar  becareful when apply to multiple
    % states
            % (1-pStates(1,:))*squeeze(sum(sum(Pmat(:,:,2:end,2:end),4),3))
            % tempvar = (#|x)
            %tempvar(count) = sum((1-pStates(1,:))*squeeze(sum(sum(Pmat(:,:,2:end,2:end),4),3))) - expBound(count)^2;
            
            % hmmjoint: p(z_i=bound,z_j=bound|z_l)
            sd = ones(size(mu))*2;
            % since symmetric only consider i>j 
            
            %m1(i,j) = p(z_j=ubound|z_i=ubound)  
            %jointm p(zi=m,zl=n) id by l i n m    can be optimized by
            %integrate into hmmCondDecode
%             jointm_temp = squeeze(Pmat(:,:,:,1)).*permute(repmat(pStates,[1 1 seqLength]),[2 3 1]);
%             jointm = squeeze(jointm_temp(:,:,1));
%             m1 = jointm./repmat(pStates(1,:)+0.000001,seqLength,1)';
%             
%             for zl = 1:numStates
%                 for l = 1:seqLength
%                     %m2(any,i) = p(z_i=ubound|z_l=zl)
%                     m2 = repmat(squeeze(Pmat(l,:,zl,1)),seqLength,1);
%                     sd(zl,l) = sumjoint(m1,m2,l);
%                 end
%             end
%             % sd(#|x,z_l)
%             sd = sqrt(sd - mu.^2);
         
            
            condLik = zeros(size(mu));
            %mu_test = zeros(size(mu));
            % condLik = sum(p(a+b#)N(y-a-b#,sigma^2)) can be approximate by
            % top k # prob
            

            for ll = 1:seqLength
                for ss = 1:numStates
                    % fast procedure only suit for ds
                    if (ll>1)&&(ss~=2)&&(ss~=len+2)&&(ss~=1)
                        condLik(ss,ll) = condLik(ss-1,ll-1);
                        continue;
                    end
                    %condition on zll = ss
                    %[currentState, logP] = hmmviterbi_topK_cond(seq,wrap_T(guessT,len),guessE,topK,ll,ss);
                    logP = condVitLogP(:,ss,ll);
                    boundNum = condVitBound(:,ss,ll);
                    % p(#|zll=s,x) = p(#,x,zll=ss)/p(x)/p(zll=ss|x)
                    pp = exp(logP-logPseq-log(pStates(ss,ll)));
                    pp = pp./sum(pp);
                    %p (y|zll=ss,x) = sum p(#|zll=ss,x)p(y|#)
                    condLik(ss,ll) = normpdf(res-alpha-beta.*boundNum',0,sigma)*pp;
                    %mu_test(ss,ll) = boundNum*pp;
                end
            end
            condLik(isnan(condLik)) = 0;
                    
                    


            % be careful of small number
            posterior = condLik.*pStates;
            % normalize 
            posterior = posterior./repmat(sum(posterior,1),size(posterior,1),1);    
            
            %p(#|y,x,theta_s,gamma_s)~N(E(#|y,x,params_s),Var(#|x,theta_s))
            %E(#|y,x,params_s)
            muConY(count)=sum(1-posterior(1,:));
            
            
          
            
            %P(z_i|z_l,z_l+1,x,theta)
            %Pmat2 = hmmCond2Decode(numStates,seq,pStates,pStates2);
            % mu2 = expected bound number|Z_l,z_l+1 
            % we assume first state is unbound.
            
            % we don't care about state other than stateSet
            stateSet = [1,2,((1:len-1)*(2*len+2))+2,len*(2*len+1)+1,2+len,((len+1:2*len-1)*(2*len+2))+2,2*len*(2*len+1)+1];
            pStates2 = pStates2(stateSet,:);
            cond2Lik = zeros(size(pStates2));
            %we only use cond1 to build cond2Lik
            
            cond2Lik(2:len+1,:) = condLik(2:len+1,:);
            cond2Lik(len+2,:) = [0,condLik(len+1,1:end-1)];
            cond2Lik(len+3:2*len+2,:) = condLik(len+2:2*len+1,:);
            cond2Lik(2*len+3,:) = [0,condLik(2*len+1,1:end-1)];
           

            %a -> b
            %[a,b] = modDecode(stateSet(ss),2*len+1);
            for ll = 1:seqLength                
                %[currentState, logP] = hmmviterbi_topK_cond2(seq,wrap_T(guessT,len),guessE,topK,ll,1,1);    
                %boundNum = sum(currentState~=1,1);
                
                pp = exp(condVitLogP2(1,:)-logPseq-log(pStates2(1,ll)));
                pp = pp./sum(pp);
                cond2Lik(1,ll) = normpdf(res-alpha-beta.*condVitBound2(1,:),0,sigma)*pp';
            end

            %normpdf(res,mu,var)=normpdf(mu,res,var)
            % posterior p(z_l,z_l+1|y,X,theta,gamma)
            % be careful of small number
%             cond2Lik(isnan(cond2Lik)) = 0;
%             cond2Lik(isinf(cond2Lik)) = 0;
            posterior2 = cond2Lik.*pStates2;       
            % normalize 
            posterior2 = posterior2./repmat(sum(posterior2,1),size(posterior2,1),1);
            %validate: sum(posterior2(1:numstate),1)==posterior(1,:)
            
%             for ss = 1:numStates
%                 posterior(ss,1:(end-1)) = sum(posterior2(1+(ss-1)*numStates:ss*numStates,:),1);
%             end
            


           
            
            
            
            %p(x,y|thetaS, gammaS)=p(x|thetaS)p(y|-)
            %p(y|-)  = sum ( p(y|gamma,#) p(#))            
            %logPres = log(normpdf(res,alpha+beta*muConY(count),sigma));
            logPres = log(normpdf(res-alpha-beta.*boundAll(count,:),0,sigma)*distAll(count,:)');
            loglik = loglik + logPseq + logPres;
            

            
            %Intensity Enhancer

%             reverseMu = (res - alpha)/beta;
%             selector = find_halfspace(MuMap(:,1),reverseMu);
%             tempP = reshape(MuMap(selector,2:end),size(posterior));
%             %plot(posterior(1,:),'r');hold on;plot(tempP(1,:));
%             posterior = posterior + tempP*intEnhancer;
%             % normalize 
%             posterior = posterior./repmat(sum(posterior,1),size(posterior,1),1);     
%             
            
            
            
            
            % update T
            assert(size(posterior2,2)==seqLength);
            if changeT
                for ii = 2:seqLength  
                    NT(ii-1) = NT(ii-1) + sum(posterior2(1,ii-1));
                    T(1,ii-1) = T(1,ii-1) + sum(posterior2(2,ii-1));
                    T(2,ii-1) = T(2,ii-1) + sum(posterior2(2+len,ii-1));    
                end
            else
                 NT = NT + sum(posterior2(1,:));
                 T(1) = T(1) + sum(posterior2(2,:));
                 T(2) = T(2) + sum(posterior2(2+len,:));
            end
            
            
            
            for k = 1:numStates
                for j = 1:numEmissions
                    E(k,j) = E(k,j) + sum(posterior(k,seq == j));
                end
            end
            
            
   
            
            
            
        else  % Viterbi training
            [estimatedStates,logPseq]  = hmmviterbi(seq,guessT,guessE);
            loglik = loglik + logPseq;
            % w = warning('off');
            [iterT, iterE] = hmmestimate(seq,estimatedStates,'pseudoe',pseudoE,'pseudoT',pseudoT);
            %warning(w);
            % deal with any possible NaN values
            iterT(isnan(iterT)) = 0;
            iterE(isnan(iterE)) = 0;
            
            T = T + iterT;
            E = E + iterE;
        end
        

        
    end
    totalEmissions = sum(E,2);

    %scatter(mutrue,muConY,'r');hold on;scatter(mutrue,expBound);
    %scatter(muConY,response,'r');hold on;scatter(expBound,response,'b');
    
%     scatter(expBound,response,'r.');
%     hold on;
%     xlabel('E(#|x,theta)');
%     ylabel('log intensity');
%     lsline;
    %scatter(expBound,muConY)
    
    %scatter(sum(distAll.*boundAll,2),muConY)
    
    % avoid divide by zero warnings
    guessE = E./(repmat(totalEmissions,1,numEmissions)); 
    guessT  = T./(sum(T,1)+NT);
    
    if avg
        guessE = modify_e(guessE);
        guessT = [sum(guessT,1)/2;sum(guessT,1)/2];
    end
    

    if ~changeT
        Etrue =  [0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000;0.271530000000000,0.111186000000000,0.145423000000000,0.471861000000000;0.289388000000000,0.117744000000000,0.137428000000000,0.455440000000000;0.265435000000000,0.114924000000000,0.131731000000000,0.487910000000000;0.125645000000000,0.555972000000000,0.0846770000000000,0.233706000000000;0.303948000000000,0.159204000000000,0.331508000000000,0.205341000000000;0.302736000000000,0.261673000000000,0.225970000000000,0.209621000000000;0.209621000000000,0.225970000000000,0.261673000000000,0.302736000000000;0.205341000000000,0.331508000000000,0.159204000000000,0.303948000000000;0.233706000000000,0.0846770000000000,0.555972000000000,0.125645000000000;0.487910000000000,0.131731000000000,0.114924000000000,0.265435000000000;0.455440000000000,0.137428000000000,0.117744000000000,0.289388000000000;0.471861000000000,0.145423000000000,0.111186000000000,0.271530000000000];
        fprintf('T : %d(%d), E->Etrue : %d, mean(mu):%d \n',guessT(1),guessT(2),norm(guessE - Etrue,inf),mean(muConY));
    end

    if smoothing
        guessT = smooth(guessT,0.3,'loess')';
        if any(guessT<0)
            guessT(guessT<0) = 0;
        end
    end
    guessT_history = [guessT_history;sum(guessT,1)];
    
    %Initiate the point
    %muConY = mutrue;
    DEBUG = false;
    if ~DEBUG
        [alpha,beta,sigma]= CalGamma(boundAll,distAll,response);
%          tempCov = cov(response,muConY);
%          beta = tempCov(1,2)/tempCov(2,2);
%          alpha = mean(response)-beta*mean(muConY);
%          sigma = sqrt(sum(((response-alpha-beta*muConY).^2))/numSeqs);
%         options = optimoptions('fminunc','Display','off','TolX',1e-8,'Algorithm','quasi-newton');
%         x0= [alpha,beta];
%         f = @(x)Qgamma(boundAll,distAll,response,x);
%         x1 = fminunc(f,x0,options);
%         alpha = x1(1);
%         beta = x1(2);
%         sigma = sqrt(f(x1)/numSeqs);    
    end
    
  
    
  
    % clean up any remaining Nans
    guessT(isnan(guessT)) = 0;
    guessE(isnan(guessE)) = 0;
    if iteration == 21
        plot(guessT);
        rl = {'b','m1','m2','m3','m4','m5','m6','c1','c2','c3','c4','c5','c6'};
        cl= {'A','C','G','T'};
        HeatMap(guessE,'RowLabels',rl,'ColumnLabels',cl,'Colormap','autumn','Symmetric',true);

    end
    if verbose
        if iteration == 1
            fprintf('%s\n',getString(message('stats:hmmtrain:RelativeChanges')));
            fprintf('   Iteration       Log Lik    Transition     Emmission  alpha  beta sigma\n');
        else 
            fprintf('  %6d      %6g  %6g  %6g   %6g   %6g   %6g\n', iteration, ...
                (abs(loglik-oldLL)./(1+abs(oldLL))), ...
                norm(guessT - oldGuessT,inf)./numStates, ...
                norm(guessE - oldGuessE,inf)./numEmissions, ...
                alpha,...
                beta,...
                sigma);
        end
    end
    % Durbin et al recommend loglik as the convergence criteria  -- we also
    % use change in TR and E. Use (undocumented) option trtol and
    % etol to set the convergence tolerance for these independently.
    %
    logliks(iteration) = loglik;
    if (abs(loglik-oldLL)/(1+abs(oldLL))) < tol
        if norm(guessT - oldGuessT,inf)/numStates < trtol
            if norm(guessE - oldGuessE,inf)/numEmissions < etol
                if verbose
                    fprintf('%s\n',getString(message('stats:hmmtrain:ConvergedAfterIterations',iteration)))
                end
                converged = true;
                break
            end
        end
    end
    E =  pseudoE;
    T = pseudoT;
    NT = pseudoNT;
end
if ~converged
    warning(message('stats:hmmtrain:NoConvergence', num2str( tol ), maxiter));
end
logliks(logliks ==0) = [];
end



function enew=modify_e(e)
    len = (size(e,1)-1)/2;
    all = (e(2:len+1,:) + e((2*len+1):-1:len+2,4:-1:1))/2;
    enew = [e(1,:);all;all(len:-1:1,4:-1:1)];
end   





