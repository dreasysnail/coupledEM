function [nBoundtable,logPtable,AllStateBound,AlllogP,nBoundtable2,logPtable2] = hmmviterbi_topK_cond_table(seq,tr,e,k) 

    % tr must be square

    numStates = size(tr,1);
    checkTr = size(tr,2);
    if checkTr ~= numStates
        error(message('stats:hmmviterbi:BadTransitions'));
    end

    % number of rows of e must be same as number of states

    checkE = size(e,1);
    if checkE ~= numStates
        error(message('stats:hmmviterbi:InputSizeMismatch'));
    end

    numEmissions = size(e,2);

    % work in log space to avoid numerical issues
    L = length(seq);
    if any(seq(:)<1) || any(seq(:)~=round(seq(:))) || any(seq(:)>numEmissions)
         error(message('stats:hmmviterbi:BadSequence', numEmissions));
    end

    if L == 0
        return
    end
    logTR = log(tr);
    logE = log(e);

 
    % allocate space  
    pTR = zeros(2,k,numStates,L);
    % assumption is that model is in state 1 at step 0
    % forward vector
    Vf = zeros(numStates,k,L);
    v = -Inf(numStates,k);
    v(1,1) = 0;
    vOld = v;

    % loop through the model
    for count = 1:L
        for state = 1:numStates
            % for each state we calculate
            % v(state) = e(state,seq(count))* max_k(vOld(:)*tr(k,state));

            % use a loop to avoid lots of calls to max
            val = vOld + repmat(logTR(:,state),1,k);
            [bestPTRK,bestValK] = topK_in_matrix(val,k);
            % save the best transition information for later backtracking
            % based on speed consideration, change the order
            pTR(:,:,state,count) = bestPTRK';
            % update v
            v(state,:) = logE(state,seq(count)) + bestValK;
        end
        Vf(:,:,count) = v;
        vOld = v;
    end
    
    %decide which of the final states is post probable
    [finalState,AlllogP] = topK_in_matrix(v,k);

    % Now back trace through the model
    currentState = zeros(L,k,2);
    currentState(L,:,:) = finalState;
    for count = L-1:-1:1
        for kk = 1:k
            currentState(count,kk,:) = pTR(:,currentState(count+1,kk,2),currentState(count+1,kk,1),count+1);
            if any(currentState(count)) == 0
                error(message('stats:hmmviterbi:ZeroTransitionProbability', currentState( count + 1 )));
            end
        end
    end
    AllStateBound = sum(currentState(:,:,1)>1,1);
     
    %%backward
    % allocate space
    currentState = zeros(L,k,2);
    BpTR = zeros(2,k,numStates,L);
    % assumption is that model is in state 1 at step 0
    % forward vector
    Vb = zeros(numStates,k,L);
    v = -Inf(numStates,k);
    v(:,1) = logE(:,seq(L));
    Vb(:,:,L) = v;
    vOld = v;

    % loop through the model
    for count = L-1:-1:1
        for state = 1:numStates
            % for each state we calculate
            % v(state) = e(state,seq(count))* max_k(vOld(:)*tr(k,state));

            % use a loop to avoid lots of calls to max
            val = vOld + repmat(logTR(state,:)',1,k);
            [bestPTRK,bestValK] = topK_in_matrix(val,k);
            % save the best transition information for later backtracking
            BpTR(:,:,state,count) = bestPTRK';
            % update v
            v(state,:) = logE(state,seq(count)) + bestValK;
        end
        Vb(:,:,count) = v;
        vOld = v;
    end
    % enforce start state is 1
    % v = v + repmat(logTR(1,:)',1,k);
    
    % decide which of the final states is post probable
%     [finalState,logP] = topK_in_matrix(v,k);
% 
%     % Now trace forward through the model
%     currentState(1,:,:) = finalState;
%     for count = 2:L
%         for kk = 1:k
%             currentState(count,kk,:) = squeeze(BpTR(currentState(count-1,kk,1),count-1,currentState(count-1,kk,2),:));
%             if any(currentState(count)) == 0
%                 error(message('stats:hmmviterbi:ZeroTransitionProbability', currentState( count + 1 )));
%             end
%         end
%     end
%     currentState = squeeze(currentState(:,:,1));
    
    % result    zl,l,(#,logp) order is for avoid squeeze
    nBoundtable = zeros(k,numStates,L);
    logPtable = zeros(k,numStates,L);
    %conditional table
    for count = 1:L
        for state = 1:numStates
            %forward
            currentState = zeros(L,k,2);
            currentState(count,:,:) = [ones(1,k)*state;1:k]';
            for fcount = count-1:-1:1
                for kk = 1:k
                    currentState(fcount,kk,:) = pTR(:,currentState(fcount+1,kk,2),currentState(fcount+1,kk,1),fcount+1);
                end
            end
            FcurrentState = currentState(:,:,1);
            FlogP = Vf(state,:,count);
            
            %backward
            currentState = zeros(L,k,2);
            currentState(count,:,:) = [ones(1,k)*state;1:k]';
            for bcount = count+1:L
                for kk = 1:k
                    currentState(bcount,kk,:) = BpTR(:,currentState(bcount-1,kk,2),currentState(bcount-1,kk,1),bcount-1);
                end
            end
            currentState(count,:,:) = zeros(k,2);
            BcurrentState = currentState(:,:,1);
            BlogP = Vb(state,:,count);
            
            %
            temp = repmat(FlogP',1,k) + repmat(BlogP,k,1);
            [finalState,logP] = topK_in_matrix(temp,k);
            catState = zeros(L,k);
            for kk = 1:k
                catState(:,kk) = FcurrentState(:,finalState(kk,1)) + BcurrentState(:,finalState(kk,2));
            end
            % # of bind
            nBoundtable(:,state,count) = sum(catState~=1,1);
            % logP
            logPtable(:,state,count) = logP-logE(state,seq(count));
            %test
%             [currentStateTest, logPTest] = hmmviterbi_topK_cond(seq,tr,e,k,count,state);
% 
%             if(~all(isinf(logP)))
%                 assert(all(all(currentStateTest==catState)));
%                 assert(all(logP-logE(state,seq(count))-logPTest<1e-4));
%             end
        end
    end 
% condition on 2 table, for motif case only, only for zl1=1 and zl2=1 
    nBoundtable2 = zeros(L,k);
    logPtable2 = zeros(L,k);
    nBoundtable2(1,:) = nBoundtable(:,1,1);
    logPtable2(1,:) = logPtable(:,1,1);
    for count = 2:L
        %forward
        currentState = zeros(L,k,2);
        currentState(count-1,:,:) = [ones(1,k);1:k]';
        for fcount = count-2:-1:1
            for kk = 1:k
                currentState(fcount,kk,:) = pTR(:,currentState(fcount+1,kk,2),currentState(fcount+1,kk,1),fcount+1);
            end
        end
        FcurrentState = currentState(:,:,1);
        FlogP = Vf(1,:,count-1);

        %backward
        currentState = zeros(L,k,2);
        currentState(count,:,:) = [ones(1,k);1:k]';
        for bcount = count+1:L
            for kk = 1:k
                currentState(bcount,kk,:) = BpTR(:,currentState(bcount-1,kk,2),currentState(bcount-1,kk,1),bcount-1);
            end
        end
        currentState(count,:,:) = zeros(k,2);
        BcurrentState = currentState(:,:,1);
        BlogP = Vb(1,:,count);

        %
        temp = repmat(FlogP',1,k) + repmat(BlogP,k,1);
        [finalState,logP] = topK_in_matrix(temp,k);
        catState = zeros(L,k);
        for kk = 1:k
            catState(:,kk) = FcurrentState(:,finalState(kk,1)) + BcurrentState(:,finalState(kk,2));
            catState(count,:) = ones(1,k);
        end
        % # of bind
        nBoundtable2(count,:) = sum(catState~=1,1);
        % logP
        logPtable2(count,:) = logP-logE(1,seq(count-1))-logE(1,seq(count));
        %test
%             [currentStateTest, logPTest] = hmmviterbi_topK_cond(seq,tr,e,k,count,state);
% 
%             if(~all(isinf(logP)))
%                 assert(all(all(currentStateTest==catState)));
%                 assert(all(logP-logE(state,seq(count))-logPTest<1e-4));
%             end
    end 
    
end

