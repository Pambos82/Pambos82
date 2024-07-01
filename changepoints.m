function changePoints = changepoints(sub, mask1, maskMinus1, sig,funcHandle)
    % Initialize loop variable
    m = 1;
    % Initialize pvalues array
   

    % Initialize the change points array
    changePoints = [];
    %mask1(:)=1;%if you wish to use one martingale also remove +funcHandle(pvalues2,15, 10, sig);
    while m <= length(sub)
        pvalues=[];
        pvalues2=[];
        for i = m:length(sub)
            if mask1(i) == 1
                pvalues(i-m+1) = (sum(sub(mask1(m:(i-1))) < sub(i)) + rand * sum(sub(mask1(m:(i))) == sub(i))) / sum(mask1(m:i));
            sub_(i-m+1)=median(sub(mask1(m:(i))));
             pvalues2(i-m+1) = (sum(sub_(mask1(m:(i-1))) < sub_(i)) + rand * sum(sub_(mask1(m:(i))) == sub_(i))) / sum(mask1(m:i));

            else
                pvalues(i-m+1) = (sum(sub(maskMinus1(m:(i-1))) < sub(i)) + rand * sum(sub(maskMinus1(m:(i))) == sub(i))) / sum(maskMinus1(m:i));
            sub_(i-m+1)=median(sub(maskMinus1(m:(i))));
             pvalues2(i-m+1) = (sum(sub_(maskMinus1(m:(i-1))) < sub_(i)) + rand * sum(sub_(maskMinus1(m:(i))) == sub_(i))) / sum(maskMinus1(m:i));

            end
        end
          
if any(isnan(pvalues))
    numNaNs = sum(isnan(pvalues));
    pvalues(isnan(pvalues)) = rand(1, numNaNs); % Replace NaNs with random values from a uniform distribution [0, 1]
    pvalues2(isnan(pvalues2)) = rand(1, numNaNs); % Replace NaNs with random values from a uniform distribution [0, 1]

end
        results1 = funcHandle(pvalues,15, 10, sig)+funcHandle(pvalues2,15, 10, sig);%adding two martingales
        
        % Update loop variable m based on the latest instruction
        if (max(results1)>log(1/sig))
            index=find(results1>log(1/sig), 1, 'first');
        
            m = m + index + 1;
        
        else
                       break;
 
        end
        % Store the updated value of m-1 as a change point
        changePoints = [changePoints, m-1];
        
        % Example condition to exit loop if m exceeds or matches the length of sub
        if m >= length(sub)
            break;
        end
    end
end
