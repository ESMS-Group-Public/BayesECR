function[outputFlag] = METROPOLIS(proposedLikely, LOGLIKELY, fact)
% METROPOLIS accepts or rejects a drawn candidate

  someNum = fact*(proposedLikely - LOGLIKELY);

  if(someNum > 0)
    outputFlag = true;
  else 
    if(someNum > log(rand()))
      outputFlag = true;
    else
      outputFlag = false; 
    end
  end
end