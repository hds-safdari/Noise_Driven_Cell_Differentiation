function b=BondCkeck(a,lower,upper)
    if a<lower
       b=lower;
    elseif a>upper
       b=upper; 
    else
       b=a;
    end