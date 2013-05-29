function y = opA_spg(x,Af, At, mode)
    if mode==1
       y = Af(x); 
    else
       y = At(x);
    end
end
