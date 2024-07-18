function x = stim( t,length )
   
    if t < length
        x = sin(t*pi/length);
    else
        x=0;
    end
    
end

