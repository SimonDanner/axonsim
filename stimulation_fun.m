function [ c ] = stimulation_fun( t , dur, type, fun)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    switch type
        case 1
            if t <= dur & t >= 0
                c = 1;
            else
                c = 0;
            end
        case 2
            if t <= dur & t >= 0
                c = 1;
            elseif t > dur & t<=2*dur
                c = -1;
            else
                c = 0;
            end
        case 3
            c = eval(fun);
    end

end

