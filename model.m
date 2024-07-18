function [ t,Y,N_nodes,b_thr] = model(mdl, dur,data,stim_dur,fun_type,custom_fun,fiberD,frq,end_on_ap,I_intra,N_intra )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if exist('I_intra')==1
        switch mdl
            case 1
                [ t,Y,N_nodes,b_thr]=mcintyre(dur,data,stim_dur,fun_type,custom_fun,fiberD,frq,end_on_ap,I_intra,N_intra);
            case 2
                [ t,Y,N_nodes,b_thr]=crrss(dur,data,stim_dur,fun_type,custom_fun,fiberD,frq,end_on_ap,I_intra,N_intra);
            case 3
                [ t,Y,N_nodes,b_thr]=mcneal(dur,data,stim_dur,fun_type,custom_fun,fiberD,frq,end_on_ap,I_intra,N_intra);
        end
    end
    if exist('I_intra')==0
       switch mdl
            case 1
                [ t,Y,N_nodes,b_thr]=mcintyre(dur,data,stim_dur,fun_type,custom_fun,fiberD,frq,end_on_ap);
            case 2
                [ t,Y,N_nodes,b_thr]=crrss(dur,data,stim_dur,fun_type,custom_fun,fiberD,frq,end_on_ap);
            case 3
                [ t,Y,N_nodes,b_thr]=mcneal(dur,data,stim_dur,fun_type,custom_fun,fiberD,frq,end_on_ap);

        end
    end
                
end

