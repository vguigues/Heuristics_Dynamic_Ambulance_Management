


function [bool]=can_answer(ambulance_type,call_type,samu_type)

bool=1;
if ((ambulance_type>call_type)&&(samu_type=='rj'))
    bool=0;
end


