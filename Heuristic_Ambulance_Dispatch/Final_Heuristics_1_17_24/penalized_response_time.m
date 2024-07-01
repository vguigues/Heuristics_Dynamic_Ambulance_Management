

function [penalized_response]=penalized_response_time(penalties,real_time,penalty_matrix,samu_type,ambulance_type,call_type,quality_coeff)

if (samu_type=='rj')
    penalized_response=penalties(call_type)*real_time;
elseif (samu_type=='us')
    penalized_response=penalties(call_type)*real_time+quality_coeff*penalty_matrix(ambulance_type,call_type);
end