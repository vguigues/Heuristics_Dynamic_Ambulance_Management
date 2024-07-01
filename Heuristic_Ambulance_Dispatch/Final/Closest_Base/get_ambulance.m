
function [index_amb]=get_ambulance(ambulances,indexes,current_time)

index_amb=-1;
initType=0;
for i=1:length(indexes)
    if (ambulances(indexes(i)).arrival_time_at_f_last_trip<=current_time)
        if (ambulances(indexes(i)).type>initType)
                initType=ambulances(indexes(i)).type;
                index_amb=indexes(i);
        end
    end         
end