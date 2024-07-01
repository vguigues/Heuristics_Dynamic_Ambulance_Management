

function [whichk,maxIndex,travelTimes,index_ambs]=foundambulance(index_ambs,calls,maxIndex,i,ambulances,index_amb,isAllocated,nb_calls,nb_ambulances,mintime,penalization)

k=maxIndex+1;
cont=1;
while (cont)
    if (k>nb_calls)
        cont=0;
    else
        if (calls(k).time>ambulances(index_amb).arrival_time_at_f_last_trip)
            cont=0;
        else
            if (isAllocated(k)==0)
                for j=1:nb_ambulances
                    if (ambulances(j).type<=calls(k).priority)
                        if (ambulances(j).arrival_time_at_b_last_trip<=calls(k).time)
                            travelTimes(k,j)=travel_time(ambulances(j).base.lat,ambulances(j).base.long,calls(k).lat,calls(k).long,ambulances(j).speed);
                        elseif (ambulances(j).arrival_time_at_f_last_trip<=calls(k).time)
                            [curr_lat,curr_long]=position_between_origin_destination(ambulances(j).free_destination.lat,ambulances(j).free_destination.long,ambulances(j).base.lat,ambulances(j).base.long,ambulances(j).arrival_time_at_f_last_trip,calls(k).time,ambulances(j).speed);
                            travelTimes(k,j)=travel_time(curr_lat,curr_long,calls(k).lat,calls(k).long,ambulances(j).speed);
                        else
                            time_to_hospital=ambulances(j).arrival_time_at_f_last_trip-calls(k).time;
                            time_from_hosp_to_call=travel_time(ambulances(j).free_destination.lat,ambulances(j).free_destination.long,calls(k).lat,calls(k).long,ambulances(j).speed);
                            travelTimes(k,j)=time_to_hospital+time_from_hosp_to_call;
                        end
                    end
                end
                mintime(k)=inf;
                for j=1:nb_ambulances
                    %First check if ambulance j can take calls of priority calls(i).priority)
                    if (ambulances(j).type<=calls(k).priority)
                        if (travelTimes(k,j)<mintime(k))
                            mintime(k)=travelTimes(k,j);
                            index_ambs{1,k}=[j];
                        elseif (travelTimes(k,j)==mintime(k))
                            index_ambs{1,k}=[index_ambs{1,k};j];
                        end
                    end
                end
            end
            k=k+1;
        end
    end
end

maxIndex=k-1;
maxmin=-1;
whichk=0;
t1=penalization(calls(i).priority)*mintime(i);
k=i+1;
if (i<nb_calls)
    if (calls(k).time<=ambulances(index_amb).arrival_time_at_f_last_trip)
        contd=1;
    else
        contd=0;
    end
    while contd
        if (isAllocated(k)==0)
            t2=penalization(calls(k).priority)*mintime(k);
            if ((mintime(k)==travelTimes(k,index_amb))&&(t2>maxmin))
                whichk=k;
                maxmin=t2;
            end
        end
        if (k+1>nb_calls)
            contd=0;
        elseif (calls(k+1).time>ambulances(index_amb).arrival_time_at_f_last_trip)
            contd=0;
        else
            k=k+1;
        end
    end
end

if (t2<=t1)
    whichk=0;
end