function [Frame_Name] = tc(t2, t1)
    time_distance = etime(t2, t1);

    if time_distance <= 60
        Frame_Name = [num2str(round(time_distance, 2)) ' s']
        time_distance = time_distance / 60;
    end

    if (time_distance > 60) && (time_distance < 60 * 60)
        Frame_Name = [num2str(round(time_distance / 60, 2)) ' min']
    end

    if time_distance >= 60 * 60
        Frame_Name = [num2str(round(time_distance / 60/60, 2)) ' h']
    end

end
