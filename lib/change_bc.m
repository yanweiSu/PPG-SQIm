function [phase] = change_bc(old_phase, bc)
    if bc >= 0 
        bc  = bc - 0.5;
    else
        bc = 0.5 + bc;
    end
    
    phase = old_phase;
    if bc >= 0
        idx1 = old_phase > -0.5 + bc & old_phase <= 0.5;
        idx2 = old_phase >= -0.5 &  old_phase <= -0.5 + bc;
        phase(idx1) = old_phase(idx1) - bc;
        phase(idx2) = old_phase(idx2) + 1 -bc;
    
    elseif bc < 0
        idx1 = old_phase <= 0.5 + bc & old_phase > -0.5;
        idx2 = old_phase <= 0.5 &  old_phase > 0.5 + bc;
        phase(idx1) = old_phase(idx1) - bc;
        phase(idx2) = old_phase(idx2) -1 -bc;
    
    end
    
end