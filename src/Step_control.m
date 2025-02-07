function [step_new] = Step_control(step_old,xx,step_ctl)
%STEP_CONTROL(XX,IS_CONVERGED,STEP_CTL); Summary of this function goes here
%   Detailed explanation goes here
%%
if ~xx(end).is_Converged
    step_new=step_old/2;
else
    switch step_ctl
        case "Automatic"
            if length(xx)<=5
                step_new=step_old;
            else
                if all([xx(end-4:end).is_Converged])
                    step_new=step_old*2;
                else
                    step_new=step_old;
                end
            end
        case "Fixed"
            step_new=step_old;
    end

end
end

