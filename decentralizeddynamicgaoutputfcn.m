


function [state,options,optchanged] = decentralizeddynamicgaoutputfcn(options,state,flag)
persistent history averagescore 
optchanged = false;
switch flag
    
    case 'init'
        
        history(:,:,1) = state.Population;
        filename = sprintf('Generation_%d_Population',1);
        save(filename,'history','-mat','v7.3')
    
    case 'iter'
        %Update and save history every 50 generations
        if rem(state.Generation,50) == 0
            ss = size(history,3);
            history(:,:,ss+1) = state.Population;
            filename = sprintf('Generation_%d_Population',(ss+1)*10);
            save(filename,'history','-mat','v7.3')
        end

        % Find the best objective function, and stop if it is low.
        ibest = state.Best(end);
        ibest = find(state.Score == ibest,1,'last');
        bestx = state.Population(ibest,:);
        bestf = decentdynoutputfeedbackcostfcn_nested(bestx);
        if bestf <= -99.8
            state.StopFlag = 'y';
            disp('Cost function evaluation is -99.5')
        end
        
        averagescore(state.Generation) = mean(state.Score);

    case 'done'
        exportbest = state.Best;
        exportaverage = averagescore;
        exportscorename = sprintf('AvgandBestScores');
        save(exportscorename,'exportbest','exportaverage','-mat','v7.3')
        
        ss = size(history,3);
        history(:,:,ss+1) = state.Population;
        filename = sprintf('Generation_%d_Population',ss+1);
        save(filename,'history','-mat','v7.3')
end