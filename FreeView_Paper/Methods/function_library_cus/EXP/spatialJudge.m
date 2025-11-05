classdef spatialJudge
    properties
        taskSpace; % all possible location, table with 2 column or N * 2 matrix
    end

    methods
        function obj = spatialJudge(taskSpace)
            if istable(taskSpace)
                taskSpace = table2array(taskSpace);
            end
            s = size(taskSpace);
            if s(1)==2 && ~s(2)==2       % correct the matrix with false size          
                taskSpace = transpose(taskSpace);
            end
            obj.taskSpace = taskSpace; 
        end

        function [judgement, Terr] = judge(obj, report, answer)
            Terr = norm(report - answer); % target error
            Aerr = repmat(report, length(obj.taskSpace), 1)-obj.taskSpace;
            Aerr = sqrt(Aerr(:,1).^2 + Aerr(:,2).^2);
            if abs(min(Aerr)-Terr) < Terr * 1e-5 
                 judgement = 1;
            else
                judgement = 0;
            end
        end
    end
end