function [x,fval,exitflag,output,population,score] = GA_code(nvars,lb,ub,PopulationSize_Data,MaxGenerations_Data,MaxStallGenerations_Data,InitialPopulationMatrix_Data)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimoptions('ga');
%% Modify options setting
options = optimoptions(options,'PopulationSize', PopulationSize_Data);
options = optimoptions(options,'MaxGenerations', MaxGenerations_Data);
options = optimoptions(options,'MaxStallGenerations', MaxStallGenerations_Data);
options = optimoptions(options,'InitialPopulationMatrix', InitialPopulationMatrix_Data);
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', { @gaplotbestf });

options = optimoptions(options,'UseParallel', false); % true


[x,fval,exitflag,output,population,score] = ...
ga(@objectivefunction,nvars,[],[],[],[],lb,ub,[],[],options);
