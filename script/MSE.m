

%% MSE
% a = a/100; a = a*100;
% b = b/100;
% a =load('C:\Users\Administrator\Desktop\wxy\deepLearning\wangxueyuNew\1.txt')*1;
Sum = length(a(:,1));
c = a - b;

MSE(1,1) = mean(c(:,1).^2);
MSE(2,1) = mean(c(:,2).^2);
MSE(3,1) = mean(c(:,3).^2);

save('Bias_voltage_5_2nm.mat','a');