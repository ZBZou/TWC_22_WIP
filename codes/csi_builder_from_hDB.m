clc
clear

% load CSI_v20ms_Ns5_raw.mat
load CSI_v20ms_Ns5_raw_new.mat

%% For Random Vehicles at Random Locations
% for i=1:length(hDB{1})
%     Z_True_Complex(:,i)=hDB{randi([1,10])}(:,i);    
% end

%% For a single vehcicle taking a U-turn repeatedly

clc
a = 0;
c=1;
for b = 1:1000
    if c==1
        a=a+1;
    else
        a=a-1;
    end
    if a==1 %rem(a/10,2)==1
        c=1;
    elseif a==10 %rem(a/10,2)==0
        c=0;
    end
%     disp(a)
    Z_True_Complex(:,b)=hDB{a}(:,b);    
end
        
%% Re-arrange Subcarriers
h = [Z_True_Complex(32:end,:);Z_True_Complex(1:31,:)];
h = h(7:59,:);                               %subcarrier indeces -26:26
Z_True_Complex=h;
