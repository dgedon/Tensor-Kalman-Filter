
clear 
close all
clc 


% Select System Identification (1: identified, 0: A=Id)
selectSysID= 0;
iter= 5e4;


% simple system identification
if selectSysID
    tempSelID= 'Identified';
else
    tempSelID= 'AIdentity';
end

clear
for cLenslets= 2:40

    temp= ['00_SystemID/Identified/AOSystemID_c',num2str(cLenslets),...
        '_iter',num2str(5e4),'.mat'];
    load(temp);
    temp= ['00_SystemID/Identified/AOSystemID_c',num2str(cLenslets),...
        '_iter',num2str(5e4),'.mat'];
    save('test.mat','data')
    a=1;
end