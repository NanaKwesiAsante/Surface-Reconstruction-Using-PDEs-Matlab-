%%COS701 ADDITIONAL PROJECT

% Vivian Montiforte
% Kwesi Acheampong

% COS 701
% Fall 2019
% Research Project

clear variables;

%% Loading/Storing Data
t = readtable('fig14b2.txt');
 g=t(:,[1,4]);
 x=t{:,1};
 y=t{:,2};
 z=t{:,3};
 scatter3(x,y,z,'b.')