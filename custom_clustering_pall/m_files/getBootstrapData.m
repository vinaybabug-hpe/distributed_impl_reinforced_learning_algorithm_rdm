function [ bootstrapData ] = getBootstrapData( X, BootIdxs, B )
%getBootstrapData - Generate bth bootstrap data of X from boot index
%
% Syntax:  bootstrapData = getBootstrapData( X, BootIdxs, B )
%
%   INPUTS ARGUMENTS
%         X  : [n,p] data
%         BootIdxs : {n,B} cell array of bootstrap data index
%         B : Bth bootstrap

%
%   OUTPUT ARGUMENTS
%   	E   : Matrix containiing data of Bth bootstrap
%
%
%
% Author:         Vinay B Gavirangaswamy
% Affiliation:    Western Michigan University, CS Dept.
% email:          vinay.b.gavirangaswamy@wmich.edu
% Website:        
% Created:        August, 2015
% Revised by/on:  person, date

%------------- BEGIN CODE --------------

bootstrapData=X(BootIdxs(:,B),:);
end

