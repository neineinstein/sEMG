function [ Z, ZX ] = zoneBuilder(a,l, data, samples)
%% input a data vector and a corresponding sample vector, 
%return segments centered at a, of length l from data and samples

% number of zones to isolate
    m = length(a);
   
% length of data and sample vectors
    n = length(data);
    
% create memory space for variables    
    Z = zeros(m,l+1);
    ZX = zeros(m,l+1);


% perform the segmentation
    for j = 1:1:m
        Z(j,:) = data(floor(a(j)-l/2):floor(a(j)+l/2));
        ZX(j,:) = samples(floor(a(j)-l/2):floor(a(j)+l/2));
    end




%% stuff i never implemnted
    
%     lowStop = (floor(a(1)-l/2));
%     highStop = (floor(a(m)+l/2));
%     
%     
%     
%     if (lowStop < 0 && highStop > n)
%         caseLogic = 1;
%         % create memory space for variables    
%         Z = zeros(m,l-1);
%         ZX = zeros(m,l-1);
%     else
%         caseLogic = 2;
%         % if just low 
%         if (lowStop < 0)
%             % create memory space for variables    
%             Z = zeros(m,l);
%             ZX = zeros(m,l);      
%         else
%             caseLogic = 3;
%             % if just high
%             % create memory space for variables    
%             Z = zeros(m,l);
%             ZX = zeros(m,l);
%         end
%     end
%     
%         switch caseLogic
%             case 1
%                   % if both out of range
%                     for j = 2:1:(m-1)
%                         Z(j,:) = data(floor(a(j)-l/2):floor(a(j)+l/2));
%                         ZX(j,:) = samples(floor(a(j)-l/2):floor(a(j)+l/2));
%                     end
%                case 2
%                         for j = 2:1:(m)
%                             Z(j,:) = data(floor(a(j)-l/2):floor(a(j)+l/2));
%                             ZX(j,:) = samples(floor(a(j)-l/2):floor(a(j)+l/2));
%                         end
%                case 3
%                         for j = 1:1:(m-1)
%                                 Z(j,:) = data(floor(a(j)-l/2):floor(a(j)+l/2));
%                                 ZX(j,:) = samples(floor(a(j)-l/2):floor(a(j)+l/2));
%                         end
%         end
%         
%         
%     end
    
%     % perform the segmentation
%     for j = 1:1:m
%             Z(j,:) = data(floor(a(j)-l/2):floor(a(j)+l/2));
%             ZX(j,:) = samples(floor(a(j)-l/2):floor(a(j)+l/2));
%     end
    
    



