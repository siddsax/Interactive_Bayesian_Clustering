function [ req_probab ] = P_h_hs( q,P_h_given_x, N) 
% q is KXN
% P_h_givn_x is a NXK
req_probab = q*P_h_given_x; %KXK   
req_probab = req_probab/N;
end

