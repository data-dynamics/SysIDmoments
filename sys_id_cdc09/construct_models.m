
function [b group] = construct_models(c,powers,x,vermap,num_sys)

[Dpn,normDpn] = cnormalize(derivative(c,powers,x)); % Dpn = norm 1

distance = abs(vermap*c)'./normDpn;
[val,index] = min(distance);
% COMPUTE NORMAL TO FIRST SUBSPACE FROM THE JACOBIAN OF THE POLYNOMIAL
b(:,1) = cnormalize(derivative(c,powers,x(:,index)));

delta = 0.0001;
% FOR THE REMAINING SUBSPACES
den = 1;
for i = 1:num_sys-1
  % PICK A POINT IN ONE SUBSPACE
  den         = den.*sqrt(sum(conj(conj(b(:,i)')*x).*(conj(b(:,i))'*x),1)); 
  [val,index] = min((distance + delta)./(den + delta));
  b(:,i+1)    = cnormalize(derivative(c,powers,x(:,index)));
end

% Segment the points: assign z to b such that (z'b) is mimimum
distance = abs(conj(b)'*cnormalize(x));
[val,group] = min(distance,[],1); 

%%
% sigma_sq = 1e-5;
% for i=1:num_sys
%     if ~isempty(find(group==i,1))
%         sigma_sq = sigma_sq + var(distance(find(group==i)));
%     end
% end
% sigma_sq = sigma_sq /num_sys;
% 
% membership = exp(-distance.^2/sigma_sq);
% membership = membership./(ones(num_sys,1)*sum(membership,1));

