function [param,lens] = arclenparam(curve)
% parameterize curve by arclength

curvediff = diff(curve');
lens = sqrt(curvediff(:,1).^2 + curvediff(:,2).^2);
param = [0 cumsum(lens')];

% param = zeros(1,size(curve,2));
% for c = 2:size(curve,2)
%     len = norm(curve(:,c)-curve(:,c-1));
%     param(c) = param(c-1)+len;
% end

lens = [lens' norm(curve(:,1)-curve(:,end))];

end