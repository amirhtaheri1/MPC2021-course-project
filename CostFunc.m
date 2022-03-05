function [J,PBpb] = CostFunc(X,U,e,data,params)
% params = [Load_up(i),psunny_up(i),tariff(i),savecontrol,saveState,PBpb]
% U = [Pgrid,PBli,Pload,Pgen]
% U
% e
sigma1=1.5;   % Weight factor Pgrid
alpha1=5e-2;  % Weight factor Pgrid
lambda1=5e-1; % Weight factor Pgrid

sigma2=5e-3;    % Weight factor SOCpb
% alpha2=;      % Weight factor SOCpb
% lambda2=5e-1; % Weight factor SOCpb

sigma3=6e-5;  % Weight factor SOCli
alpha2=2e-2;  % Weight factor SOCli
lambda2=2e-2; % Weight factor SOCli

sigma4=8e-2;  % Weight factor Ploads
alpha3=1e-2;  % Weight factor Ploads
lambda3=8e-2; % Weight factor Ploads


global PBpb
PBpb=U(:,1)+U(:,4)-U(:,3)-U(:,2);
% U(:,1)
% % testing
% if sum(U(:,1)>1)
%     disp('')
% end
sell_index = U(:,1)<0;
buy_index = U(:,1)>=0;
% J1=1000*sigma1*sum( U(:,1).*U(:,1)*params{3} ) + sum( alpha1*(U(:,1).*U(:,1))...
%     + lambda1*( U(:,1)-repmat(params{4}(1),[length(U),1]) ).*(U(:,1)-repmat(params{4}(1),[length(U),1])));
J1=1*sigma1*sum( U(buy_index,1).*U(buy_index,1)*params{3} ) +...
   sigma1*sum( U(sell_index,1).*U(sell_index,1)*params{3} ) +...
   + sum( alpha1*(U(:,1).*U(:,1))...
   + lambda1*( U(:,1)-repmat(params{4}(1),[length(U),1]) ).*(U(:,1)-repmat(params{4}(1),[length(U),1])));
J2=sigma2*sum(sum(X.*X));
J3=sigma3*sum(sum(X.*X)) + sum( alpha2* (U(:,2).*U(:,2)) + lambda2*(U(:,2)...
    -repmat(params{4}(2),[length(U),1])).*(U(:,2)-repmat(params{4}(2),[length(U),1])));
J4=sum( PBpb.*PBpb*(1661.9*30)/(3000*0.1*48*.94) );
J5=sum( U(:,2)*(1060.2*30*.92)/(2*1500*367) + 10^(-6)*U(:,2).*U(:,2)) ;
J6=1e5*sigma4*sum(sum( (U(:,3) - params{1}) .* (U(:,3) - params{1}) )) +...
    sum( alpha3*U(:,3).*U(:,3) + lambda3*(U(:,3)-repmat(params{4}(3),[length(U),1])).*(U(:,3)-repmat(params{4}(3),[length(U),1])));
% PBpb
% [J1 J2 J3 J4 J5 J6]
PBpb=PBpb(1);
J=J1+J2+J3+J4+J5+J6;
end
