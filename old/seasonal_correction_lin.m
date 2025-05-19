function [lin_out,ml_out] = seasonal_correction_lin(tin,pin)
%
% Add linear seasonal fit to pre-existing save structure
% 

for l=1:length(pin)
    t=tin{l};
    p=pin{l};
    if isempty(p)
        continue
    else
        dim=size(t);
        if dim(1)<dim(2) % transpose as needed
            t=t';
            p=p';
        end
    end
    lamda_list=linspace(1/max(t-t(1)),180/max(t-t(1)),1000);
    tinv=(t-t(1))/max(t-t(1)); % better time basis

    %----- exponential only
    for jj=1:1000
        %construct matrix for inversion
        Gl=[tinv,ones(size(tinv))];
        gexp=exp(-tinv/lamda_list(jj)); gexp(gexp<10^-7)=0;
        Gl=[Gl,gexp];
        % rcond(Gl'*Gl)
        ml{l}(:,jj)=inv(Gl'*Gl)*Gl'*p;
        ple_fit{l}(:,jj)=Gl*ml{l}(:,jj);
        stds{l}(jj)=std(p-ple_fit{l}(:,jj));
    end
    [~,imin]=min(stds{l});
    lin_out{l}=p-ple_fit{l}(:,imin);
    ml_out{l}=[ml{l}(:,imin);lamda_list(imin)];
    keyboard % intervene as necessary
end