
% Min energy multivariate test
% p-value < 0.05  (*). 
% p-value < 0.01  (**) 
% p-value < 0.001 (***)
% p returns the p value

function p = multivartest(xmang,dims)

Ang  = length(xmang);
uS   = size(xmang{1},4);
T    = size(xmang{1},2);

ndim = 1; % for multivariate tests, one p-value across all dims

Period = 1:T;
lp = length(Period);

p = zeros(Ang,uS-1,ndim,lp);

for d=1:ndim
    for t=1:lp
        ti = Period(t);
        for a = 1:Ang
            for u = 1:uS              
                dist{a,u,t} = squeeze((xmang{a}(dims,ti,:,u)));
                
                [~,nanidxj] = find(isnan(dist{a,u,t}));
                dist{a,u,t}(:,nanidxj) = [];
                
                if u>1
                    [p(a,u-1,d,t), ~] = ...
                        minentest(dist{a,1,t}',dist{a,u,t}','sr',100);
                end
            end
        end
    end
end

end
