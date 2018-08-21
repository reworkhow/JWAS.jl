

#deal with the problem that some fixed effects are not estimated from ssBayes.
#Because some fixed effects only appeared for observations without phenotypes.

function get_correlation(file) #ID, fiexed, phenotype, ebv
    d=readdlm(file);
    fixed=int(d[:,2]);
    X,effects=QTL.mkmat_incidence_factor(fixed);
    X=[ones(size(X,1)) X];
    y=[d[:,3] d[:,4]]
    betaHat=pinv(X'X)*X'y;
    res=y-X*betaHat
    varcov=res'res/(size(res,1)-length(effects))
    cor=varcov[1,2]/sqrt(varcov[1,1]*varcov[2,2])
    return cor,betaHat
end
