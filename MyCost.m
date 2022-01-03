function [z, sol]=MyCost(q,model)

    sol=ParseSolution(q,model);
    
    eta=model.eta;
    
    z1=eta*sol.TotalD+(1-eta)*sol.MaxD;
    
    beta=5;
    
    z=z1*(1+beta*sol.MeanCV);

end