clc;
clear;
close all;
model=SelectModel();      
model.eta=0.1;
CostFunction=@(q) MyCost(q,model);    
MaxIt=1200;   
MaxIt2=80;    
T0=100;        
alpha=0.98;   
x.Position=CreateRandomSolution(model);
[x.Cost, x.Sol]=CostFunction(x.Position);
BestSol=x;
BestCost=zeros(MaxIt,1);
T=T0;
[SBest, Smin]=NBA();
for it=1:MaxIt
    for it2=1:MaxIt2
        xnew.Position=CreateNeighbor(x.Position);
        [xnew.Cost, xnew.Sol]=CostFunction(xnew.Position);        
        if xnew.Cost<=x.Cost
            x=xnew;
        else
            delta=xnew.Cost-x.Cost;
            p=exp(-delta/T);
            if rand<=p
                x=xnew;
            end
            
        end
        if x.Cost<=BestSol.Cost
            BestSol=x;
        end
    end
    BestCost(it)=BestSol.Cost;
    if BestSol.Sol.IsFeasible
        FLAG=' *';
    else
        FLAG='';
    end
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it)) FLAG]);
    
    T=alpha*T;
    figure(1);
    PlotSolution(BestSol.Sol,model);
    pause(0.01);
    
end
figure;
plot(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

