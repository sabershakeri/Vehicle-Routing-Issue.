function [ bestX, fMin ] = NBA( FitFunc, M, pop, dim, G, gamma, alpha, ...
    r0Max, r0Min, AMax, AMin, freqDMax, freqDMin, probMax, probMin, ...
    CMax, CMin, thetaMax, thetaMin, wMax, wMin )
if nargin < 1
    % 1)The parameters in the basic Bat Algorithm (BA)
    FitFunc = @Sphere;
    M = 1000;   
    pop = 30;   
    dim = 20;   
    gamma = 0.9;
    alpha = 0.99;
    r0Max = 1;
    r0Min = 0;
    AMax = 2;
    AMin = 1;
    freqDMax = 1.5;
    freqDMin = 0;
    
    % 2)The additional parameters in Novel Bat Algorithm (NBA)
    G = 10;   
    probMax = 0.9;
    probMin = 0.6;
    thetaMax = 1;
    thetaMin = 0.5;
    wMax = 0.9;
    wMin = 0.5;
    CMax = 0.9;
    CMin = 0.1;
end
lb= -100 * ones( 1,dim );   % Lower bounds
ub= 100 * ones( 1,dim );    % Upper bounds
vLb = 0.6 * lb;
vUb = 0.6 * ub;

r = rand( pop, 1 ) .* 0.2 + 0;
r0 = rand( pop, 1 ) .* ( r0Max - r0Min ) + r0Min;
A = rand( pop, 1 ) .* ( AMax - AMin ) + AMin;
for i = 1 : pop
    x( i, : ) = lb + (ub - lb) .* rand( 1, dim ); 
    v( i, : ) = rand( 1, dim ); 
    fit( i ) = FitFunc( x( i, : ) ); 
end
pFit = fit; 
pX = x;    
[ fMin, bestIndex ] = min( fit ); 
bestX = x( bestIndex, : );   
bestIter = 1;
 for iteration = 1 : M
    C = rand( pop, 1 ) .* ( CMax - CMin ) + CMin;
    prob = rand( pop, 1 ) .* ( probMax - probMin ) + probMin;
    theta=( thetaMax - thetaMin ) * ( M - iteration )/(1.0 * M) + thetaMin;
    
    freqD = rand( pop, dim ) .* ( freqDMax - freqDMin ) + freqDMin;  
    w = (wMax - wMin) * ( M - iteration )/(1.0 * M) + wMin;
    meanP = mean( pX );
    meanA = mean( A );
    
    for i = 1 : pop
        if rand < prob
            if rand < 0.5
                x( i, : ) = bestX + theta * abs( meanP - pX(i, :) ) *...
                    log( 1.0/rand );
            else
                x( i, : ) = bestX - theta * abs( meanP - pX(i, :) ) *...
                    log( 1.0/rand );
            end
        else
            freqD( i, :) = freqD(i, :) .* ( 340 + v( i, : ) )./( 340 + ...
                v( bestIndex, : ) + realmin );
            v( i, : ) = w .* v( i, : ) + ( bestX - pX(i, :) ) .* ...
                freqD(i,:) .* ( 1 + C(i) .* ( bestX - pX(i, :) ) ./...
                ( abs( bestX - pX(i, :) ) + realmin ) );
            
            v( i, : ) = Bounds( v( i, : ), vLb, vUb );  
            x( i, : ) = x( i, : ) + v( i, : );
        end    
    
        if rand > r( i )
            randnValueA = randn( 1,dim ).* ( abs( A(i) - meanA )+ realmin);
            x( i, : ) = bestX .* ( 1 + randnValueA );
        end   
        
        x( i, : ) = Bounds( x( i, : ), lb, ub );  
        fit( i ) = FitFunc( x( i, : ) );
    end
    for i = 1 : pop 
        if fit( i ) < pFit( i )
            pFit( i ) = fit( i );
            pX( i, : ) = x( i, : );
        end
        
        if( pFit( i ) < fMin && rand < A(i) )
            fMin = pFit( i );
            bestX = pX( i, : );
            bestIndex = i;
            bestIter = iteration;
            
            A(i) = A(i) * alpha;
            r(i) = r0(i) * ( 1 - exp( -gamma * iteration ) );
        end
    end
    
    if( iteration - bestIter > G )                  
         r = rand( pop, 1 ) .* 0.05 + 0.85;
         A = rand( pop, 1 ) .* ( AMax - AMin ) + AMin;
    end

 end
function y = Sphere( x )
y = sum( x .^ 2 );

function s = Bounds( s, Lb, Ub)
  temp = s;
  I = temp < Lb;
  temp(I) = Lb(I);
  J = temp > Ub;
  temp(J) = Ub(J);
  s = temp;