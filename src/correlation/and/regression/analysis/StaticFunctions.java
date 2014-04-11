package correlation.and.regression.analysis;

public class StaticFunctions {
    public static double pairCorrelationCoef(OrderedSeries X, OrderedSeries Y){
        double xMean = X.mean();
        double yMean = Y.mean();
        double xyMean = xyMean(X, Y);
        double xSquareMean = X.meanSquare(xMean);
        double ySquareMean = Y.meanSquare(yMean);
        return (xyMean-(xMean*yMean))/(xSquareMean*ySquareMean);
    }
    
    public static double significancePairCor(OrderedSeries X, OrderedSeries Y){
        double N = X.countOfNumbers;
        double pairCor = pairCorrelationCoef(X, Y);
        return (pairCor*Math.pow(N-2, 0.5))/(Math.pow(1-Math.pow(pairCor, 2), 0.5));
    }
    
    private static double xyMean(OrderedSeries X, OrderedSeries Y){
        double result=0;
        double N = X.countOfNumbers;
        for(int i=0; i<N; i++){
            result+=X.getNumber(i)*Y.getNumber(i);
        }
        return result/N;
    }
}
