package correlation.and.regression.analysis;

import java.text.NumberFormat;
import java.util.ArrayList;
import jsc.correlation.KendallCorrelation;
import jsc.datastructures.PairedData;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

public class StaticFunctions {
   // static specialFunctions spFunc = new specialFunctions();
    
    public static double pairCorrelationCoef(OrderedSeries X, OrderedSeries Y){
        double xMean = X.mean();
        double yMean = Y.mean();
        double xyMean = xyMean(X, Y);
        double xSquareMean = X.meanSquare(xMean);
        double ySquareMean = Y.meanSquare(yMean);
        return (xyMean-(xMean*yMean))/(xSquareMean*ySquareMean);
    }
    
    public static double significancePairCor(OrderedSeries X, OrderedSeries Y){
        double N = X.size();
        double pairCor = pairCorrelationCoef(X, Y);
        return (pairCor*Math.pow(N-2, 0.5))/(Math.pow(1-Math.pow(pairCor, 2), 0.5));
    }
    
    public static double markCorrelationRatio(OrderedSeries X, OrderedSeries Y){
        ArrayList<ArrayList<Integer>> indexs = getArray(X, Y);
        double numerator = 0;
        for(int i=0; i<indexs.size(); i++){
            numerator+=indexs.get(i).size()*Math.pow(Yi_Mean(Y, indexs.get(i))-Y_Mean(Y, indexs), 2);
        }
        
        double denominator=0;
        for(int i=0; i<indexs.size(); i++){
            for(int j=0; j<indexs.get(i).size(); j++){
                denominator+=Math.pow(Y.getNumber(indexs.get(i).get(j))-Y_Mean(Y, indexs), 2);
            }
        }
        
        return Math.pow(numerator/denominator, 0.5);
    }
    
    public static int getK(int N){
        return (int)Math.round(1+3.2*Math.log10(N));
    }
    
    public static double statisticsCorrelationRatio(double markCorr, OrderedSeries X){
        int k=(int) (1+3.2*Math.log10(X.size()));
        int N=X.size();
        double numerator = Math.pow(markCorr, 2)*(N-k);
        double denominator = (k-1)*(1-Math.pow(markCorr, 2));
        return numerator/denominator;
    }
    
    public static double spirmenCoef(double[] X, double[] Y){
        SpearmansCorrelation correlation = new SpearmansCorrelation();
        return correlation.correlation(X, Y);
    }
    
    public static double statisticsSpirmen(double[] X, double[] Y){
        double spirmen = spirmenCoef(X, Y);
        return (spirmen*Math.pow(X.length-2, 0.5))/(Math.pow(1-Math.pow(spirmen, 2), 0.5));
    }
    
    public static double kendallsCoef(double[] X, double[] Y){
        KendallCorrelation correlation = new KendallCorrelation(new PairedData(X, Y)); 
        return correlation.getR();
    }
    
    public static double staticticsKendall(double[] X, double[] Y){
        int N = X.length;
        double kendall = kendallsCoef(X, Y);
        return (3*kendall*Math.pow(N*N-N, 0.5))/(Math.pow(4*N+10, 0.5));
    }
    
    public static String getIntervalForPairCorr(double r, int N){
        NumberFormat f = NumberFormat.getInstance();
        f.setGroupingUsed(false);
        double res = r+(r-Math.pow(r, 3))/(2*N);
        double sec = Quantiles.Norm(1-0.05/2)*(1-r*r)/(Math.sqrt(N-1));
        double fir = res-sec;
        double s = res+sec;
        return new String("["+f.format(fir)+":"+f.format(s)+"]");
    }
    
    public static double estimateA(OrderedSeries X, OrderedSeries Y){
        double numerator = y_mean(Y)*X.meanSq()-X.mean()*xy_mean(X, Y);
        double denominator = X.meanSq()-Math.pow(X.mean(), 2);
        return numerator/denominator;
    }
    
    public static double estimateB(OrderedSeries X, OrderedSeries Y){
        double numerator = xy_mean(X, Y)-X.mean()*y_mean(Y);
        double denominator = X.meanSq()-Math.pow(X.mean(), 2);
        return numerator/denominator;
    }
    
    private static double funcY(double x){
        return Math.log(x);
    }
    
    private static double y_mean(OrderedSeries Y){
        double res = 0;
        for(int i=0; i<Y.size(); i++){
            res+=funcY(Y.getNumber(i));
        }
        return res/Y.size();
    }
    
    private static double xy_mean(OrderedSeries X, OrderedSeries Y){
        double res=0;
        for(int i=0; i<X.size(); i++){
            res+=X.getNumber(i)*funcY(Y.getNumber(i));
        }
        return res/X.size();
    }
    
    private static double Yi_Mean(OrderedSeries Y, ArrayList<Integer> index){
        double res=0;
        for(Integer i : index){
            res+=Y.getNumber(i);
        }
        return res/index.size();
    }
    
    private static double Y_Mean(OrderedSeries Y, ArrayList<ArrayList<Integer>> indexs){
        double res=0;
        for(ArrayList<Integer> index : indexs){
            res+=index.size()*Yi_Mean(Y, index);
        }
        return res/Y.size();
    }
    
    private static double xyMean(OrderedSeries X, OrderedSeries Y){
        double result=0;
        double N = X.size();
        for(int i=0; i<N; i++){
            result+=X.getNumber(i)*Y.getNumber(i);
        }
        return result/N;
    }
    
    private static ArrayList<ArrayList<Integer>> getArray(OrderedSeries X, OrderedSeries Y){
        int k=(int) (1+3.2*Math.log10(X.size()));
        //System.out.println(""+k);
        ArrayList<ArrayList<Integer>> res = new ArrayList<>();
        double x_min=X.min();
        //System.out.println(""+x_min);
        double x_max=X.max();
        //System.out.println(""+x_max);
        double h=(x_max-x_min)/k;
        //System.out.println(""+h);
        for(int i=0; i<k; i++){
            ArrayList<Integer> row = new ArrayList<>();
            for(int j=0; j<X.size(); j++){
                if((X.getNumber(j)>=x_min+i*h)&&(X.getNumber(j)<=x_min+(i+1)*h)){
                    row.add(j);
                }
            }
            res.add(row);
        }
        return res;
    }
    
    
    
}
