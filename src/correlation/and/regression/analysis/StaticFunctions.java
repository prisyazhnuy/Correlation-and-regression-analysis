package correlation.and.regression.analysis;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import jsc.correlation.KendallCorrelation;
import jsc.datastructures.PairedData;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

import org.jfree.data.xy.XYSeries;

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
    
    public static double tDistrib(double nu){
        TDistribution distribution = new TDistribution(nu);
        return distribution.cumulativeProbability(1-0.05/2);
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
        return Math.log(Math.exp(numerator/denominator));
    }
    
    public static double estimateB(OrderedSeries X, OrderedSeries Y){
        double numerator = xy_mean(X, Y)-X.mean()*y_mean(Y);
        double denominator = X.meanSq()-Math.pow(X.mean(), 2);
        return numerator/denominator;
    }
    
    public static XYSeries drawRegressionLine(OrderedSeries X, OrderedSeries Y){
        double a = estimateA(X, Y);
        double b = estimateB(X, Y);
        //double point = X.min();
        //double max = X.max();
        //double step = (max-point)/100.0;
        XYSeries res = new XYSeries("Regresion");
       // while(point<=max){
            for(int i=0; i<X.size(); i++){
                res.add(X.getNumber(i), regression(X.getNumber(i), a, b));
            }
            //res.add(point, regression(point, a, b));
           // point+=step;
       // }
        res.setDescription("Regression Line");
        return res;
    }
    
    public static XYSeries drawConfidenceIntervalMin(OrderedSeries X, OrderedSeries Y){
        double a = estimateA(X, Y);
        double b = estimateB(X, Y);
        double st = Quantiles.Student(1-0.05/2, X.size()-2);
        //double point = X.min();
        //double max = X.max();
        //double step = (max-point)/100.0;
        XYSeries res = new XYSeries("ConfidenceIntervalForRegMin");
       // while(point<=max){
            for(int i=0; i<X.size(); i++){
                //res.add(X.getNumber(i), regression(X.getNumber(i), a, b)-st*getSy(X.getNumber(i), X, Y));
                res.add(X.getNumber(i), getZSmin(X.getNumber(i), a, b, st, getSy(X.getNumber(i), X, Y)));
            }
            //res.add(point, regression(point, a, b));
           // point+=step;
       // }
        res.setDescription("ConfidenceIntervalMin");
        return res;
    }
    
           private static double getZSmin(double x, double a, double b, double t, double S){
           return Math.exp(a+b*x-t*S);
       }
       
       private static double getZSmax(double x, double a, double b, double t, double S){
           return Math.exp(a+b*x+t*S);
       }
    
        public static XYSeries drawConfidenceIntervalMax(OrderedSeries X, OrderedSeries Y){
        double a = estimateA(X, Y);
        double b = estimateB(X, Y);
        double st = Quantiles.Student(1-0.05/2, X.size()-2);
        //double point = X.min();
        //double max = X.max();
        //double step = (max-point)/100.0;
        XYSeries res = new XYSeries("ConfidenceIntervalForRegMax");
       // while(point<=max){
            for(int i=0; i<X.size(); i++){
                //res.add(X.getNumber(i), regression(X.getNumber(i), a, b)+st*getSy(X.getNumber(i), X, Y));
                res.add(X.getNumber(i), getZSmax(X.getNumber(i), a, b, st, getSy(X.getNumber(i), X, Y)));
            }
            //res.add(point, regression(point, a, b));
           // point+=step;
       // }
        res.setDescription("ConfidenceIntervalForRegMax");
        return res;
    }
        
            public static XYSeries drawConfidence2IntervalMin(OrderedSeries X, OrderedSeries Y){
        double a = estimateA(X, Y);
        double b = estimateB(X, Y);
        double st = Quantiles.Student(1-0.05/2, X.size()-2);
        //double point = X.min();
        //double max = X.max();
        //double step = (max-point)/100.0;
        XYSeries res = new XYSeries("ConfidenceForNewValueIntervalMin");
       // while(point<=max){
            for(int i=0; i<X.size(); i++){
                res.add(X.getNumber(i), getZSymin(X.getNumber(i), a, b, st, getS(X.getNumber(i), X, Y)));
                //res.add(X.getNumber(i), regression(X.getNumber(i), a, b)-st*getS(X.getNumber(i), X, Y));
            }
            //res.add(point, regression(point, a, b));
           // point+=step;
       // }
        //res.setDescription("ConfidenceIntervalMin");
        return res;
    }
            private static double getZSymin(double x, double a, double b, double t, double S){
           return Math.exp(a+b*x-t*S);
       }
       
       private static double getZSymax(double x, double a, double b, double t, double S){
           return Math.exp(a+b*x+t*S);
       }
    
        public static XYSeries drawConfidence2IntervalMax(OrderedSeries X, OrderedSeries Y){
        double a = estimateA(X, Y);
        double b = estimateB(X, Y);
        double st = Quantiles.Student(1-0.05/2, X.size()-2);
        //double point = X.min();
        //double max = X.max();
        //double step = (max-point)/100.0;
        XYSeries res = new XYSeries("ConfidenceIntervalForNewValueMax");
    
       // while(point<=max){
            for(int i=0; i<X.size(); i++){
                //res.add(X.getNumber(i), regression(X.getNumber(i), a, b)+st*getS(X.getNumber(i), X, Y));
                res.add(X.getNumber(i), getZSymax(X.getNumber(i), a, b, st, getS(X.getNumber(i), X, Y)));
            }
            //res.add(point, regression(point, a, b));
           // point+=step;
       // }
       // res.setDescription("ConfidenceIntervalMax");
        return res;
    }
        
          public static XYSeries drawTolerantIntervalMax(OrderedSeries X, OrderedSeries Y){
        double a = estimateA(X, Y);
        double b = estimateB(X, Y);
        double st = Quantiles.Student(1-0.05/2, X.size()-2);
        double S = S_residual(X, Y);
        //double t = st*Math.sqrt(S_residual(X, Y));
        //double point = X.min();
        //double max = X.max();
        //double step = (max-point)/100.0;
        XYSeries res = new XYSeries("TolerantIntervalMax");
            for(int i=0; i<X.size(); i++){
                res.add(X.getNumber(i), getZmax(X.getNumber(i), a, b, st, S));
                //res.add(X.getNumber(i), regression(X.getNumber(i), a, b)+t);
            }
      
        res.setDescription("TolerantIntervalMax");
        return res;
    }
          
       private static double getZmin(double x, double a, double b, double t, double S){
           return Math.exp(a+b*x-t*Math.sqrt(S));
       }
       
       private static double getZmax(double x, double a, double b, double t, double S){
           return Math.exp(a+b*x+t*Math.sqrt(S));
       }
          
            public static XYSeries drawTolerantIntervalMin(OrderedSeries X, OrderedSeries Y){
        double a = estimateA(X, Y);
        double b = estimateB(X, Y);
        double st = Quantiles.Student(1-0.05/2, X.size()-2);
        double S = S_residual(X, Y);
        double t = st*Math.sqrt(S_residual(X, Y));
        //double point = X.min();
        //double max = X.max();
        //double step = (max-point)/100.0;
        XYSeries res = new XYSeries("TolerantIntervalMin");
       // while(point<=max){
            for(int i=0; i<X.size(); i++){
                res.add(X.getNumber(i), getZmin(X.getNumber(i),a,b,st, S));
                //res.add(X.getNumber(i), regression(X.getNumber(i), a, b)-t);
            }
            //res.add(point, regression(point, a, b));
           // point+=step;
       // }
        res.setDescription("TolerantIntervalMin");
        return res;
    }
    
    private static double getS(double x, OrderedSeries X, OrderedSeries Y){
        return Math.sqrt(S_residual(X, Y)*(1+1/(double)X.size())+Math.pow(Sb(X, Y), 2)*Math.pow(x-X.mean(), 2));
    }
    
    private static double getSy(double x, OrderedSeries X, OrderedSeries Y){
           
        return Math.sqrt(S_residual(X, Y)/(double)X.size()+Math.pow(Sb(X, Y), 2)*Math.pow(x-X.mean(), 2));//*(1+1/(double)X.size()+Math.pow(x-X.mean(), 2)/(double)X.size()*(X.meanSq()-Math.pow(X.mean(), 2)));
    }
    
    public static double Sa(OrderedSeries X, OrderedSeries Y){
      
        double res = Math.sqrt(S_residual(X, Y))*Math.sqrt((1/(double)X.size())+(Math.pow(X.mean(), 2)/(X.unbiasedEstimator()*(X.size()-1))));
        //System.out.println("Sa "+res);
        return res;
    }
    
    public static double Sb(OrderedSeries X, OrderedSeries Y){
             double res = Math.sqrt(S_residual(X, Y))/(Math.sqrt(X.unbiasedEstimator())*Math.sqrt(X.size()-1));
        //System.out.println("Sb "+res);
        return res;
    }
    
    public static double dispersiaA(OrderedSeries X, OrderedSeries Y){  
        return S_residual(X, Y)*(1.0/(double)X.size()+Math.pow(X.mean(), 2)/(X.size()*(X.meanSq()-Math.pow(X.mean(), 2))));
    }
    
    public static double dispersiaB(OrderedSeries X, OrderedSeries Y){  
        return S_residual(X, Y)/(double)(X.size()*(X.meanSq()-Math.pow(X.mean(), 2)));
    }
    
    public static String getIntervalForEstemateA(OrderedSeries X, OrderedSeries Y){
        NumberFormat f = NumberFormat.getInstance();
        f.setGroupingUsed(false);
        double a = estimateA(X, Y);
        double sec = Quantiles.Student(1-0.05/2, X.size()-2)*Sa(X, Y);
        double fir = a-sec;
        double s = a+sec;
        return new String("["+f.format(fir)+":"+f.format(s)+"]");
    }
    
    public static String getIntervalForEstemateB(OrderedSeries X, OrderedSeries Y){
        NumberFormat f = NumberFormat.getInstance();
        f.setGroupingUsed(false);
        double b = estimateB(X, Y);
        double sec = Quantiles.Student(1-0.05/2.0, X.size()-2)*Sb(X, Y);
        double fir = b-sec;
        double s = b+sec;
        return new String("["+f.format(fir)+":"+f.format(s)+"]");
    }
    
    public static double S_residual(OrderedSeries X, OrderedSeries Y){
        double sum = 0;
        double a = estimateA(X, Y);
        double b = estimateB(X, Y);
        for(int i=0; i<X.size(); i++){
            sum+=Math.pow(funcY(Y.getNumber(i))-a-b*X.getNumber(i), 2);
            //sum+=Math.pow(Y.getNumber(i)-a-b*X.getNumber(i), 2);
        }
        //System.out.println("S_resideul "+sum/(double)(X.size()-2));
        return roundResult(sum/(double)(X.size()-2), 6);
    }
    
    public static String significanceA(OrderedSeries X, OrderedSeries Y){
        if(statisticA(X, Y)<=Quantiles.Student(1-0.05/2, X.size()-2)){
            return "Не значимый";
        }else{
            return "Значимый";
        }
    }
    
    public static double statisticA(OrderedSeries X, OrderedSeries Y){
        return Math.abs(estimateA(X, Y))/Sa(X, Y);
    }
    
    public static double statisticB(OrderedSeries X, OrderedSeries Y){
        return Math.abs(estimateB(X, Y))/Sb(X, Y);
    }
    
    public static String significanceB(OrderedSeries X, OrderedSeries Y){
        if(statisticB(X, Y)<=Quantiles.Student(1-0.05/2, X.size()-2)){
            return "Не значимый";
        }else{
            return "Значимый";
        }
    }
    
    public static double regression(double x, double a, double b){
        //return a*Math.exp(b*x);
        return Math.exp(a+b*x);
    }
    
    private static double funcY(double x){
        //return Math.exp(x);
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
    
    public static double getFtest(OrderedSeries X, OrderedSeries Y){
        double sum=0;
        double a = estimateA(X, Y);
        double b = estimateB(X, Y);
        double y_m = Y.mean();
        for(int i=0; i<X.size(); i++){
            sum+=Math.pow(regression(X.getNumber(i), a, b)-y_m, 2);
        }
        return sum/S_residual(X, Y);
    }
    
    
    static double roundResult (double d, int precise) {
        precise = 10^precise;
        d = d*precise;
        int i = (int) Math.round(d);
        return (double) i/precise;
    }
    
    
    
}
