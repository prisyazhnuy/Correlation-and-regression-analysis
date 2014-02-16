package correlation.and.regression.analysis;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Vector;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;


public class OrderedSeries {
    
    ArrayList<Count> data = new ArrayList<>();
    Vector<Double> array = new Vector<>();
    int countOfNumbers = 0;
    static final int countCol=10;
    static final int countRow=5;
    
    private double D(int index){
        double max=0;
        for(int i=0; i<data.size(); i++){
            double func;
            if(i-index<0){
                func=0;
            }else
            if(data.get(i).number<markA()){
                func=0;
            }else{
                if(data.get(i).number>=markB()){
                    func=1;
                }else{
                    func=(data.get(i).number-markA())/(markB()-markA());
                }
            }
            if(Math.abs(data.get(i).emprDistrib-func)>max){
                max=Math.abs(data.get(i).emprDistrib-func);
            }
        }
        System.out.println("Max= "+max);
        return max;
    }
    
    private double Z(){
        double max=D(0);
        if(D(1)>max){
            max=D(1);
        }
        return Math.pow(countOfNumbers, 0.5)*max;
    }
    
    private double f1(int k){
        return Math.pow(k, 2)-0.5*(1-Math.pow(-1, k));
    }
    
    private double f2(int k){
        return 5*Math.pow(k, 2)+22-7.5*(1-Math.pow(-1, k));
    }
    
    private double K(){
        double z=Z();
        double sum=0;
        for(int k=0;k<6;k++){
            sum+=Math.pow(-1, k)*Math.exp(-2*Math.pow(k, 2)*Math.pow(z, 2))*(1-
                    (2*k*k*z)/(3*Math.pow(countOfNumbers, 0.5))-(1/(double)(18*countOfNumbers))*
                    ((f1(k)-4*(f1(k)+3))*k*k*z*z+8*k*k*k*k*z*z*z*z)+(k*k*z/(27*
                    Math.pow(countOfNumbers, (double)3/2)))*(((f2(k)*f2(k))/5)-
                    (4*(f2(k)+45)*k*k*z*z)/(15)+8*k*k*k*k*z*z*z*z));
        }
        return 1+2*sum;
    }
    
    private double P(){
        return 1-K();
    }
    
    public DefaultTableModel showKolmogorov(DefaultTableModel tm){
        NumberFormat f = NumberFormat.getInstance();
        f.setGroupingUsed(false);
        tm.setValueAt(f.format(Z()), 0, 0);
        tm.setValueAt(f.format(K()), 0, 1);
        tm.setValueAt(f.format(P()), 0, 2);
        //tm.addRow(new String[]{""+K(), ""+P()});
        return tm;
    }
    
    public double[] toArrayDouble(int countAnom){
//        int myNULL=0;
//        for(Count c : data){
//            if(c==null){
//                myNULL++;
//            }
//        }
//        double[] result = new double[countOfNumbers-myNULL];
//        int index = 0;
//        for(Count c : data){
//            if(c!=null){
//            for(int i=0; i<c.count; i++){
//                
//                result[index]=c.number;
//                index++;
//                }
//            }
//        }
        double[] result = new double[countOfNumbers-countAnom];
        int index = 0;
        for(Count c : data){
            for(int i=0; i<c.count; i++){
                
                result[index]=c.number;
                index++;
            }
        }
        return result;
    }
    
     public JFreeChart createStatFunkDensity(float a, float b){
        XYSeriesCollection dataset = new XYSeriesCollection();
        int t=0;
         XYSeries tmp = new XYSeries("1");
        for(Count c : data){
           
            if((c.number<a)||(c.number>=b)){
                tmp.add(c.number, 0);
            }
            else{
                tmp.add(c.number, (double)1/(double)(b-a));
            }
            //tmp.add(c.number, c.emprDistrib);
            //tmp.add(c.number, c.emprDistrib);
            
            t++;
        }
        dataset.addSeries(tmp);
        JFreeChart chart = ChartFactory.createXYLineChart(
            "Static Density", "X", "Fn(x)", dataset, PlotOrientation.VERTICAL, false, false, false);
        XYPlot plot = chart.getXYPlot();
        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
        //for(int i=0;i<t;i++)
            renderer.setSeriesPaint(0, Color.BLUE);
        plot.setBackgroundPaint(Color.WHITE);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
        return chart;
    }
     
      public JFreeChart createStatFunkDestrib(float a, float b){
        XYSeriesCollection dataset = new XYSeriesCollection();
         XYSeries tmp = new XYSeries("1");
        for(Count c : data){
           
            if(c.number<a){
                tmp.add(c.number, 0);
            }
            else{
                if(c.number>=b){
                    tmp.add(c.number, 1);
                }else
                tmp.add(c.number, (double)(c.number-a)/(double)(b-a));
            }
        }
        dataset.addSeries(tmp);
        JFreeChart chart = ChartFactory.createXYLineChart(
            "Static Distribution", "X", "Fn(x)", dataset, PlotOrientation.VERTICAL, false, false, false);
        XYPlot plot = chart.getXYPlot();
        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
            renderer.setSeriesPaint(0, Color.BLUE);
        plot.setBackgroundPaint(Color.WHITE);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
        return chart;
    }
    
    public double markA(){
        double result = 3*(meanSq()-Math.pow(mean(), 2));
        return mean()-Math.pow(result, 0.5);
    }
    
    public double markB(){
        double result = 3*(meanSq()-Math.pow(mean(), 2));
        return mean()+Math.pow(result, 0.5);
    }

    public void loadFromFile(String path) throws FileNotFoundException, IOException {
        BufferedReader buffile;
        buffile = new BufferedReader(new FileReader(path));
        String tmp;                        //рабочая строка
        Boolean flag = false;
        while(buffile.ready()){
            tmp=buffile.readLine();
            String[] str = tmp.split("\\s+");
            for(int i=0; i<str.length; i++){
                array.add(Double.parseDouble(str[i]));
                for(int j=0; j<=data.size(); j++){
                    if(j==data.size()){
                        data.add(new Count(Float.parseFloat(str[i])));
                        countOfNumbers++;
                        break;
                    }else{
                        if(data.get(j).number==Float.parseFloat(str[i])){
                            data.get(j).count++;
                            countOfNumbers++;
                            break;
                        }
                    }
                }
                
            }
        }
        sort(0,data.size()-1);
        relatFreg();
        emprDistrib();
    }
    
    public void loadFromArray(double[] arr) {
        
            for(int i=0; i<arr.length; i++){
                for(int j=0; j<=data.size(); j++){
                    if(j==data.size()){
                        data.add(new Count((float) arr[i]));
                        countOfNumbers++;
                        break;
                    }else{
                        if(data.get(j).number==arr[i]){
                            data.get(j).count++;
                            countOfNumbers++;
                            break;
                        }
                    }
                }
                
            }
        
        sort(0,data.size()-1);
        relatFreg();
        emprDistrib();
    }
    
    public JFreeChart createGraph(){
        XYSeriesCollection dataset = new XYSeriesCollection();
        int t=0;
        for(Count c : data){
            XYSeries tmp = new XYSeries(""+c.number);
            tmp.add(c.number, c.emprDistrib);
            tmp.add(c.number, c.emprDistrib);
            dataset.addSeries(tmp);
            t++;
        }
        
        XYSeries tmp = new XYSeries("tmp");
        for(Count c : data){
           
            if(c.number<markA()){
                tmp.add(c.number, 0);
            }
            else{
                if(c.number>=markB()){
                    tmp.add(c.number, 1);
                }else
                tmp.add(c.number, (double)(c.number-markA())/(double)(markB()-markA()));
            }
        }
        dataset.addSeries(tmp);
        
        JFreeChart chart = ChartFactory.createXYLineChart(
            "Empire for series", "X", "Fn(x)", dataset, PlotOrientation.VERTICAL, false, false, false);
        XYPlot plot = chart.getXYPlot();
        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
        for(int i=0;i<t;i++)
            renderer.setSeriesPaint(i, Color.BLUE);
        renderer.setSeriesPaint(t, Color.RED);
        plot.setBackgroundPaint(Color.WHITE);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
        return chart;
    }
    
    public JFreeChart buildGrid(){
        double Xmin = data.get(0).number;
        double Xmax = data.get(data.size()-1).number;
        double Fmin = data.get(0).emprDistrib;
        double Fmax = data.get(data.size()-1).emprDistrib;
        double stepCol = (Xmax-Xmin)/(double)countCol;
        double stepRow = (Fmax-Fmin)/(double)countRow;
        
        XYSeriesCollection dataset = new XYSeriesCollection();
        double t=Xmin;
        for(int c=0; c<countCol; c++){
            XYSeries tmp = new XYSeries(""+c);
            tmp.add(t, 0);
            tmp.add(t, 1);
            dataset.addSeries(tmp);
            t+=stepCol;
        }
        
        t=Fmin;
        for(int r=0; r<countRow; r++){
            XYSeries tmp = new XYSeries("A"+r);
            tmp.add(Xmin, t);
            tmp.add(Xmax, t);
            dataset.addSeries(tmp);
            t+=stepRow;
        }
        int k=0;
        for(Count c : data){
            XYSeries tmp = new XYSeries("F"+c.number);
            tmp.add(c.number, c.emprDistrib);
            tmp.add(c.number, c.emprDistrib);
            dataset.addSeries(tmp);
            k++;
        }
        
        JFreeChart chart = ChartFactory.createXYLineChart(
            "Grid", "X", "Fn(x)", dataset, PlotOrientation.VERTICAL, false, false, false);
        XYPlot plot = chart.getXYPlot();
        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
        for(int i=0;i<15;i++)
            renderer.setSeriesPaint(i, Color.RED);
        for(int i=15;i<15+t;i++)
            renderer.setSeriesPaint(i, Color.BLUE);
        plot.setBackgroundPaint(Color.WHITE);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
        return chart;
    }
    
    public void showOrderedSeries(DefaultTableModel model){
        NumberFormat f = NumberFormat.getInstance();
        f.setGroupingUsed(false);
        int n = model.getRowCount();
        for(int j=0; j<n; j++)
            model.removeRow(0);
        int i=0;
        for(Count c : data){
            model.addRow(new String[] {String.valueOf(i+1), String.valueOf(c.number),String.valueOf(c.count), f.format(c.relativeFreq), f.format(c.emprDistrib)});
            i++;
        }
    }
    
    public TableModel showMark(TableModel tm, double alpha){
        NumberFormat f = NumberFormat.getInstance();
        f.setGroupingUsed(false);
        double mean = mean();
         tm.setValueAt("a", 0, 0);
        tm.setValueAt("b", 1, 0);
        tm.setValueAt(f.format(markA()), 0, 1);
        tm.setValueAt(f.format(markB()), 1, 1);
        tm.setValueAt(markA()-Quantiles.GetStudentsQuantile(alpha, countOfNumbers-2)*variance(markA(),markB()), 0, 2);
        tm.setValueAt(markB()-Quantiles.GetStudentsQuantile(alpha, countOfNumbers-2)*variance(markA(),markB()), 1, 2);
        tm.setValueAt(markA()+Quantiles.GetStudentsQuantile(alpha, countOfNumbers-2)*variance(markA(),markB()), 0, 3);
        tm.setValueAt(markB()+Quantiles.GetStudentsQuantile(alpha, countOfNumbers-2)*variance(markA(),markB()), 1, 3);
        return tm;
    }
    
    private double variance(double a, double b){
        return 2*Math.pow((b-a), 2)/(double)(15*countOfNumbers);
    }
    
    public TableModel showStatistics(TableModel tm, double alpha){
        NumberFormat f = NumberFormat.getInstance();
        f.setGroupingUsed(false);
        double mean = mean();
        tm.setValueAt(f.format(mean), 0, 1);
        tm.setValueAt(f.format(mediana()), 1, 1);
        tm.setValueAt(f.format(meanSquare(mean)), 2, 1);
        tm.setValueAt("-", 1, 2);
        tm.setValueAt("-", 1, 3);
        tm.setValueAt(f.format(meanSquareDeviationArithm(mean)), 0, 2);
        tm.setValueAt(f.format(meanSquareDeviationSquare(mean)), 2, 2);
        tm.setValueAt(f.format(asymmetryCoefficient(mean)), 3, 1);
        tm.setValueAt(f.format(ekstsessaCoefficient(mean)), 4, 1);
        tm.setValueAt(f.format(kontrekstsessaCoefficient(mean)), 5, 1);
        tm.setValueAt(f.format(variationCoefficient(mean)), 6, 1);
        tm.setValueAt(f.format(meanSquareDeviationAsymmetry()), 3, 2);
        tm.setValueAt(f.format(meanSquareDeviationEkstsessa()), 4, 2);
        tm.setValueAt(f.format(meanSquareDeviationKontrekstsessa(mean)), 5, 2);
        tm.setValueAt(f.format(meanSquareDeviationVariation(mean)), 6, 2);
        tm.setValueAt(""+f.format(confidenceIntervalMin(mean, meanSquareDeviationArithm(mean), alpha))+" : "+f.format(confidenceIntervalMax(mean, meanSquareDeviationArithm(mean), alpha)), 0, 3);
        tm.setValueAt(""+f.format(confidenceIntervalMin(meanSquare(mean), meanSquareDeviationSquare(mean), alpha))+" : "+f.format(confidenceIntervalMax(meanSquare(mean), meanSquareDeviationSquare(mean), alpha)), 2, 3);
        tm.setValueAt(""+f.format(confidenceIntervalMin(asymmetryCoefficient(mean), meanSquareDeviationAsymmetry(), alpha))+" : "+f.format(confidenceIntervalMax(asymmetryCoefficient(mean), meanSquareDeviationAsymmetry(), alpha)), 3, 3);
        tm.setValueAt(""+f.format(confidenceIntervalMin(ekstsessaCoefficient(mean), meanSquareDeviationEkstsessa(), alpha))+" : "+f.format(confidenceIntervalMax(ekstsessaCoefficient(mean), meanSquareDeviationEkstsessa(), alpha)), 4, 3);
        tm.setValueAt(""+f.format(confidenceIntervalMin(kontrekstsessaCoefficient(mean), meanSquareDeviationKontrekstsessa(mean), alpha))+" : "+f.format(confidenceIntervalMax(kontrekstsessaCoefficient(mean), meanSquareDeviationKontrekstsessa(mean), alpha)), 5, 3);
        tm.setValueAt(""+f.format(confidenceIntervalMin(variationCoefficient(mean), meanSquareDeviationVariation(mean), alpha))+" : "+f.format(confidenceIntervalMax(variationCoefficient(mean), meanSquareDeviationVariation(mean), alpha)), 6, 3);
        return tm;
    }
     
    public int deleteAnomaly(Classes cl, double gamma1, double gamma2){
        float[] range = cl.anomaly(gamma1, gamma2);
        //int i=0;
        int count = 0;
        for(int i=0; i<data.size(); i++){
            if(data.get(i).number<range[0]){
                count=count+data.get(i).count;
                data.remove(i);
                i--;
            }else{
                if(data.get(i).number>range[1]){
                    count=count+data.get(i).count;
                    data.remove(i);
                    i--;
                }
            }
        }
//        while(data.get(i).number<range[0]){
//            //data.set(i, null);
//            data.remove(i);
//            i++;
//        }
//        for(int t=i; t<data.size(); t++){
//            if(data.get(t).number>range[1]){
//                //data.set(t, null);
//            }
//        }
        return count;
    }
        
    private double confidenceIntervalMin(double statistic,double standardDeviation,double alpha){
        return statistic - Quantiles.GetStudentsQuantile(alpha, countOfNumbers -2) * standardDeviation;
    }
    
    private double confidenceIntervalMax(double statistic,double standardDeviation,double alpha){
        return statistic + Quantiles.GetStudentsQuantile(alpha, countOfNumbers -2) * standardDeviation;
    }
    
    private double mean(){
        double middle = 0;
        int i=0;
        int j=data.size()-1;
        while(i<=j){
            if(i==j){
                middle+=data.get(i).number*data.get(i).count;
            }
            else{
                middle+=data.get(i).number*data.get(i).count;
                middle+=data.get(j).number*data.get(j).count;
            }
            i++;
            j--;
        }
        return middle/(double)countOfNumbers;
    }
    
    private double meanSq(){
        double middle = 0;
        int i=0;
        int j=data.size()-1;
        while(i<=j){
            if(i==j){
                middle+=Math.pow(data.get(i).number*data.get(i).count,2);
            }
            else{
                middle+=Math.pow(data.get(i).number*data.get(i).count,2);
                middle+=Math.pow(data.get(j).number*data.get(j).count,2);
            }
            i++;
            j--;
        }
        return middle/(double)countOfNumbers;
    }
    
    private double meanSquare(double mean){
        double result = 0;
        for(Count c : data){
            result+=Math.pow(c.number-mean,2)*c.count;
        }
        return Math.pow(result/(double)(countOfNumbers-1),0.5);
    }
    
    private double meanSquareDeviationArithm(double mean){
        return meanSquare(mean)/Math.pow(countOfNumbers, 0.5);
    }
    
    private double meanSquareDeviationSquare(double mean){
        return meanSquare(mean)/Math.pow(2*countOfNumbers, 0.5);
    }
    
    private double meanSquareDeviationVariation(double mean){
         double variationCoe = variationCoefficient(mean);
         return variationCoe*Math.pow(((1+2*Math.pow(variationCoe, 2))/(double)(2*countOfNumbers)), 0.5);
     }
    
    private double variance(double mean){
        double result = 0;
        for(Count c : data){
            result+=Math.pow(c.number, 2)*c.count;
        }
        return Math.pow(result/(double)countOfNumbers-Math.pow(mean, 2), 0.5);
    }
    
    private double asymmetryCoefficient(double mean){
        double result = 0;
        for(Count c : data){
            result+=Math.pow(c.number - mean, 3)*c.count;
        }
        return (result/(double)(Math.pow(variance(mean), 3)*countOfNumbers))*(Math.pow(countOfNumbers*(countOfNumbers-1), 0.5)/(double)(countOfNumbers-2));
    }
    
    private double ekstsessaCoefficientOffset(double mean){
        double result = 0;
        for(Count c : data){
            result+=Math.pow(c.number-mean, 4)*c.count;
        }
        return result/(Math.pow(variance(mean), 4)*countOfNumbers);
    }
    
    private double ekstsessaCoefficient(double mean){
        double result = (Math.pow(countOfNumbers, 2)-1)/(double)((countOfNumbers-2)*(countOfNumbers-3));
        return result*((ekstsessaCoefficientOffset(mean)-3)+6/(double)(countOfNumbers+1));
    }
    
    private double kontrekstsessaCoefficient(double mean){
        return Math.pow(Math.abs(ekstsessaCoefficient(mean)), -0.5);
    }
    
    private double variationCoefficient(double mean){   
        return meanSquare(mean)/mean;
    }
    
    private int less(double mediana){
        int count=0;
        for(Count c : data){
            if(c.number<mediana){
                count+=c.count;
            }
            else return count;
        }
        return count;
    }
    
    private int more(double mediana){
        int count=0;
        for(Count c : data){
            if(c.number>mediana){
                //count+=c.count;
                return countOfNumbers-count;
            }
            else count+=c.count;
        }
        return countOfNumbers-count;
    }
    
    public boolean odnorodnost2(DefaultTableModel model){
        int N1M = (int) model.getValueAt(0, 1);
        int N2M = (int) model.getValueAt(1, 1);
        int N1B = (int) model.getValueAt(0, 2);
        int N2B = (int) model.getValueAt(1, 2);
        int NM = (int) model.getValueAt(2, 1);
        int NB = (int) model.getValueAt(2, 2);
        int N1 = (int) model.getValueAt(0, 3);
        int N2 = (int) model.getValueAt(1, 3);
        int N = (int) model.getValueAt(2, 3);
        double U = (Math.abs(N1M*N2B-N2M*N1B)-((double)N/(double)2))/(Math.pow(NM*NB*N1*N2, 0.5));
        U=Math.abs(U*Math.pow(N, 0.5));
        return U<=Distributions.NormInv(1-0.05/2.0,0,1);
    }
    
    private double mediana(){
        if(countOfNumbers%2==0){
            int count = 0;
            for(int i=0; i<data.size()-1; i++){
                int precount=count;
                count+=data.get(i).count;
                if(count==countOfNumbers/2){
                    return (data.get(i).number+data.get(i+1).number)/(double)2;
                }
                if(precount<countOfNumbers/2&&countOfNumbers/2<count){
                    return (data.get(i).number+data.get(i+1).number)/(double)2;
                }
                
            }
        }else{
            int count = 0;
            for(Count c : data){
                int precount=count;
                count+=c.count;
                if(count == (countOfNumbers+1)/2){
                    return c.number;
                }
                 if(precount<(countOfNumbers+1)/2&&(countOfNumbers+1)/2<count){
                    return c.number;
                }
            }
        }
        return 2;
    }
    
    private double meanSquareDeviationAsymmetry(){
        return Math.pow((6*countOfNumbers-12)/(double)((countOfNumbers+1)*(countOfNumbers+3)), 0.5);
    }
            
    private double meanSquareDeviationEkstsessa(){
        double result = 1-(225/(double)(15*countOfNumbers+124));
        return Math.pow(24*result/countOfNumbers, 0.5);
    }
    
    private double meanSquareDeviationKontrekstsessa(double mean){
        double factor_1 = Math.pow(Math.abs(ekstsessaCoefficientOffset(mean)/(double)(29*countOfNumbers)), 0.5);
        double factor_2 = Math.pow(Math.pow(Math.abs(Math.pow(ekstsessaCoefficientOffset(mean), 2)-1), 3), 0.25);
        return factor_1*factor_2;
    }        
            
    private void relatFreg(){
        for(Count c : data){
            c.relativeFreq=(double)c.count/(double)countOfNumbers;
        }
    }
    
    private void emprDistrib(){
        double tmp = 0;
        for(Count c : data){
            tmp+=c.relativeFreq;
            c.emprDistrib = tmp;
        }
    }

    private void sort(int start, int end) {
        if (start >= end)
            return;
	int i = start, j = end;
	int cur = i - (i - j) / 2;
	while (i < j) {
            while (i < cur && (data.get(i).number <= data.get(cur).number)) {
                i++;
            }
	while (j > cur && (data.get(cur).number <= data.get(j).number)) {
            j--;
	}
	if (i < j) {
            Count temp = data.get(i);
            data.set(i, data.get(j));
            data.set(j, temp);
            if (i == cur)
                cur = j;
            else if (j == cur)
                cur = i;
        }
        }
        sort(start, cur);
        sort(cur+1, end);
    }
    
    public int size() {
        return data.size();
    }

    public float getNumber(int i) {
        return data.get(i).number;
    }

    public int getCount(int i) {
        return data.get(i).count;
    }

    public double getRelativeFreg(int i) {
        return data.get(i).relativeFreq;
    }

    public double getEmpiricalDistribution(int i) {
        return data.get(i).emprDistrib;
    }

    public int getCountOfNumbers() {
        return countOfNumbers;
    }

    public void addCount(float f, int count, double relFreg, double emprDis) {
        data.add(new Count(f, count, relFreg, emprDis));
    }

    void clear() {
        data.clear();
    }

    DefaultTableModel getMedians(DefaultTableModel modelThird, int index) {
        if(index==1){
            modelThird.setValueAt(less(mediana()), 0, 1);
            modelThird.setValueAt(more(mediana()), 0, 2);
            modelThird.setValueAt(countOfNumbers, 0, 3);
            if(!modelThird.getValueAt(1, 1).equals("")){
                int m=(int) modelThird.getValueAt(0, 1)+(int) modelThird.getValueAt(1, 1);
                modelThird.setValueAt(m, 2, 1);
                m=(int) modelThird.getValueAt(0, 2)+(int) modelThird.getValueAt(1, 2);
                modelThird.setValueAt(m, 2, 2);
                m=(int) modelThird.getValueAt(0, 3)+(int) modelThird.getValueAt(1, 3);
                modelThird.setValueAt(m, 2, 3);
            }
        }
        else{
            modelThird.setValueAt(less(mediana()), 1, 1);
            modelThird.setValueAt(more(mediana()), 1, 2);
            modelThird.setValueAt(countOfNumbers, 1, 3);
            if(!modelThird.getValueAt(0, 1).equals("")){
                int m=(int) modelThird.getValueAt(0, 1)+(int) modelThird.getValueAt(1, 1);
                modelThird.setValueAt(m, 2, 1);
                m=(int) modelThird.getValueAt(0, 2)+(int) modelThird.getValueAt(1, 2);
                modelThird.setValueAt(m, 2, 2);
                m=(int) modelThird.getValueAt(0, 3)+(int) modelThird.getValueAt(1, 3);
                modelThird.setValueAt(m, 2, 3);
            }
        }
        return modelThird;
    }

    void showArray(DefaultTableModel modelFridman) {
        modelFridman.addColumn("", array);
    }

    void showStatistics(DefaultTableModel modelStatistic) {
        NumberFormat f = NumberFormat.getInstance();
        f.setGroupingUsed(false);
        Vector<String> vect = new Vector<>();
        double mean = mean();
        vect.add(f.format(mean));
        vect.add(f.format(meanSquare(mean)));
        vect.add(f.format(asymmetryCoefficient(mean)));
        vect.add(f.format(ekstsessaCoefficient(mean)));
        vect.add(f.format(variationCoefficient(mean)));
        modelStatistic.addColumn(null, vect);
    }

    
    
}
class Count{
    int count = 0;
    float number;
    double relativeFreq;
    double emprDistrib;
    
    Count(float num){
        number = num;
        count = 1;
    }
    
    Count(float f, int count, double relFreg, double emprDis){
        number = f;
        this.count = count;
        relativeFreq = relFreg;
        emprDistrib = emprDis;
    }
}