package correlation.and.regression.analysis;

import java.awt.Color;
import javax.swing.table.DefaultTableModel;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.general.AbstractDataset;
import org.jfree.data.general.Dataset;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;
import org.jfree.data.xy.DefaultTableXYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class Classes {
    private double step;
    private OrderedSeries data;
    private int M = 0;
    private OrderedSeries classes;
    
    public JFreeChart createGraph(double a, double b){
        XYSeriesCollection dataset = new XYSeriesCollection();
        for(int i=0; i<M; i++){
            XYSeries tmp = new XYSeries(""+i);
            tmp.add(classes.getNumber(i)-step, classes.getEmpiricalDistribution(i));
            tmp.add(classes.getNumber(i), classes.getEmpiricalDistribution(i));
            dataset.addSeries(tmp);
        }
        
         XYSeries tmp = new XYSeries("tmp");
        for(Count c : data.data){
           
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
            "Empire for classes", "X", "Fn(x)", dataset, PlotOrientation.VERTICAL, false, false, false);
        XYPlot plot = chart.getXYPlot();
        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
        for(int i=0;i<M;i++)
            renderer.setSeriesPaint(i, Color.BLUE);
        renderer.setSeriesPaint(M, Color.RED);
        plot.setBackgroundPaint(Color.WHITE);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
        return chart;
    }
    
    public float[] anomaly(double gamma1, double gamma2){
        float[] result = new float[2];
        double anom = 0;
        result[0] = (float) (classes.getNumber(0)-step);
        for(Count c : classes.data){
             if((anom+c.emprDistrib)<=gamma1){
                 result[0] = c.number;
                 anom+= c.emprDistrib;
             }else {
                 break;
             }
        }
             anom = 0;
             result[1] = classes.getNumber(classes.size()-1);
             float sumEmpr = (float) classes.getRelativeFreg(classes.size()-1);
        for(int i = classes.size()-2; i>0; i--){
            if(anom+sumEmpr<=gamma2){
                result[1] = (float) (classes.getNumber(i+1)-step);
                anom+=sumEmpr;
                sumEmpr+=classes.getRelativeFreg(i);
            }else{
                break;
            }
        }
        
        return result;
    }
    
    public void showAnomaly(DefaultTableModel model, double gamma1, double gamma2){
        int n = model.getRowCount();
        for(int j=0; j<n; j++)
            model.removeRow(0);
        float[] range = anomaly(gamma1, gamma2);
        for(int i=0; i<data.size(); i++){
            if(data.getNumber(i)<range[0]){
                model.addRow(new String[] {""+data.getNumber(i)});
            }
            if(data.getNumber(i)>range[1]){
                model.addRow(new String[] {""+data.getNumber(i)});
            }
        }
    }
    
    public JFreeChart createHistogram(){
        HistogramDataset tmp = new HistogramDataset();
        DefaultTableXYDataset ty = new DefaultTableXYDataset();
        XYSeries sd = new XYSeries("dsv");
        tmp.setType(HistogramType.RELATIVE_FREQUENCY);
        double[] db = data.toArrayDouble(0);
        tmp.addSeries("", db, M, data.getNumber(0), data.getNumber(data.size()-1));
        DefaultTableXYDataset dataset = new DefaultTableXYDataset();
        String plotTitle = "Histogram with classes "+M;
        String xaxis = "X";
        String yaxis = "f(x)";
        PlotOrientation orientation = PlotOrientation.VERTICAL;
        boolean show = true;
        boolean toolTips = true;
        boolean urls = false;
        return  ChartFactory.createHistogram(plotTitle, xaxis, yaxis,
                tmp, orientation, show, toolTips, urls);
    }
    
    public int getm(){
        return M;
    }

    public int getM() {
        if(M==0 && data.getCountOfNumbers()<100){
            M=(int) Math.sqrt(data.getCountOfNumbers());
        }
        if(M==0 && data.getCountOfNumbers()>=100){
            M=(int) Math.pow(data.getCountOfNumbers(),1.0/3.0);
        }
        double tmp = (double)data.getNumber(data.size()-1)-(double)data.getNumber(0);
        step = (double)(tmp)/(double)M;
        process();
        return M;
    }
    
    public OrderedSeries getData(){
        return classes;
    }
    
    private void process(){
        int t=0;
        int count=0;
        double emprDis = 0.0;
        for(int i=1; i<=M; i++){
            count=0;
            double tmp = data.getNumber(0)+(double)i*step;
            for(int j=t; j<data.size();j++){
                if(data.getNumber(j)<=tmp){
                    count += data.getCount(j);
                }
                else{
                    t=j;
                    break;
                }
            }
            double relFreg = (double)count/(double)data.getCountOfNumbers();
            emprDis += relFreg;
            classes.addCount((float) tmp, count, relFreg, emprDis);
        }
    }

    public void setM(int M) {
        this.M = M;
        double tmp = (double)data.getNumber(data.size()-1)-(double)data.getNumber(0);
        step = (double)(tmp)/(double)M;
        classes = new OrderedSeries();
        process();
    }
    
    public Classes(OrderedSeries data){
        this.data = data;
        classes = new OrderedSeries();
    }
}
