package correlation.and.regression.analysis;
    public class Quantiles
    {
        static double Phi(double alpha)
        {
            double c0 = 2.515517;
            double c1 = .802853;
            double c2 = .010328;
            double t = Math.sqrt(-2 * Math.log(alpha));
            double d1 = 1.432788;
            double d2 = .1892659;
            double d3 = .001308;
            return t - (c0 + c1 * t + c2 * t * t) / (1 + d1 * t + d2 * t * t + d3 * t * t * t);
        }
        public static double Norm(double p)
        {
            if (p <= .5)
            {
                return -1 * Phi(p);
            }
            else
            {
                return Phi(1 - p);
            }
        }
        public static double Student(double p, int nu)
        {
            double norm = Norm(p);
            double g1 = .25 * (Math.pow(norm, 3) + norm);
            double g2 = (5 * Math.pow(norm, 5) + 16 * Math.pow(norm, 3) + 3 * norm) / 96;
            double g3 = (3 * Math.pow(norm, 7) + 19 * Math.pow(norm, 5) + 17 * Math.pow(norm, 3) - 15 * norm)/384;
            double g4 = (79 * Math.pow(norm, 9) + 779 * Math.pow(norm, 7) + 1482 * Math.pow(norm, 5) - 1920 * Math.pow(norm, 3) - 945 * norm) / 92160;
            return norm + g1 / nu + g2 / Math.pow(nu, 2) + g3 / Math.pow(nu, 3) + g4 / Math.pow(nu, 4);
        }
        
        public static double Fisher(double p, int nu1, int nu2){
            System.out.println("nu1"+nu1);
            System.out.println("nu2"+nu2);
            double g = G(nu1, nu2);
            System.out.println("g"+g);
            double j = J(nu1, nu2);
            System.out.println("j"+j);
            double u = Norm(p);
            System.out.println("u "+u);
            double a = u*Math.sqrt(g/2)-(j*Math.pow(u, 2)+2*j)/6;
            double b = Math.sqrt(g/2)*((g*u*u+g*3*u)/24+j*j*(u*u*u+11*u)/72*g);
            double c = g*j*(u*u*u*u+9*u*u+8)/120;
            double d = j*j*j*(3*u*u*u*u+7*u*u-16)/(3240*g);
            double e = Math.sqrt(g/2)*((g*g*(u*u*u*u*u+20*u*u*u+15*u)/1920)+(j*j*j*j*(u*u*u*u*u+44*u*u*u+18*u)/2880)+(j*j*j*j*(9*u*u*u*u*u-284*u*u*u-1513-u)/(155520*g*g)));
            return Math.pow(Math.E, 2*(a+b-c+d+e));
        }
        
        private static double G(int nu1, int nu2){
            return (1/(double)nu1)+(1/(double)nu2);
        }
        
        private static double J(int nu1, int nu2){
            return (1/(double)nu1)-(1/(double)nu2);
        }
    }