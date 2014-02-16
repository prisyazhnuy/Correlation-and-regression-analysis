package correlation.and.regression.analysis;

public class Quantiles {
    //Сначала реализуем дополнительную функцию ошибок, которая понадобится при расчете обратной функции распределения стандартного нормального распределения.
        public static double Errorfunction(double x)
        {
            double Result = 0;
            double a = 0;

            a = 8 * (Math.PI - 3) / (3 * Math.PI * (4 - Math.PI));
            Result = x / Math.abs(x) * Math.sqrt(1 - Math.exp(-x * x * (4 / Math.PI + a * x * x) / (1 + a * x * x)));

            return 1 - Result;
        }
        //Теперь реализуем обратную функцию распределения стандартного нормального распределения.
        public static double InverseNormalCDF(double p)
        {
            double a1 = -3.969683028665376e+01;
            double a2 = 2.209460984245205e+02;
            double a3 = -2.759285104469687e+02;
            double a4 = 1.383577518672690e+02;
            double a5 = -3.066479806614716e+01;
            double a6 = 2.506628277459239e+00;
            double b1 = -5.447609879822406e+01;
            double b2 = 1.615858368580409e+02;
            double b3 = -1.556989798598866e+02;
            double b4 = 6.680131188771972e+01;
            double b5 = -1.328068155288572e+01;
            double c1 = -7.784894002430293e-03;
            double c2 = -3.223964580411365e-01;
            double c3 = -2.400758277161838e+00;
            double c4 = -2.549732539343734e+00;
            double c5 = 4.374664141464968e+00;
            double c6 = 2.938163982698783e+00;
            double d1 = 7.784695709041462e-03;
            double d2 = 3.224671290700398e-01;
            double d3 = 2.445134137142996e+00;
            double d4 = 3.754408661907416e+00;
            double p_low = 0.02425;
            double p_high = 1.0 - p_low;
            double q = 0, x = 0;

            if (0.0 < p && p < p_low)
            {
                q = Math.sqrt(-2.0 * Math.log(p));
                x = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
            }
            else if (p_low <= p && p <= p_high)
            {
                double r;
                q = p - 0.5;
                r = q * q;
                x = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q / (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0);
            }
            else
            {
                q = Math.sqrt(-2.0 * Math.log(1.0 - p));
                x = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
            }
            if (0.0 < p && p < 1.0)
            {
                double u, e;
                e = 0.5 * Errorfunction(-x / Math.sqrt(2.0)) - p;
                u = e * Math.sqrt(2.0 * Math.PI) * Math.exp(x * x / 2.0);
                x = x - u / (1.0 + x * u / 2.0);
            }

            return x;
        }
        //Осталось дело за квантилью распределения Стьюдента.
        public static double GetStudentsQuantile(double p, int n)
        {
            double a = 0, b = 0, c = 0, d = 0, x = 0, y = 0;

            a = 1.0 / (n - 0.5);
            b = 48.0 / (a * a);
            c = ((20700.0 * a / b - 98.0) * a - 16.0) * a + 96.36;
            d = ((94.5 / (b + c) - 3.0) / b + 1.0) * Math.sqrt(a * Math.PI / 2.0) * (double)n;
            x = d * p;
            y = Math.pow(x, 2.0 / (double)n);
            if (y > 0.05 + a)
            {
                x = InverseNormalCDF(p * 0.5);
                y = x * x;
                if (n < 5)
                    c += 0.3 * ((double)n - 4.5) * (x + 0.6);
                c += (((0.05 * d * x - 5) * x - 7.0) * x - 2.0) * x + b;
                y = (((((0.4 * y + 6.3) * y + 36.0) * y + 94.5) / c - y - 3.0) / b + 1.0) * x;
                y = a * y * y;
                if (y > 0.002)
                    y = Math.exp(y) - 1.0;
                else
                    y += 0.5 * y * y;

            }
            else
                y = ((1.0 / ((((double)n + 6.0) / ((double)n * y) - 0.089 * d - 0.822) * ((double)n + 2.0) * 3.0) + 0.5 / ((double)n + 4.0)) * y - 1.0) * ((double)n + 1.0) / ((double)n + 2.0) + 1.0 / y;

            return Math.sqrt((double)n * y);
        }
        
        public static double compute(double p, double mu, double sigma) {
		if(p < 0 || p > 1) 
			throw new RuntimeException("The probality p must be bigger than 0 and smaller than 1");		
		if(sigma < 0)
			throw new RuntimeException("The standard deviation sigma must be positive");
		if(p == 0)
			return Double.NEGATIVE_INFINITY;		
		if(p == 1)
			return Double.POSITIVE_INFINITY;		
		if(sigma == 0)
			return mu;		
		double  q, r, val;
 
		q = p - 0.5;
 
		/* 0.075 <= p <= 0.925 */
		if(Math.abs(q) <= .425) {
			r = .180625 - q * q;
			val = q * (((((((r * 2509.0809287301226727 +
		                    33430.575583588128105) * r + 67265.770927008700853) * r +
		                  45921.953931549871457) * r + 13731.693765509461125) * r +
		                1971.5909503065514427) * r + 133.14166789178437745) * r +
		              3.387132872796366608)
		         / (((((((r * 5226.495278852854561 +
		                  28729.085735721942674) * r + 39307.89580009271061) * r +
		                21213.794301586595867) * r + 5394.1960214247511077) * r +
		              687.1870074920579083) * r + 42.313330701600911252) * r + 1);
		}
		/* closer than 0.075 from {0,1} boundary */
		else {
	      /* r = min(p, 1-p) < 0.075 */
	      if (q > 0)
	          r = 1 - p;
	      else
	          r = p;
		      r = Math.sqrt(-Math.log(r));
		      /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */
 
		      if (r <= 5) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
		    	  r += -1.6;
		          val = (((((((r * 7.7454501427834140764e-4 +
		                     .0227238449892691845833) * r + .24178072517745061177) *
		                   r + 1.27045825245236838258) * r +
		                  3.64784832476320460504) * r + 5.7694972214606914055) *
		                r + 4.6303378461565452959) * r +
		               1.42343711074968357734)
		              / (((((((r *
		                       1.05075007164441684324e-9 + 5.475938084995344946e-4) *
		                      r + .0151986665636164571966) * r +
		                     .14810397642748007459) * r + .68976733498510000455) *
		                   r + 1.6763848301838038494) * r +
		                  2.05319162663775882187) * r + 1);
		      } else { /* very close to  0 or 1 */
		          r += -5;
		          val = (((((((r * 2.01033439929228813265e-7 +
		                     2.71155556874348757815e-5) * r +
		                    .0012426609473880784386) * r + .026532189526576123093) *
		                  r + .29656057182850489123) * r +
		                 1.7848265399172913358) * r + 5.4637849111641143699) *
		               r + 6.6579046435011037772)
		              / (((((((r *
		                       2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
		                      r + 1.8463183175100546818e-5) * r +
		                     7.868691311456132591e-4) * r + .0148753612908506148525)
		                   * r + .13692988092273580531) * r +
		                  .59983220655588793769) * r + 1);
		      }
 
		      if (q < 0.0) {
		          val = -val;
		      }
		  }
 
		  return mu + sigma * val;		
	}
}
