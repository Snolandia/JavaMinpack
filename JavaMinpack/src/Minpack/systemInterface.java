package Minpack;

public class systemInterface {
		public interface func {
			double[] operation(double[] a);
		}
		public interface jacob {
			double[][] operation(double[] a);
		}
}
