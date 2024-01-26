namespace FEM.Solvers
{
    internal abstract class IterativeSolver : Solver
    {
        protected int maximumIterationNumber;
        protected int iterationtsNumber = 1;

        public double Accuracy { get; init; }
        protected IterativeSolver() : base()
        {
            try
            {
                using (StreamReader sr = new StreamReader("SLAE\\kuslau.txt"))
                {
                    matrixDimension = Int32.Parse(sr.ReadLine());
                    Accuracy = Double.Parse(sr.ReadLine());
                    maximumIterationNumber = Int32.Parse(sr.ReadLine());
                    sr.Close();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("The file kuslau.txt could not be read:");
                Console.WriteLine(e.Message);
            }
        }

        protected List<double> GetQStart()
        {
            List<double> q0 = new List<double>(new double[matrixDimension]);
            for (int i = 0; i < matrixDimension; i++)
            {
                q0[i] = 0.1;
            }
            return q0;
        }
    }
}
