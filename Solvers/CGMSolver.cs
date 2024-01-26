namespace FEM.Solvers
{
    internal class CGMSolver : IterativeSolver
    {
        protected List<double> rk;
        protected List<double> zk;
        public CGMSolver() : base()
        {
            rk = new List<double>(new double[matrixDimension]);
            zk = new List<double>(new double[matrixDimension]);
        }
        
        public override void Solve()
        {
            double alphaK = 0;
            double betaK = 0;
            x = GetQStart();

            List<double> Ax0 = matrixMultiplicationByVector(x);
            for (int i = 0; i < matrixDimension; i++)
                rk[i] = zk[i] = b[i] - Ax0[i];

            double rkrk = ScalarProduct(rk, rk);
            int iterationtsNumber = 0;

            double bEucl = EuclideanNorm(b);
            while (Accuracy < EuclideanNorm(rk) / bEucl && maximumIterationNumber > iterationtsNumber++)
            {
                List<double> Azk = matrixMultiplicationByVector(zk);
                double Azkpk = ScalarProduct(Azk, zk);
                alphaK = ScalarProduct(rk, rk) / Azkpk;

                for (int i = 0; i < matrixDimension; i++)
                {
                    x[i] = x[i] + alphaK * zk[i];
                    rk[i] = rk[i] - alphaK * Azk[i];
                }
                List<double> Ark = matrixMultiplicationByVector(rk);
                betaK = ScalarProduct(rk, rk) / rkrk;
                rkrk = ScalarProduct(rk, rk);

                for (int i = 0; i < matrixDimension; i++)
                    zk[i] = rk[i] + betaK * zk[i];

            }
            Console.WriteLine($"Number of iterations {iterationtsNumber}.");
            Console.WriteLine($"Accuracy is {EuclideanNorm(rk) / bEucl}.");
        }
    }
}
