namespace FEM.Solvers
{
    internal sealed class LOSSolver : IterativeSolver
    {
        public LOSSolver()
            : base()
        { }
        public override void Solve()
        {
            List<double> rk = new List<double>(new double[matrixDimension]);
            List<double> zk = new List<double>(new double[matrixDimension]);
            List<double> pk = new List<double>(new double[matrixDimension]);

            double alphaK = 0;
            double betaK = 0;
            x = GetQStart();

            List<double> Ax0 = matrixMultiplicationByVector(x);
            for (int i = 0; i < matrixDimension; i++)
                rk[i] = zk[i] = b[i] - Ax0[i];

            List<double> Ar0 = matrixMultiplicationByVector(rk);
            for (int i = 0; i < matrixDimension;i++)
                pk[i] = Ar0[i];
            double r0r0 = ScalarProduct(rk, rk);
            double rkrk = r0r0;
            int iterationtsNumber = 0;

            while (Accuracy * Accuracy < rkrk / r0r0 && maximumIterationNumber > iterationtsNumber++)
            {
                double pkpk = ScalarProduct(pk, pk);
                alphaK = ScalarProduct(pk, rk) / pkpk;
                for (int i = 0; i < matrixDimension; i++) {
                    x[i] = x[i] + alphaK * zk[i];
                    rk[i] = rk[i] - alphaK * pk[i];
                }
                List<double> Ark = matrixMultiplicationByVector(rk);
                betaK = -ScalarProduct(pk, Ark) / pkpk;
                for (int i = 0; i < matrixDimension; i++)
                { 
                    zk[i] = rk[i] + betaK * zk[i];
                    pk[i] = Ark[i] + betaK * pk[i];
                }

                rkrk = ScalarProduct(rk, rk);
            }
            Console.WriteLine($"Number of iterations {iterationtsNumber}.");
            Console.WriteLine($"Accuracy is {Math.Sqrt(rkrk/r0r0)}.");
        }
    }
}
