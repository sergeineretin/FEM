
using System;
using System.Runtime.Serialization.Formatters;

namespace FEM.Solvers
{
    internal class CGMIncompleteLLTSolver : CGMSolver
    {
        List<double> gglCholesky = new List<double>();
        List<double> diCholesky = new List<double>();
        public CGMIncompleteLLTSolver() : base() 
        {
            gglCholesky.AddRange(ggl);
            diCholesky.AddRange(di);
            ShowProgress = true;
        }

        public override void Solve()
        {
            double alphaK = 0;
            double betaK = 0;
            double bEucl = EuclideanNorm(b);
            GenerateCholeskyMatrix();
            x = GetQStart();

            List<double> Ax0 = matrixMultiplicationByVector(x);
            for (int i = 0; i < matrixDimension; i++)
                rk[i] = b[i] - Ax0[i];

            List<double> inverseMrkPrev = Cholesky(rk);

            for (int i = 0; i < matrixDimension; i++)
                zk[i] = inverseMrkPrev[i];

            double rrk = ScalarProduct(zk, rk);

            double accuracy = 100000;
            while (Accuracy < accuracy && maximumIterationNumber > iterationtsNumber++)
            {
                List<double> rkPrev = new List<double>(new double[matrixDimension]);
                for (int i = 0; i < matrixDimension; i++)
                    rkPrev[i] = rk[i];

                List<double> Azk = matrixMultiplicationByVector(zk);

                double rrk1 = rrk;
                alphaK = rrk1 / ScalarProduct(Azk, zk);

                for (int i = 0; i < matrixDimension; i++)
                {
                    x[i] = x[i] + alphaK * zk[i];
                    rk[i] = rk[i] - alphaK * Azk[i];
                }

                List<double> inverseMrk = Cholesky(rk);
                rrk = ScalarProduct(inverseMrk, rk);
                betaK = rrk / rrk1;

                for (int i = 0; i < matrixDimension; i++) {
                    zk[i] = inverseMrk[i] + betaK * zk[i];
                }
                accuracy = EuclideanNorm(rk) / bEucl;
            }
            if (ShowProgress)
            {
                Console.WriteLine($"Number of iterations {iterationtsNumber}.");
                Console.WriteLine($"Accuracy is {accuracy}.");
            }
        }

        private void GenerateCholeskyMatrix()
        {
            double sum_d, sum_l;
            for (int k = 0; k < matrixDimension; k++)
            {
                sum_d = 0;
                int i_s = ig[k], i_e = ig[k + 1];
                for (int i = i_s; i < i_e; i++)
                {
                    sum_l = 0;
                    int j_s = ig[jg[i]], j_e = ig[jg[i] + 1];

                    for (int m = i_s; m < i; m++)
                    {
                        for (int j = j_s; j < j_e; j++)
                        {
                            if (jg[m] == jg[j])
                            {
                                sum_l += gglCholesky[m] * gglCholesky[j];
                                j_s++;
                            }
                        }
                    }
                    gglCholesky[i] = (gglCholesky[i] - sum_l) / diCholesky[jg[i]];

                    sum_d += gglCholesky[i] * gglCholesky[i];
                }
                diCholesky[k] = Math.Sqrt(diCholesky[k] - sum_d);
            }
        }

        private double SumIJ(int i, int j)
        {
            double sum = 0;
            for (int jj = ig[j]; jj < ig[j + 1]; jj++)
            {
                for (int ii = ig[i]; ii < ig[i + 1] && jg[ii] <= jg[jj]; ii++)
                {
                    if (jg[ii] == jg[jj])
                    {
                        sum += gglCholesky[ii] * gglCholesky[jj];
                        break;
                    }
                }
            }
            return sum;
        }

        private List<double> Cholesky(List<double> rk)
        {
            return Reverse(Forward(rk));
        }

        private List<double> Forward(List<double> b)
        { 
            List<double> y = new List<double>(new double[matrixDimension]);
            for (int i = 0; i < matrixDimension; i++)
            {
                double sum = 0;
                for (int j = ig[i]; j < ig[i + 1]; j++)
                    sum += y[jg[j]] * gglCholesky[j];

                y[i] = (b[i] - sum) / diCholesky[i];
            }
            return y;
        }
        private List<double> Reverse(List<double> y) 
        {
            for (int ii = matrixDimension - 1; ii >= 0; ii--) {
                double yii = y[ii] / diCholesky[ii];
                y[ii] = yii;
                for (int jj = ig[ii]; jj < ig[ii+1]; jj++)
                {
                    int nst = jg[jj];
                    y[nst] -= gglCholesky[jj] * yii;
                }
            }
            return y;
        }
    }
}
