
namespace FEM.Solvers
{
    internal abstract class Solver
    {
        protected int matrixDimension;

        protected List<int> ig = new List<int>();
        protected List<int> jg = new List<int>();
        protected List<double> ggl = new List<double>();
        protected List<double> ggu = new List<double>();
        protected List<double> di = new List<double>();
        protected List<double> b = new List<double>();
        protected List<double> x;
        public bool ShowProgress { get; set; }

        public List<double> Q { get { return x; } }
        public Solver() 
        {
            ReadListFromFile("SLAE\\ig.txt", ig);
            matrixDimension = ig.Count - 1;
            x = new List<double>(new double[matrixDimension]);
            ReadListFromFile("SLAE\\jg.txt", jg);
            ReadListFromFile("SLAE\\ggl.txt", ggl);
            ReadListFromFile("SLAE\\ggu.txt", ggu);
            ReadListFromFile("SLAE\\di.txt", di);
            ReadListFromFile("SLAE\\pr.txt", b);
            ggl = ggu;
        }

        private void ReadListFromFile(string fileName, List<int> list)
        {
            try
            {
                using (StreamReader sr = new StreamReader(fileName))
                {
                    string line;
                    while ((line = sr.ReadLine()) != null)
                    {
                        list.Add(Int32.Parse(line));
                    }
                    sr.Close();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("The file " + fileName + " could not be read:");
                Console.WriteLine(e.Message);
            }
        }
        private void ReadListFromFile(string fileName, List<double> list)
        {
            try
            {
                using (StreamReader sr = new StreamReader(fileName))
                {
                    string line;
                    while ((line = sr.ReadLine()) != null)
                    {
                        list.Add(Double.Parse(line));
                    }
                    sr.Close();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("The file " + fileName + " could not be read:");
                Console.WriteLine(e.Message);
            }
        }


        public abstract void Solve();


        public void WriteSolution(string filePath)
        {
            if (Q.Count != 0)
            {
                try
                {
                    File.WriteAllText(filePath, string.Empty);
                    using (StreamWriter sw = new StreamWriter(filePath, false))
                    {
                        sw.WriteLine(Q.Count);
                        foreach (var q in Q)
                            sw.WriteLine(q.ToString());
                        sw.Close();
                    }
                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                }
            }
        }
        protected double ScalarProduct(List<double> v1, List<double> v2)
        {
            double result = 0;
            if (v1.Count == v2.Count)
            {
                for (int i = 0; i < v1.Count; i++)
                {
                    result += v1[i] * v2[i];
                }
                return result;
            }
            else
                return 0;
        }

        protected List<double> matrixMultiplicationByVector(List<double> vector)
        { 
            List<double> result = new List<double>(new double[matrixDimension]);
            for (int i = 0; i < matrixDimension; i++)
            {
                result[i] += di[i] * vector[i];
                for (int j = ig[i]; j < ig[i + 1]; j++)
                {
                    result[i] += ggl[j] * vector[jg[j]];
                    result[jg[j]] += ggu[j] * vector[i];
                }
            }
            return result;
        }

        protected double EuclideanNorm(List<double> vector)
        {
            double norm = 0;
            for (int i = 0; i < vector.Count; i++)
            {
                norm += vector[i] * vector[i];
            }
            return Math.Sqrt(norm);
        }

        public List<double> GetSolution()
        { 
            return x;
        }

        public void PrintSolutionToConsole()
        {
            String line = "";
            for (int i = 0; i < x.Count; i++) { line += $"{x[i]:N3}\numberOfReseivers"; }
            Console.WriteLine(line);
        }
    }
}
