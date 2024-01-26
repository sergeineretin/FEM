using FEM.Points;
using FEM.Elements;
using FEM.FiniteElementsMethods;
using FEM.BoundaryСonditions;

namespace FEM
{
    public static class DataReader
    {
        public static List<TriangularElement> ReadElemsFromFile()
        {
            List<TriangularElement> elems = new List<TriangularElement>();
            try
            {
                using (StreamReader elemSR = new StreamReader("elems.txt"))
                {
                    int n = Int32.Parse(elemSR.ReadLine());
                    string line;
                    while ((line = elemSR.ReadLine()) != null)
                    {
                        List<int> ints = new List<int>(line.Trim().Split(" ").ToList().Select(e => int.Parse(e)));
                        elems.Add(new TriangularElement(ints[0], ints[1], ints[2], ints[3]));
                    }
                    elemSR.Close();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("The file elems.txt could not be read:");
                Console.WriteLine(e.Message);
            }
            return elems;
        }
        public static List<Point2D> ReadNodesFromFile()
        {
            List<Point2D> nodes = new List<Point2D>();
            try
            {
                using (StreamReader nodeSR = new StreamReader("nodes.txt"))
                {
                    int n = Int32.Parse(nodeSR.ReadLine());
                    string line;
                    while ((line = nodeSR.ReadLine()) != null)
                    {
                        List<double> coords = new List<double>(line.Trim().Split(" ").ToList().Select(e => double.Parse(e)));
                        nodes.Add(new Point2D(coords[0], coords[1]));
                    }
                    nodeSR.Close();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("The file nodes.txt could not be read:");
                Console.WriteLine(e.Message);
            }
            return nodes;
        }

        public static List<DirichletBoundaryCondition> ReadBC1FromFile()
        {
            List<DirichletBoundaryCondition> bc1 = new List<DirichletBoundaryCondition>();
            try
            {
                using (StreamReader bc1SR = new StreamReader("bc1.txt"))
                {
                    int n = Int32.Parse(bc1SR.ReadLine());
                    if (n != 0)
                    {
                        string line;
                        while ((line = bc1SR.ReadLine()) != null)
                        {
                            List<double> bc = new List<double>(line.Replace('\t', ' ').Trim().Split(" ").ToList().Select(e => double.Parse(e)));
                            bc1.Add(new DirichletBoundaryCondition((int)bc[0], (int)bc[1], (int)bc[2]));
                        }
                    }
                    bc1SR.Close();

                }
            }
            catch (Exception e)
            {
                Console.WriteLine("The file bc1.txt could not be read:");
                Console.WriteLine(e.Message);
            }
            return bc1;
        }

        public static List<NeumannBoundaryCondition> ReadBC2FromFile()
        {
            List<NeumannBoundaryCondition> bc2 = new List<NeumannBoundaryCondition>();
            try
            {
                using (StreamReader bc2SR = new StreamReader("bc2.txt"))
                {
                    int n = Int32.Parse(bc2SR.ReadLine());
                    if (n != 0)
                    {
                        string line;
                        while ((line = bc2SR.ReadLine()) != null)
                        {
                            List<double> bc = new List<double>(line.Replace('\t', ' ').Split(" ").ToList().Select(e => double.Parse(e)));
                            bc2.Add(new NeumannBoundaryCondition((int)bc[0], (int)bc[1], (int)bc[2]));
                        }
                    }
                    bc2SR.Close();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("The file bc2.txt could not be read:");
                Console.WriteLine(e.Message);
            }
            return bc2;
        }
        static public List<Point2D> ReadPointSourcesFromFile()
        {
            List<Point2D> sources = new List<Point2D>();
            try
            {
                using (StreamReader sr = new StreamReader("pointSources.txt"))
                {
                    int n = Int32.Parse(sr.ReadLine());
                    if (n != 0)
                    {
                        string line;
                        while ((line = sr.ReadLine()) != null)
                        {
                            List<double> x = new List<double>(line.Replace('\t', ' ').Split(" ").ToList().Select(e => double.Parse(e)));
                            sources.Add(new(x[0], x[1]));
                        }
                    }
                    sr.Close();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("The file pointSources.txt could not be read:");
                Console.WriteLine(e.Message);
            }
            return sources;
        }
        public static List<BC3> ReadBC3FromFile()
        {
            List<BC3> bc3 = new List<BC3>();

            try
            {
                using (StreamReader bc3SR = new StreamReader("bc3.txt"))
                {
                    int n = Int32.Parse(bc3SR.ReadLine());
                    if (n != 0)
                    {
                        string line;
                        while ((line = bc3SR.ReadLine()) != null)
                        {
                            List<double> bc = new List<double>(line.Split(" ").ToList().Select((e) => double.Parse(e)));
                            bc3.Add(new BC3((int)bc[0], (int)bc[1], (int)bc[2]));
                        }
                    }
                    bc3SR.Close();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("The file bc3.txt could not be read:");
                Console.WriteLine(e.Message);
            }
            return bc3;
        }

    }
}
