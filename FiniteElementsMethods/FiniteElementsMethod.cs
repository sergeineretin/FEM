using FEM.Points;
using FEM.Elements;
using FEM.Problems;
using FEM.BoundaryСonditions;
using System.Text.Json.Serialization;
using System.Xml.Linq;

namespace FEM.FiniteElementsMethods
{
    public struct BC1
    {
        public BC1(int node, int value)
        {
            Node = node;
            Subarea = value;
        }
        public int Node { get; private set; }
        public int Subarea { get; private set; }

        public override string ToString() { return $"{Node} {Subarea}"; }
    }

    public struct BC2
    {
        int node0 = 0;
        int node1 = 0;
        int subarea = 0;
        public BC2(int node0, int node1, int subarea)
        {
            this.node0 = node0;
            this.node1 = node1;
            this.subarea = subarea;
        }
        public int Node0 { get { return node0; } }
        public int Node1 { get { return node1; } }
        public int Subarea { get { return subarea; } }
    }

    public struct BC3
    {
        int node0 = 0;
        int node1 = 0;
        int subarea = 0;
        public BC3(int node0, int node1, int subarea)
        {
            this.node0 = node0;
            this.node1 = node1;
            this.subarea = subarea;
        }
        public int Node0 { get { return node0; } }
        public int Node1 { get { return node1; } }

        public int Subaria { get { return subarea; } }
    }

    public abstract class FiniteElementsMethod
    {

        protected List<TriangularElement> elems = new List<TriangularElement>();
        protected List<Point2D> nodes = new List<Point2D>();
        protected List<int> igEdge = new List<int>();
        protected List<int> jgEdge = new List<int>();
        protected List<List<Point2D>> ggEdge = new List<List<Point2D>>();
        protected List<List<int>> unknownsOfEdges = new List<List<int>>();
        protected List<List<int>> unknownsOfElements = new List<List<int>>();
        protected List<List<int>> unknownsOfNodes = new List<List<int>>();



        protected List<int> ig = new List<int>();
        protected List<int> jg = new List<int>();
        protected List<double> ggl;
        protected List<double> ggu;
        protected List<double> di;
        protected List<double> b;

        protected Problem problem;

        public List<Point2D> Nodes { get { return nodes; } }
        public List<TriangularElement> Elems { get { return elems; } }
        public List<double> Q { get; set; }

        public int N;

        GridGenerator gridGenerator;

        protected abstract int GetNumberOfGlobalBasisFunctions();
        protected abstract void GetGrid();

        public void GenSLAE(Problem problem)
        {
            this.problem = problem;
            ClearSLAE();
            //GetGrid();
            BuildPortraitOfMatrix();
            AddElementsToGlobalSLAE();
            // AddConcentratedSources();
            AddBc2();
            //AddBc3();
            AddBc1();
            WriteSLAE();

        }

        public void PrintSolution()
        {
            double rr = 0;
            double r = 0;
            for (int i = 0; i < Nodes.Count; i++)
            {
                Console.WriteLine($"{nodes[i]}  {Q[unknownsOfNodes[i][0]]}  {problem.Bc1(nodes[i],1)}  {Math.Abs(Q[unknownsOfNodes[i][0]] - problem.Bc1(nodes[i], 1))}");
                rr += (Q[unknownsOfNodes[i][0]] - problem.Bc1(nodes[i], 1)) * (Q[unknownsOfNodes[i][0]] - problem.Bc1(nodes[i], 1));
                r += (problem.Bc1(nodes[i], 1)) * problem.Bc1(nodes[i], 1);
            }

            for (int i = 0; i < igEdge.Count - 1; i++)
            {
                for (int j = igEdge[i]; j < igEdge[i+1]; j++)
                {
                    for (int k = 0; k < unknownsOfEdges[j].Count; k++)
                    {
                        Console.WriteLine($"{ggEdge[j][k]}  {Q[unknownsOfEdges[j][k]]}  {problem.Bc1(ggEdge[j][k], 1)}  {Math.Abs(Q[unknownsOfEdges[j][k]] - problem.Bc1(ggEdge[j][k], 1))}");
                        rr += (Q[unknownsOfEdges[j][k]] - problem.Bc1(ggEdge[j][k], 1)) * (Q[unknownsOfEdges[j][k]] - problem.Bc1(ggEdge[j][k], 1));
                        r += problem.Bc1(ggEdge[j][k], 1) * problem.Bc1(ggEdge[j][k], 1);
                    }
                }
            }

            for (int i = 0; i < elems.Count; i++)
            {
                Point2D p = new Point2D(0, 0);
                p.X = (nodes[elems[i].Node0].X + nodes[elems[i].Node1].X + nodes[elems[i].Node2].X) / 3;
                p.Y = (nodes[elems[i].Node0].Y + nodes[elems[i].Node1].Y + nodes[elems[i].Node2].Y) / 3;
                Console.WriteLine($"{p} {Q[unknownsOfElements[i][0]]} {problem.Bc1(p, elems[i].Subarea)} {Math.Abs(Q[unknownsOfElements[i][0]] - problem.Bc1(p, elems[i].Subarea))}");
                rr += (Q[unknownsOfElements[i][0]] - problem.Bc1(p, elems[i].Subarea)) * (Q[unknownsOfElements[i][0]] - problem.Bc1(p, elems[i].Subarea));
                r += problem.Bc1(p, elems[i].Subarea) * problem.Bc1(p, elems[i].Subarea);
            }
            Console.WriteLine($"relative residual = {Math.Sqrt(rr/r)}");
        }

        public void GenerateGrid(string mediaGeometryFileName, string boundaryFileName)
        {
            gridGenerator = new GridGenerator(mediaGeometryFileName, boundaryFileName);
            gridGenerator.GenerateGrid();
            gridGenerator.WriteGrid();
            gridGenerator.DefineGlobalIndicesOfUnknownsOnElements();
            nodes = gridGenerator.Nodes;
            elems = gridGenerator.Elems;
            igEdge = gridGenerator.IgEdge;
            jgEdge = gridGenerator.JgEdge;
            ggEdge = gridGenerator.GgEdge;
            unknownsOfEdges = gridGenerator.UnknownsOfEdges;
            unknownsOfElements = gridGenerator.UnknownsOfElements;
            unknownsOfNodes = gridGenerator.UnknownsOfNodes;
            N = gridGenerator.NumberOfGlobalBasisFunctions;
        }
        private void ClearSLAE()
        {
            if (b != null)
                b.Clear();
            if (ggl!= null)
                ggl.Clear(); 
            if (ggu != null)
                ggu.Clear();
            if (di != null)
                di.Clear();
            if (ig != null)
                ig.Clear();
            if (jg != null)
                jg.Clear();
        }

        protected abstract void AddConcentratedSources();
        public void WriteKuslau(double accuracy, int maximumIterationsNumber)
        {
            try
            {
                File.WriteAllText("SLAE\\kuslau.txt", string.Empty);
                using (StreamWriter sw = new StreamWriter("SLAE\\kuslau.txt", false))
                {
                    sw.WriteLine(ig.Count - 1);
                    sw.WriteLine(accuracy);
                    sw.WriteLine(maximumIterationsNumber);
                    sw.Close();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
        }
        public abstract double GetSolution(Point2D point, TriangularElement element);
        private void WriteSLAE()
        {
            WriteListToFile("SLAE\\ig.txt", ig);
            WriteListToFile("SLAE\\jg.txt", jg);
            WriteListToFile("SLAE\\ggl.txt", ggl);
            WriteListToFile("SLAE\\ggu.txt", ggu);
            WriteListToFile("SLAE\\di.txt", di);
            WriteListToFile("SLAE\\pr.txt", b);
        }
        private void WriteListToFile(string fileName, List<int> list)
        {
            try
            {
                File.WriteAllText(fileName, string.Empty);
                using (StreamWriter sw = new StreamWriter(fileName, false))
                {
                    foreach (var item in list)
                        sw.WriteLine(item.ToString());
                    sw.Close();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
        }
        private void WriteListToFile(string fileName, List<double> list)
        {
            try
            {
                File.WriteAllText(fileName, string.Empty);
                using (StreamWriter sw = new StreamWriter(fileName, false))
                {
                    foreach (var item in list)
                        sw.WriteLine(item.ToString());
                    sw.Close();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
        }

        protected abstract void AddElementsToGlobalSLAE();

        protected abstract void AddLocalStiffnessMatrix(Element element);

        protected abstract void AddLocalMassMatrix(Element element);

        protected abstract void AddToVector(Element element);



        protected void AddElementToGlobal(int i, int j, double a)
        {
            if (i == j)
                di[i] += a;
            else if (i > j)
            {
                for (int ii = ig[i]; ii < ig[i + 1]; ii++)
                    if (jg[ii] == j) { ggl[ii] += a; break; }
            }
            else
            {
                for (int jj = ig[j]; jj < ig[j + 1]; jj++)
                    if (jg[jj] == i) { ggu[jj] += a; break; }
            }
        }

        protected void AddBc1()
        {
            List<DirichletBoundaryCondition> bc1 = DataReader.ReadBC1FromFile();

            List<int> consideredIndexes = new List<int>();

            for (int i = 0; i < bc1.Count; i++)
            {
                if (!consideredIndexes.Contains(unknownsOfNodes[bc1[i].NodeNumber1][0]))
                {
                    AddDirichletPoint(unknownsOfNodes[bc1[i].NodeNumber1][0], nodes[bc1[i].NodeNumber1], bc1[i].FormulaNumber);
                    consideredIndexes.Add(unknownsOfNodes[bc1[i].NodeNumber1][0]);
                }
                if (!consideredIndexes.Contains(unknownsOfNodes[bc1[i].NodeNumber2][0]))
                {
                    AddDirichletPoint(unknownsOfNodes[bc1[i].NodeNumber2][0], nodes[bc1[i].NodeNumber2], bc1[i].FormulaNumber);
                    consideredIndexes.Add(unknownsOfNodes[bc1[i].NodeNumber2][0]);
                }

                List<int> edgeUnknowns = GetUnknownsOfEdge(bc1[i].NodeNumber1, bc1[i].NodeNumber2);
                List<Point2D> edgePoints = GetPointsOfEdge(bc1[i].NodeNumber1, bc1[i].NodeNumber2);
                for (int j = 0; j < edgeUnknowns.Count; j++)
                {
                    AddDirichletPoint(edgeUnknowns[j], edgePoints[j], bc1[i].FormulaNumber);
                }

            }
        }
        protected List<int> GetUnknownsOfEdge(int iNode, int jNode)
        {
            if (iNode > jNode)
            {
                int k = iNode;
                iNode = jNode;
                jNode = k;
            }

            for (int jj = igEdge[iNode]; jj < igEdge[iNode + 1]; jj++)
            {
                if (jgEdge[jj] == jNode)
                {
                    return unknownsOfEdges[jj];
                }
            }
            return null;
        }

        protected List<Point2D> GetPointsOfEdge(int iNode, int jNode)
        {
            if (iNode > jNode)
            {
                int k = iNode;
                iNode = jNode;
                jNode = k;
            }

            for (int jj = igEdge[iNode]; jj < igEdge[iNode + 1]; jj++)
            {
                if (jgEdge[jj] == jNode)
                {
                    return ggEdge[jj];
                }
            }
            return null;
        }
        protected void AddDirichletPoint(int unknownIndex, Point2D point, int formulaNumber)
        {
            di[unknownIndex] = 1;
            b[unknownIndex] = problem.Bc1(point, formulaNumber);
            for (int ii = ig[unknownIndex]; ii < ig[unknownIndex + 1]; ii++)
            {
                ggl[ii] = 0;

                b[jg[ii]] -= ggu[ii] * problem.Bc1(point, formulaNumber);
                ggu[ii] = 0;
            }
            for (int ii = unknownIndex + 1; ii < N; ii++)
            {
                for (int jj = ig[ii]; jj < ig[ii + 1]; jj++)
                {
                    if (jg[jj] == unknownIndex)
                    {
                        ggu[jj] = 0;

                        b[ii] -= ggl[jj] * problem.Bc1(point, formulaNumber);
                        ggl[jj] = 0;
                        break;
                    }
                }
            }
        }

        protected abstract void AddBc2();
        protected abstract void AddBc3();
        protected double GetLenght(int inode1, int inode2)
        {
            return Math.Sqrt((((Point2D)nodes[inode1]).X - ((Point2D)nodes[inode2]).X) * (((Point2D)nodes[inode1]).X - ((Point2D)nodes[inode2]).X) +
                    (((Point2D)nodes[inode1]).Y - ((Point2D)nodes[inode2]).Y) * (((Point2D)nodes[inode1]).Y - ((Point2D)nodes[inode2]).Y));
        }
        protected void BuildPortraitOfMatrix()
        {
            List<List<int>> list = new List<List<int>>();

            for (int inodes = 0; inodes < N; inodes++)
                list.Add(new List<int>());

            for (int ielem = 0; ielem < elems.Count; ielem++)
            {
                for (int i = 0; i < (elems[ielem]).GetNumberOfUnknowns(); i++)
                {
                    for (int j = 0; j < (elems[ielem]).GetNumberOfUnknowns(); j++)
                    {
                        if ((elems[ielem]).GetIndexOfUnknown(i) > (elems[ielem]).GetIndexOfUnknown(j) &&
                            !list[(elems[ielem]).GetIndexOfUnknown(i)].Contains((elems[ielem]).GetIndexOfUnknown(j)))
                            list[(elems[ielem]).GetIndexOfUnknown(i)].Add((elems[ielem]).GetIndexOfUnknown(j));
                    }
                }
            }

            int curSum = 0;
            ig.Add(curSum);
            for (int inodes = 0; inodes < N; inodes++)
            {
                list[inodes].Sort();
                for (int i = 0; i < list[inodes].Count; i++)
                {
                    jg.Add(list[inodes][i]);
                    curSum++;
                }
                ig.Add(curSum);
            }

            di = new List<double>(new double[ig.Count - 1]);
            b = new List<double>(new double[ig.Count - 1]);
            ggl = new List<double>(new double[jg.Count]);
            ggu = new List<double>(new double[jg.Count]);
        }
        public void MatrixRenderer()
        {
            for (int i = 0; i < ig.Count - 1; i++)
            {
                string line = "";
                for (int j = 0; j < ig.Count - 1; j++)
                {
                    bool finded = false;

                    if (i == j)
                        line += $"{di[i]:N5}";
                    else if (i > j)
                    {
                        for (int ii = ig[i]; ii < ig[i + 1]; ii++)
                            if (jg[ii] == j) { line += $"{ggl[ii]:N5}"; finded = true; break; }

                    }
                    else
                    {
                        for (int jj = ig[j]; jj < ig[j + 1]; jj++)
                            if (jg[jj] == i) { line += $"{ggu[jj]:N5}"; finded = true; break; }
                    }
                    if (!finded) { line += "0"; }
                    line += " ";
                }
                Console.WriteLine(line);
            }
        }
    }
}
