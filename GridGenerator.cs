using FEM.Elements;
using FEM.FiniteElementsMethods;
using FEM.Points;
using FEM.BoundaryСonditions;

namespace FEM
{ 
    struct Subarea
    {
        public int SubariaNumber { get; private set; }
        public int X0 { get; private set; }
        public int X1 { get; private set; }
        public int Y0 { get; private set; }
        public int Y1 { get; private set; }
        public Subarea(int subariaNumber, int x0, int x1, int y0, int y1)
        {
            SubariaNumber = subariaNumber;
            X0 = x0;
            X1 = x1;
            Y0 = y0;
            Y1 = y1;
        }
        public override string ToString()
        {
            return $"{SubariaNumber}\t{X0}\t{X1}\t{Y0}\t{Y1}";
        }


    }
    struct Bound
    {
        public int BoundaryType { get; private set; }
        public int FormulaNumber { get; private set; }
        public int X0 { get; private set; }
        public int X1 { get; private set; }
        public int Y0 { get; private set; }
        public int Y1 { get; private set; }
        public Bound(int boundaryType, int formulaNumber, int x0, int x1, int y0, int y1)
        {
            BoundaryType = boundaryType;
            FormulaNumber = formulaNumber;
            X0 = x0;
            X1 = x1;
            Y0 = y0;
            Y1 = y1;
        }
        public override string ToString()
        {
            return $"{BoundaryType}\t{FormulaNumber}\t{X0}\t{X1}\t{Y0}\t{Y1}";
        }

        public override bool Equals(object? obj)
        {
            return obj is Bound bound &&
                   X0 == bound.X0 &&
                   X1 == bound.X1 &&
                   Y0 == bound.Y0 &&
                   Y1 == bound.Y1;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(X0, X1, Y0, Y1);
        }
    }
    struct Bc1
    {
        public int NodeNumber { get; set; }
        public int FormulaNumber { get; private set; }

        public Bc1(int nodeNumber, int formulaNumber)
        {
            NodeNumber = nodeNumber;
            FormulaNumber = formulaNumber;
        }
        public override string ToString()
        {
            return $"{NodeNumber}\t{FormulaNumber}";
        }
    }
    struct Bc2
    {
        public int NodeNumber1 { get; set; }
        public int NodeNumber2 { get; set; }
        public int FormulaNumber { get; set; }

        public Bc2(int nodeNumber1, int nodeNumber2, int formulaNumber)
        {
            NodeNumber1 = nodeNumber1;
            NodeNumber2 = nodeNumber2;
            FormulaNumber = formulaNumber;
        }
        public override string ToString()
        {
            return $"{NodeNumber1}\t{NodeNumber2}\t{FormulaNumber}";
        }
    }
    struct QuadrupleInteger
    { 
        public int Yi1 { get; init; }
        public int Xi1 { get; init; }
        public int Yi2 { get; init; }
        public int Xi2 { get; init; }
        public QuadrupleInteger(int yi1, int xi1, int yi2, int xi2)
        {
            Yi1 = yi1;
            Xi1 = xi1;
            Yi2 = yi2;
            Xi2 = xi2;
        }

        public override string ToString()
        {
            return $"{Yi1}\t{Xi1}\t{Yi2}\t{Xi2}";
        }

        public override bool Equals(object? obj)
        {
            return obj is QuadrupleInteger integer &&
                   Yi1 == integer.Yi1 &&
                   Xi1 == integer.Xi1 &&
                   Yi2 == integer.Yi2 &&
                   Xi2 == integer.Xi2;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(Yi1, Xi1, Yi2, Xi2);
        }
    }
    struct IntegerPair
    {
        public int Yi { get; init; }
        public int Xi { get; init; }
        public IntegerPair(int yi, int xi)
        {
            Yi = yi;
            Xi = xi;
        }

        public override string ToString()
        {
            return $"{Yi}\t{Xi}";
        }

        public override bool Equals(object? obj)
        {
            return obj is IntegerPair pairs &&
                   Yi == pairs.Yi &&
                   Xi == pairs.Xi;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(Yi, Xi);
        }
    }
    internal class GridGenerator
    {
        int Kx, Ky;
        List<List<Point2D>> yCoordinateLines = new List<List<Point2D>>();
        int nO;

        List<double> X = new List<double>();
        List<double> Y = new List<double>();
        List<Subarea> subareas = new List<Subarea>();
        List<Subinterval> xIntervalSplits = new List<Subinterval>();
        List<Subinterval> yIntervalSplits = new List<Subinterval>();
        List<Point2D> nodes = new List<Point2D>();
        List<TriangularElement> elems = new List<TriangularElement>();
        List<Bound> bounds = new List<Bound>();
        List<DirichletBoundaryCondition> bc1s = new List<DirichletBoundaryCondition>();
        List<NeumannBoundaryCondition> bc2s = new List<NeumannBoundaryCondition>();
        List<RobinBoundaryCondition> bc3s = new List<RobinBoundaryCondition>();
        private int gridСrushingX;
        private int gridСrushingY;

        Dictionary<QuadrupleInteger, List<int>> subareaBoundaries = new Dictionary<QuadrupleInteger, List<int>>();
        Dictionary<IntegerPair, int> subareaAnglePoints = new Dictionary<IntegerPair, int>();

        List<int> igEdge = new List<int>();
        List<int> jgEdge = new List<int>();
        List<List<Point2D>> ggEdge = new List<List<Point2D>>();

        List<List<int>> unknownsOfEdges = new List<List<int>>();
        List<List<int>> unknownsOfElements = new List<List<int>>();
        List<List<int>> unknownsOfNodes = new List<List<int>>();
        private bool isGridGenerated = false;
        public int NumberOfGlobalBasisFunctions { get; private set; }


        public List<TriangularElement> Elems { get => elems; }
        public List<Point2D> Nodes { get => nodes; }
        public List<int> IgEdge { get => igEdge; }
        public List<int> JgEdge { get => jgEdge; }
        public List<List<Point2D>> GgEdge { get => ggEdge; }


        public List<List<int>> UnknownsOfEdges { get
            {
                try
                {
                    if (isGridGenerated)
                        return unknownsOfEdges;
                    else
                    {
                        throw new Exception("An attempt is found to determine the global indices of the unknowns on an element prior to mesh generation.");
                    }

                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                }
                return null;
            }
        }
        public List<List<int>> UnknownsOfNodes
        {
            get
            {
                try
                {
                    if (isGridGenerated)
                        return unknownsOfNodes;
                    else
                    {
                        throw new Exception("An attempt is found to determine the global indices of the unknowns on an element prior to mesh generation.");
                    }

                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                }
                return null;
            }
        }
        public List<List<int>> UnknownsOfElements
        {
            get
            {
                try
                {
                    if (isGridGenerated)
                        return unknownsOfElements;
                    else
                    {
                        throw new Exception("An attempt is found to determine the global indices of the unknowns on an element prior to mesh generation.");
                    }

                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                }
                return null;
            }
        }
        public GridGenerator(string mediaGeometryFileName, string boundaryFileName)
        {
            try
            {
                StreamReader sr = new StreamReader(mediaGeometryFileName);

                List<int> x = sr.ReadLine().Trim().Split(" ").Select(x => Int32.Parse(x)).ToList();
                Kx = x[0]; Ky = x[1];

                for (int i = 0; i < Ky; i++)
                {
                    yCoordinateLines.Add(new List<Point2D>());
                    List<double> values = sr.ReadLine().Trim().Split(" ").Select(x => Double.Parse(x)).ToList();
                    for (int j = 0; j < Kx; j++)
                    {
                        yCoordinateLines[i].Add(new(values[j * 2], values[j * 2 + 1]));
                    }
                }

                nO = Int32.Parse(sr.ReadLine());

                for (int i = 0; i < nO; i++)
                {
                    List<int> subareaParams = sr.ReadLine().Replace('\t', ' ').Trim().Split(" ").Select(x => Int32.Parse(x)).ToList();
                    subareas.Add(new Subarea(subareaParams[0], subareaParams[1], subareaParams[2], subareaParams[3], subareaParams[4]));
                }

                List<double> intervalSplitsValues = sr.ReadLine().Trim().Split(" ").Select(x => Double.Parse(x)).ToList();
                for (int i = 0;i < Kx - 1;i++) 
                {
                    xIntervalSplits.Add(new((int)intervalSplitsValues[2 * i], intervalSplitsValues[2 * i + 1]));
                }

                intervalSplitsValues = sr.ReadLine().Trim().Split(" ").Select(x => Double.Parse(x)).ToList();
                for (int i = 0; i < Ky - 1; i++)
                {
                    yIntervalSplits.Add(new((int)intervalSplitsValues[2 * i], intervalSplitsValues[2 * i + 1]));
                }
                sr.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("The file " + mediaGeometryFileName + " could not be read:");
                Console.WriteLine(e.Message);
            }

            try
            {
                using (StreamReader sr = new StreamReader(boundaryFileName))
                {
                    string line;
                    while ((line = sr.ReadLine()) != null)
                    {
                        List<int> ints = new List<int>(line.Trim().Split(" ").ToList().Select(e => int.Parse(e)));
                        bounds.Add(new Bound(ints[0], ints[1], ints[2], ints[3], ints[4], ints[5]));
                    }
                    sr.Close();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine($"The file {boundaryFileName} could not be read:");
                Console.WriteLine(e.Message);
            }
        }
        public void GenerateGrid()
        {
            AddAnglePoints();
            for (int i = 0; i < nO; i++)
            {
                ProcessBoundary(subareas[i]);
                AddSubarea(subareas[i]);
            }
            NodeOrder();
            EdgeOrder();
            ElementsOrder();
            BasisFunctionsOrder();

            for (int i = 0; i < elems.Count; i++)
            {
                elems[i].BasisFunctionsType = 2; // Cubic Lagrange
            }
            isGridGenerated = true;
        }
        internal void DefineGlobalIndicesOfUnknownsOnElements()
        {
            for (int iElem = 0; iElem < elems.Count; iElem++)
            {
                List<int> globalIndices = new List<int>(new int[10]);
                globalIndices[0] = unknownsOfNodes[elems[iElem].Node0][0];
                globalIndices[1] = unknownsOfNodes[elems[iElem].Node1][0];
                globalIndices[2] = unknownsOfNodes[elems[iElem].Node2][0];
                List<int> edge1unknowns = GetUnknownsOfEdge(elems[iElem].Node0, elems[iElem].Node1);
                globalIndices[3] = edge1unknowns[0];
                globalIndices[4] = edge1unknowns[1];
                List<int> edge2unknowns = GetUnknownsOfEdge(elems[iElem].Node1, elems[iElem].Node2);
                globalIndices[5] = edge2unknowns[0];
                globalIndices[6] = edge2unknowns[1];
                List<int> edge3unknowns = GetUnknownsOfEdge(elems[iElem].Node0, elems[iElem].Node2);
                globalIndices[7] = edge3unknowns[1]; // | на третьем ребре локальные базисные функции в аналитических формулах нумеруются от вершины 3 к вершине 1
                globalIndices[8] = edge3unknowns[0]; // | в то же время в предложенном варианте в edge3unknowns базисные функции на ребре от меньшего к большему локальному индексу (от 1 к 3).
                globalIndices[9] = unknownsOfElements[iElem][0];
                
                elems[iElem].SetGlobalIndicesOfUnknowns(globalIndices);
            }
        }
        private void BasisFunctionsOrder()
        {
            for (int i = 0; i < jgEdge.Count; i++)
                unknownsOfEdges.Add(new List<int>());

            for (int i = 0; i < elems.Count; i++)
                unknownsOfElements.Add(new List<int>());

            for (int i = 0; i < nodes.Count; i++)
                unknownsOfNodes.Add(new List<int>());

            int unknownsCounter = 0;
            // !!!!!!!!!!!! это только для кубических лагранжевых функций на треугольниках
            for (int i = 0; i < nodes.Count; i++)
            {
                // для треугольников в кубическом лагранжевом
                // базисе с вершиной ассоциированна только
                // одна базисная функция. 
                for (int j = 0; j < 1; j++)
                {
                    unknownsOfNodes[i].Add(unknownsCounter);
                    unknownsCounter++;
                }

                for (int edgeIndex = igEdge[i]; edgeIndex < igEdge[i+1]; edgeIndex++)
                    for (int k = 0; k < ggEdge[edgeIndex].Count; k++)
                    {
                        unknownsOfEdges[edgeIndex].Add(unknownsCounter);
                        unknownsCounter++;
                    }

                foreach (var elem in (from element in elems where element.Node0 == i select element).ToList())
                {
                    // для треугольников в кубическом лагранжевом
                    // базисе с элементом ассоциированна только
                    // одна базисная функция. 
                    for (int j = 0; j < 1; j++) 
                    {
                        unknownsOfElements[elems.IndexOf(elem)].Add(unknownsCounter);
                        unknownsCounter++;
                    }
                }
            }
            NumberOfGlobalBasisFunctions = unknownsCounter;

            //for (int ii = 0; ii < elems.Count; ii++)
            //{
            //    if (ii == 1466 || ii == 1468)
            //    {
            //        Console.WriteLine($"ELEM NUMBER = {ii}");
            //        Console.WriteLine($"node0  ({elems[ii].Node0}): {unknownsOfNodes[elems[ii].Node0][0]} node1 ({elems[ii].Node1}): {unknownsOfNodes[elems[ii].Node1][0]} node2 ({elems[ii].Node2}): {unknownsOfNodes[elems[ii].Node2][0]}");

            //        List<int> edge1unknowns = GetUnknownsOfEdge(elems[ii].Node0, elems[ii].Node1);
            //        Console.Write($"Edge 1 unknowns ({elems[ii].Node0} {elems[ii].Node1}): ");
            //        for (int i = 0; i < edge1unknowns.Count; i++)
            //        {
            //            Console.Write($"{edge1unknowns[i]}\t");
            //        }
            //        Console.WriteLine();

            //        List<int> edge2unknowns = GetUnknownsOfEdge(elems[ii].Node1, elems[ii].Node2);
            //        Console.Write($"Edge 2 unknowns ({elems[ii].Node1} {elems[ii].Node2}): ");
            //        for (int i = 0; i < edge2unknowns.Count; i++)
            //        {
            //            Console.Write($"{edge2unknowns[i]}\t");
            //        }
            //        Console.WriteLine();

            //        List<int> edge3unknowns = GetUnknownsOfEdge(elems[ii].Node0, elems[ii].Node2);
            //        Console.Write($"Edge 3 unknowns ({elems[ii].Node0} {elems[ii].Node2}): ");
            //        for (int i = 0; i < edge3unknowns.Count; i++)
            //        {
            //            Console.Write($"{edge3unknowns[i]}\t");
            //        }
            //        Console.WriteLine();

            //        Console.Write("Inner element unknowns: ");
            //        Console.Write($"{unknownsOfElements[ii][0]}");

            //        Console.WriteLine();
            //        Console.WriteLine();
            //    }

            //}

        }

        private List<int> GetUnknownsOfEdge(int i, int j)
        {
            if (i > j)
            {
                int k = i;
                i = j;
                j = k;
            }

            for (int jj = igEdge[i]; jj < igEdge[i+1]; jj++)
            {
                if (jgEdge[jj] == j)
                {
                    return unknownsOfEdges[jj];
                }
            }
            return null;
        }

        private void ElementsOrder()
        {
            elems.Sort();
        }

        private void EdgeOrder()
        {
            List<List<int>> list = new List<List<int>>();
            int N = nodes.Count;

            for (int inodes = 0; inodes < N; inodes++)
                list.Add(new List<int>());

            for (int ielem = 0; ielem < elems.Count; ielem++)
            {
                for (int i = 0; i < ((Element)elems[ielem]).GetNumberOfVertexes(); i++)
                {
                    for (int j = 0; j < ((Element)elems[ielem]).GetNumberOfVertexes(); j++)
                    {
                        if (((Element)elems[ielem]).GetIndexOfVertex(i) < ((Element)elems[ielem]).GetIndexOfVertex(j) &&
                            !list[((Element)elems[ielem]).GetIndexOfVertex(i)].Contains(((Element)elems[ielem]).GetIndexOfVertex(j)))
                            list[((Element)elems[ielem]).GetIndexOfVertex(i)].Add(((Element)elems[ielem]).GetIndexOfVertex(j));
                    }
                }
            }

            int curSum = 0;
            igEdge.Add(curSum);
            for (int inodes = 0; inodes < N; inodes++)
            {
                list[inodes].Sort();
                for (int i = 0; i < list[inodes].Count; i++)
                {
                    jgEdge.Add(list[inodes][i]);
                    curSum++;
                }
                igEdge.Add(curSum);
            }

            for (int i = 0; i < jgEdge.Count; i++)
            {
                ggEdge.Add(new List<Point2D>());
            }

            for (int i = 0; i < N; i++)
            {
                for (int j = igEdge[i]; j < igEdge[i+1]; j++)
                {
                    List<Point2D> points = SegmentSeparator.GetPointsOfLine(nodes[i], nodes[jgEdge[j]], new(3, 1));
                    points.Remove(points[0]);
                    points.Remove(points.Last());
                    ggEdge[j].Add(points[0]);
                    ggEdge[j].Add(points[1]);
                }
            }

            //for (int i = 0; i < N; i++)
            //{
            //    Console.Write($"i = {i}:\t");
            //    for (int j = igEdge[i]; j < igEdge[i + 1]; j++)
            //    {
            //        Console.Write($"{jgEdge[j]}\t\t");
            //    }
            //    Console.WriteLine();
            //}

        }

        private void NodeOrder()
        {
            int digitNumber = 10;

            List<Point2D> newNodes = new List<Point2D>();
            for (int i = 0; i < nodes.Count; i++)
                nodes[i] = new Point2D(Math.Round(nodes[i].X, digitNumber), Math.Round(nodes[i].Y, digitNumber));
            
            newNodes.AddRange(nodes);
            newNodes.Sort();


            Dictionary<int, int> keyValuePairs = new Dictionary<int, int>();

            for (int i = 0; i < nodes.Count; i++)
                keyValuePairs.Add(i, newNodes.IndexOf(nodes[i]));

            foreach (var item in elems)
            {
                List<int> list = new List<int>();
                list.Add(keyValuePairs[item.Node0]);
                list.Add(keyValuePairs[item.Node1]);
                list.Add(keyValuePairs[item.Node2]);
                list.Sort();
                item.Node0 = list[0];
                item.Node1 = list[1];
                item.Node2 = list[2];
            }
            nodes.Clear();
            nodes.AddRange(newNodes);

            for (int i = 0; i < bc1s.Count; i++)
            {
                bc1s[i].NodeNumber1 = keyValuePairs[bc1s[i].NodeNumber1];
                bc1s[i].NodeNumber2 = keyValuePairs[bc1s[i].NodeNumber2];
            }
            for (int i = 0; i < bc2s.Count; i++)
            {
                bc2s[i].NodeNumber1 = keyValuePairs[bc2s[i].NodeNumber1];
                bc2s[i].NodeNumber2 = keyValuePairs[bc2s[i].NodeNumber2];
            }
            for (int i = 0; i < bc3s.Count; i++)
            {
                bc3s[i].NodeNumber1 = keyValuePairs[bc3s[i].NodeNumber1];
                bc3s[i].NodeNumber2 = keyValuePairs[bc3s[i].NodeNumber2];
            }
        }
 
        private void AddAnglePoints()
        {
            for (int i = 0; i < yCoordinateLines.Count; i++)
            {
                for (int j = 0; j < yCoordinateLines[i].Count; j++)
                {
                    nodes.Add(yCoordinateLines[i][j]);
                    subareaAnglePoints.Add(new(i,j), nodes.Count - 1);
                }
            }
        }
        private void AddSubarea(Subarea subarea)
        {

            Point2D ldPoint = yCoordinateLines[subarea.Y0][subarea.X0];
            Point2D rdPoint = yCoordinateLines[subarea.Y0][subarea.X1];
            Point2D luPoint = yCoordinateLines[subarea.Y1][subarea.X0];
            Point2D ruPoint = yCoordinateLines[subarea.Y1][subarea.X1];

            int ldPointIndex = subareaAnglePoints[new(subarea.Y0, subarea.X0)];
            int rdPointIndex = subareaAnglePoints[new(subarea.Y0, subarea.X1)];
            int luPointIndex = subareaAnglePoints[new(subarea.Y1, subarea.X0)];
            int ruPointIndex = subareaAnglePoints[new(subarea.Y1, subarea.X1)];
            Subinterval xIntervalSplit = xIntervalSplits[subarea.X0]; 
            Subinterval yIntervalSplit = yIntervalSplits[subarea.Y0];

            List<Point2D> leftLinePoints = new List<Point2D>();
            List<Point2D> rightLinePoints = new List<Point2D>();
            List<Point2D> upLinePoints = new List<Point2D>();
            List<Point2D> downLinePoints = new List<Point2D>();
            List<int> leftLineIndexes = new List<int>();
            List<int> rightLineIndexes = new List<int>();
            List<int> upLineIndexes = new List<int>();
            List<int> downLineIndexes = new List<int>();
            downLineIndexes.AddRange(subareaBoundaries.GetValueOrDefault(new(subarea.Y0, subarea.X0, subarea.Y0, subarea.X1)));
            upLineIndexes.AddRange(subareaBoundaries.GetValueOrDefault(new(subarea.Y1, subarea.X0, subarea.Y1, subarea.X1)));
            leftLineIndexes.AddRange(subareaBoundaries.GetValueOrDefault(new(subarea.Y0, subarea.X0, subarea.Y1, subarea.X0)));
            rightLineIndexes.AddRange(subareaBoundaries.GetValueOrDefault(new(subarea.Y0, subarea.X1, subarea.Y1, subarea.X1)));

            downLinePoints.AddRange(downLineIndexes.Select(x => nodes[x]));
            upLinePoints.AddRange(upLineIndexes.Select(x => nodes[x]));
            leftLinePoints.AddRange(leftLineIndexes.Select(x => nodes[x]));
            rightLinePoints.AddRange(rightLineIndexes.Select(x => nodes[x]));

            for (int i = 0; i < downLineIndexes.Count; i++)
            {
                Point2D downPointOfCurrentRightLine = downLinePoints[i];
                Point2D upPointOfCurrentRightLine = upLinePoints[i];
                List<Point2D> curRightLinePoints = SegmentSeparator.GetPointsOfLine(downPointOfCurrentRightLine, upPointOfCurrentRightLine, yIntervalSplit);
                curRightLinePoints.RemoveAt(0);
                curRightLinePoints.RemoveAt(curRightLinePoints.Count - 1);

                for (int j = 0; j < leftLineIndexes.Count; j++)
                {
                    nodes.Add(curRightLinePoints[j]);
                    if (j == 0 && i == 0)
                    {
                        AddTwoTriangularElements(ldPointIndex, downLineIndexes[i], leftLineIndexes[j], nodes.Count - 1, subarea.SubariaNumber);
                    }
                    else if (i == 0)
                    {
                        AddTwoTriangularElements(leftLineIndexes[j - 1], nodes.Count - 2, leftLineIndexes[j], nodes.Count - 1, subarea.SubariaNumber);
                    }
                    else if (j == 0)
                    {
                        AddTwoTriangularElements(downLineIndexes[i - 1], downLineIndexes[i], nodes.Count - 1 - curRightLinePoints.Count, nodes.Count - 1, subarea.SubariaNumber);
                    }
                    else
                    {
                        AddTwoTriangularElements(nodes.Count - 1 - curRightLinePoints.Count - 1, nodes.Count - 1 - 1, nodes.Count - 1 - curRightLinePoints.Count, nodes.Count - 1, subarea.SubariaNumber);
                    }
                }
                if (i == 0)
                {
                    AddTwoTriangularElements(leftLineIndexes[leftLineIndexes.Count - 1], nodes.Count - 1, luPointIndex, upLineIndexes[0], subarea.SubariaNumber);
                }
                else 
                {
                    AddTwoTriangularElements(nodes.Count - 1 - curRightLinePoints.Count, nodes.Count - 1, upLineIndexes[i - 1], upLineIndexes[i], subarea.SubariaNumber);
                }
            }

            AddTwoTriangularElements(downLineIndexes[downLineIndexes.Count - 1], rdPointIndex, nodes.Count - rightLineIndexes.Count, rightLineIndexes[0], subarea.SubariaNumber);
            for (int j = 1; j < rightLineIndexes.Count; j++)
            {
                AddTwoTriangularElements(nodes.Count - rightLinePoints.Count + j - 1, rightLineIndexes[j - 1], nodes.Count - rightLineIndexes.Count + j, rightLineIndexes[j], subarea.SubariaNumber);
            }
            AddTwoTriangularElements(nodes.Count - 1, rightLineIndexes[rightLineIndexes.Count - 1], upLineIndexes[upLineIndexes.Count - 1], ruPointIndex, subarea.SubariaNumber);
        }
        private void ProcessBoundary(Subarea subarea)
        {
            ProcessSubareaBoundaryPoints(subarea.Y0, subarea.X0, subarea.Y0, subarea.X1, xIntervalSplits[subarea.X0], subarea.SubariaNumber); 
            ProcessSubareaBoundaryPoints(subarea.Y1, subarea.X0, subarea.Y1, subarea.X1, xIntervalSplits[subarea.X0], subarea.SubariaNumber);
            ProcessSubareaBoundaryPoints(subarea.Y0, subarea.X0, subarea.Y1, subarea.X0, yIntervalSplits[subarea.Y0], subarea.SubariaNumber); 
            ProcessSubareaBoundaryPoints(subarea.Y0, subarea.X1, subarea.Y1, subarea.X1, yIntervalSplits[subarea.Y0], subarea.SubariaNumber); 
        }
        private void ProcessSubareaBoundaryPoints(int y1index, int x1index, int y2index, int x2index, Subinterval intervalSplit, int subareaNumber)
        {
            List<Point2D> linePoints = new List<Point2D>();
            if (!subareaBoundaries.ContainsKey(new(y1index, x1index, y2index, x2index)))
            {
                linePoints.AddRange(SegmentSeparator.GetPointsOfLine(yCoordinateLines[y1index][x1index], yCoordinateLines[y2index][x2index], intervalSplit));
                List<int> indexes = new List<int>();
                for (int i = 1; i < linePoints.Count - 1; i++)
                {
                    nodes.Add(linePoints[i]);
                    indexes.Add(nodes.Count - 1);
                }
                subareaBoundaries.Add(new(y1index, x1index, y2index, x2index), indexes);
                AddBoundaryCondition(subareaNumber, x1index, x2index, y1index, y2index);
            }
        }

        private void AddBoundaryCondition(int subareaNumber, int x1index, int x2index, int y1index, int y2index)
        {
            foreach (var bound in bounds)
            {
                if (bound.Equals(new Bound(0, 0, x1index, x2index, y1index, y2index)))
                {
                    List<int> indexes = new List<int>();
                    indexes.Add(subareaAnglePoints[new IntegerPair(y1index, x1index)]);
                    indexes.AddRange(subareaBoundaries[new(y1index, x1index, y2index, x2index)]);
                    indexes.Add(subareaAnglePoints[new IntegerPair(y2index, x2index)]);

                    if (bound.BoundaryType == 1)
                        AddBc1(indexes, bound.FormulaNumber);
                    if (bound.BoundaryType == 2)
                        AddBc2(indexes, bound.FormulaNumber);
                    if (bound.BoundaryType == 3)
                        AddBc3(indexes, bound.FormulaNumber);

                    bounds.Remove(bound);
                    break;
                }
            }
        }
        private void AddBc3(List<int> indexes, int formulaNumber)
        {
            for (int i = 0; i < indexes.Count - 1; i++)
            {
                bc3s.Add(new RobinBoundaryCondition(indexes[i], indexes[i + 1], formulaNumber));
            }
        }
        private void AddBc2(List<int> indexes, int formulaNumber)
        {
            for (int i = 0; i < indexes.Count - 1; i++)
            {
                bc2s.Add(new NeumannBoundaryCondition(indexes[i], indexes[i + 1], formulaNumber));
            }
        }
        private void AddBc1(List<int> indexes, int formulaNumber)
        {
            for (int i = 0; i < indexes.Count - 1; i++)
            {
                bc1s.Add(new DirichletBoundaryCondition(indexes[i], indexes[i + 1], formulaNumber));
            }
        }
        private void AddTwoTriangularElements(int leftDownNodesIndex, int RigthDownNodesIndex, int leftUpNodesIndex, int RightUpNodesIndex, int subareaNumber)
        {
            if (GetLength(nodes[leftDownNodesIndex], nodes[RightUpNodesIndex]) > GetLength(nodes[RigthDownNodesIndex], nodes[leftUpNodesIndex]))
            {
                elems.Add(new(leftDownNodesIndex, RigthDownNodesIndex, leftUpNodesIndex, subareaNumber));
                elems.Add(new(RigthDownNodesIndex, RightUpNodesIndex, leftUpNodesIndex, subareaNumber));
            }
            else
            { 
                elems.Add(new(leftDownNodesIndex, RigthDownNodesIndex, RightUpNodesIndex, subareaNumber));
                elems.Add(new(leftDownNodesIndex, RightUpNodesIndex, leftUpNodesIndex, subareaNumber));
            }
        }
        private double GetLength(Point2D p1, Point2D p2)
        {
            return Math.Sqrt((p1.X - p2.X) * (p1.X - p2.X) + (p1.Y - p2.Y) * (p1.Y - p2.Y));
        }
        public void WriteGrid()
        {
            WriteElements();
            WritePoints();
            WriteBound();
        }

        private void WriteElements()
        {
            if (elems.Count != 0)
            {
                try
                {
                    File.WriteAllText("elems.txt", string.Empty);
                    using (StreamWriter sw = new StreamWriter("elems.txt", false))
                    {
                        sw.WriteLine(elems.Count);
                        foreach (var element in elems)
                            sw.WriteLine(element.ToString());
                        sw.Close();
                    }
                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                }
            }
        }

        private void WritePoints()
        {
            if (nodes.Count != 0)
            {
                try
                {
                    File.WriteAllText("nodes.txt", string.Empty);
                    using (StreamWriter sw = new StreamWriter("nodes.txt", false))
                    {
                        sw.WriteLine(nodes.Count);
                        foreach (var node in nodes)
                            sw.WriteLine(node.ToString());
                        sw.Close();
                    }
                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                }
            }
        }
        private void WriteBound()
        {
            try
            {
                File.WriteAllText("bc1.txt", string.Empty);
                using (StreamWriter sw = new StreamWriter("bc1.txt", false))
                {
                    sw.WriteLine(bc1s.Count);
                    foreach (var bc in bc1s)
                        sw.WriteLine(bc.ToString());
                    sw.Close();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
            try
            {
                File.WriteAllText("bc2.txt", string.Empty);
                using (StreamWriter sw = new StreamWriter("bc2.txt", false))
                {
                    sw.WriteLine(bc2s.Count);
                    foreach (var bc in bc2s)
                        sw.WriteLine(bc.ToString());
                    sw.Close();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
            try
            {
                File.WriteAllText("bc3.txt", string.Empty);
                using (StreamWriter sw = new StreamWriter("bc3.txt", false))
                {
                    sw.WriteLine(bc3s.Count);
                    foreach (var bc in bc3s)
                        sw.WriteLine(bc.ToString());
                    sw.Close();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
        }


    }
}
