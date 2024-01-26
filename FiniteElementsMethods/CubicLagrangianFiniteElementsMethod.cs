using FEM.BoundaryСonditions;
using FEM.Elements;
using FEM.Points;
using System.ComponentModel.Design;
using System.Dynamic;
using System.Reflection.Metadata.Ecma335;

namespace FEM.FiniteElementsMethods
{
    public struct BasisFunctionTerm
    {
        public double NumericalMultiplier { get; init; }
        public int L1Degree { get; init; }
        public int L2Degree { get; init; }
        public int L3Degree { get; init; }
        public int Alpha1Degree { get; init; }
        public int Alpha2Degree { get; init; }
        public int Alpha3Degree { get; init; }

        public BasisFunctionTerm(double numericalMultiplier, int L1Degree, int L2Degree, int L3Degree, int alpha1Degree, int alpha2Degree, int alpha3Degree)
        { 
            NumericalMultiplier = numericalMultiplier;
            this.L1Degree = L1Degree;
            this.L2Degree = L2Degree;
            this.L3Degree = L3Degree;
            Alpha1Degree = alpha1Degree;
            Alpha2Degree = alpha2Degree;
            Alpha3Degree = alpha3Degree;
        }
    }

    internal class CubicLagrangianFiniteElementsMethod : FiniteElementsMethod
    {
        List<List<BasisFunctionTerm>> psi = new List<List<BasisFunctionTerm>>();
        List<List<BasisFunctionTerm>> derivativePsi = new List<List<BasisFunctionTerm>>();

        public CubicLagrangianFiniteElementsMethod()
        {
            List<BasisFunctionTerm> Psi1 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> Psi2 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> Psi3 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> Psi4 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> Psi5 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> Psi6 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> Psi7 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> Psi8 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> Psi9 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> Psi10 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> DerivativePsi1 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> DerivativePsi2 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> DerivativePsi3 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> DerivativePsi4 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> DerivativePsi5 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> DerivativePsi6 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> DerivativePsi7 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> DerivativePsi8 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> DerivativePsi9 = new List<BasisFunctionTerm>();
            List<BasisFunctionTerm> DerivativePsi10 = new List<BasisFunctionTerm>();

            Psi1.Add(new(9.0 / 2.0, 3, 0, 0, 0, 0, 0)); Psi1.Add(new(-9.0 / 2.0, 2, 0, 0, 0, 0, 0)); Psi1.Add(new(1.0, 1, 0, 0, 0, 0, 0));
            Psi2.Add(new(9.0 / 2.0, 0, 3, 0, 0, 0, 0)); Psi2.Add(new(-9.0 / 2.0, 0, 2, 0, 0, 0, 0)); Psi2.Add(new(1.0, 0, 1, 0, 0, 0, 0));
            Psi3.Add(new(9.0 / 2.0, 0, 0, 3, 0, 0, 0)); Psi3.Add(new(-9.0 / 2.0, 0, 0, 2, 0, 0, 0)); Psi3.Add(new(1.0, 0, 0, 1, 0, 0, 0));
            Psi4.Add(new(27.0 / 2.0, 2, 1, 0, 0, 0, 0)); Psi4.Add(new(-9.0 / 2.0, 1, 1, 0, 0, 0, 0));
            Psi5.Add(new(27.0 / 2.0, 1, 2, 0, 0, 0, 0)); Psi5.Add(new(-9.0 / 2.0, 1, 1, 0, 0, 0, 0));
            Psi6.Add(new(27.0 / 2.0, 0, 2, 1, 0, 0, 0)); Psi6.Add(new(-9.0 / 2.0, 0, 1, 1, 0, 0, 0));
            Psi7.Add(new(27.0 / 2.0, 0, 1, 2, 0, 0, 0)); Psi7.Add(new(-9.0 / 2.0, 0, 1, 1, 0, 0, 0));
            Psi8.Add(new(27.0 / 2.0, 1, 0, 2, 0, 0, 0)); Psi8.Add(new(-9.0 / 2.0, 1, 0, 1, 0, 0, 0));
            Psi9.Add(new(27.0 / 2.0, 2, 0, 1, 0, 0, 0)); Psi9.Add(new(-9.0 / 2.0, 1, 0, 1, 0, 0, 0));
            Psi10.Add(new(27.0, 1, 1, 1, 0, 0, 0));
            psi.Add(Psi1);
            psi.Add(Psi2);
            psi.Add(Psi3);
            psi.Add(Psi4);
            psi.Add(Psi5);
            psi.Add(Psi6);
            psi.Add(Psi7);
            psi.Add(Psi8); 
            psi.Add(Psi9); 
            psi.Add(Psi10);

            DerivativePsi1.Add(new(27.0 / 2.0, 2, 0, 0, 1, 0, 0)); DerivativePsi1.Add(new(-9.0, 1, 0, 0, 1, 0, 0)); DerivativePsi1.Add(new(1.0, 0, 0, 0, 1, 0, 0));
            DerivativePsi2.Add(new(27.0 / 2.0, 0, 2, 0, 0, 1, 0)); DerivativePsi2.Add(new(-9.0, 0, 1, 0, 0, 1, 0)); DerivativePsi2.Add(new(1.0, 0, 0, 0, 0, 1, 0));
            DerivativePsi3.Add(new(27.0 / 2.0, 0, 0, 2, 0, 0, 1)); DerivativePsi3.Add(new(-9.0, 0, 0, 1, 0, 0, 1)); DerivativePsi3.Add(new(1.0, 0, 0, 0, 0, 0, 1));
            DerivativePsi4.Add(new(27.0, 1, 1, 0, 1, 0, 0)); DerivativePsi4.Add(new(27.0 / 2.0, 2, 0, 0, 0, 1, 0)); DerivativePsi4.Add(new(-9.0 / 2.0, 0, 1, 0, 1, 0, 0)); DerivativePsi4.Add(new(-9.0 / 2.0, 1, 0, 0, 0, 1, 0));
            DerivativePsi5.Add(new(27.0, 1, 1, 0, 0, 1, 0)); DerivativePsi5.Add(new(27.0 / 2.0, 0, 2, 0, 1, 0, 0)); DerivativePsi5.Add(new(-9.0 / 2.0, 0, 1, 0, 1, 0, 0)); DerivativePsi5.Add(new(-9.0 / 2.0, 1, 0, 0, 0, 1, 0));
            DerivativePsi6.Add(new(27.0, 0, 1, 1, 0 , 1, 0)); DerivativePsi6.Add(new(27.0 / 2.0, 0, 2, 0 ,0,0,1)); DerivativePsi6.Add(new(-9.0 / 2.0, 0, 0, 1, 0, 1, 0)); DerivativePsi6.Add(new(-9.0 / 2.0, 0, 1, 0, 0, 0, 1));
            DerivativePsi7.Add(new(27.0, 0, 1, 1, 0, 0, 1)); DerivativePsi7.Add(new(27.0 / 2.0, 0, 0, 2, 0, 1, 0)); DerivativePsi7.Add(new(-9.0 / 2.0, 0, 0, 1, 0, 1, 0)); DerivativePsi7.Add(new(-9.0 / 2.0, 0, 1, 0, 0, 0, 1));
            DerivativePsi8.Add(new(27.0, 1, 0, 1, 0, 0, 1)); DerivativePsi8.Add(new(27.0 / 2.0, 0, 0, 2, 1, 0, 0)); DerivativePsi8.Add(new(-9.0 / 2.0, 0, 0, 1, 1, 0, 0)); DerivativePsi8.Add(new(-9.0 / 2.0, 1, 0, 0, 0, 0, 1));
            DerivativePsi9.Add(new(27.0, 1, 0, 1, 1, 0, 0)); DerivativePsi9.Add(new(27.0 / 2.0, 2, 0, 0, 0, 0, 1)); DerivativePsi9.Add(new(-9.0 / 2.0, 0, 0, 1, 1, 0, 0)); DerivativePsi9.Add(new(-9.0 / 2.0, 1, 0, 0, 0, 0, 1));
            DerivativePsi10.Add(new(27.0, 0, 1, 1, 1, 0, 0)); DerivativePsi10.Add(new(27.0, 1, 0, 1, 0, 1, 0)); DerivativePsi10.Add(new(27.0, 1, 1, 0, 0, 0, 1));
            derivativePsi.Add(DerivativePsi1);
            derivativePsi.Add(DerivativePsi2);
            derivativePsi.Add(DerivativePsi3);
            derivativePsi.Add(DerivativePsi4);
            derivativePsi.Add(DerivativePsi5);
            derivativePsi.Add(DerivativePsi6);
            derivativePsi.Add(DerivativePsi7);
            derivativePsi.Add(DerivativePsi8); 
            derivativePsi.Add(DerivativePsi9); 
            derivativePsi.Add(DerivativePsi10);
        }

        protected override void AddElementsToGlobalSLAE()
        {
            for (int ielem = 0; ielem < elems.Count; ielem++)
            {
                TriangularElement element = elems[ielem];

                AddLocalStiffnessMatrix(element);
                AddLocalMassMatrix(element);
                AddToVector(element);
            }
        }

        protected override void AddLocalMassMatrix(Element element)
        {
            TriangularElement triangularElement = (TriangularElement)element;
            for (int i = 0; i < triangularElement.GetNumberOfUnknowns(); i++)
            {
                for (int j = 0; j < triangularElement.GetNumberOfUnknowns(); j++)
                {
                    AddElementToGlobal(triangularElement.GetIndexOfUnknown(i), triangularElement.GetIndexOfUnknown(j), problem.Gamma(element.Subarea) * CMatrixComponent(i, j, triangularElement));
                }
            }
        }

        private double CMatrixComponent(int iLocal, int jLocal, TriangularElement triangularElement)
        {
            double component = 0;
            for (int i = 0; i < psi[iLocal].Count; i++)
            {
                for (int j = 0; j < psi[jLocal].Count; j++)
                {
                    double numericalMultiplier = psi[iLocal][i].NumericalMultiplier * psi[jLocal][j].NumericalMultiplier;
                    int L1 = psi[iLocal][i].L1Degree + psi[jLocal][j].L1Degree;
                    int L2 = psi[iLocal][i].L2Degree + psi[jLocal][j].L2Degree;
                    int L3 = psi[iLocal][i].L3Degree + psi[jLocal][j].L3Degree;
                    int alpha1 = psi[iLocal][i].Alpha1Degree + psi[jLocal][j].Alpha1Degree;
                    int alpha2 = psi[iLocal][i].Alpha2Degree + psi[jLocal][j].Alpha2Degree;
                    int alpha3 = psi[iLocal][i].Alpha3Degree + psi[jLocal][j].Alpha3Degree;
                    component += GetIntegral(new BasisFunctionTerm(numericalMultiplier, L1, L2, L3, alpha1, alpha2, alpha3), triangularElement);
                }
            }
            return component;
        }

        protected override void AddLocalStiffnessMatrix(Element element)
        {
            TriangularElement triangularElement = (TriangularElement)element;
            for (int i = 0; i < triangularElement.GetNumberOfUnknowns(); i++)
            {
                for (int j = 0; j < triangularElement.GetNumberOfUnknowns(); j++)
                {
                    AddElementToGlobal(triangularElement.GetIndexOfUnknown(i), triangularElement.GetIndexOfUnknown(j), problem.Lambda(element.Subarea) * StiffnessMatrixComponent(i, j, triangularElement));
                }
            }
        }

        private double StiffnessMatrixComponent(int iLocal, int jLocal, TriangularElement triangularElement)
        {
            double component = 0;

            for (int i = 0; i < derivativePsi[iLocal].Count; i++)
            {
                for (int j = 0; j < derivativePsi[jLocal].Count; j++)
                {
                    double numericalMultiplier = derivativePsi[iLocal][i].NumericalMultiplier * derivativePsi[jLocal][j].NumericalMultiplier;
                    int L1 = derivativePsi[iLocal][i].L1Degree + derivativePsi[jLocal][j].L1Degree;
                    int L2 = derivativePsi[iLocal][i].L2Degree + derivativePsi[jLocal][j].L2Degree;
                    int L3 = derivativePsi[iLocal][i].L3Degree + derivativePsi[jLocal][j].L3Degree;
                    int alpha1 = derivativePsi[iLocal][i].Alpha1Degree + derivativePsi[jLocal][j].Alpha1Degree;
                    int alpha2 = derivativePsi[iLocal][i].Alpha2Degree + derivativePsi[jLocal][j].Alpha2Degree;
                    int alpha3 = derivativePsi[iLocal][i].Alpha3Degree + derivativePsi[jLocal][j].Alpha3Degree;
                    component += GetIntegral(new BasisFunctionTerm(numericalMultiplier, L1, L2, L3, alpha1, alpha2, alpha3), triangularElement, true) +
                        GetIntegral(new BasisFunctionTerm(numericalMultiplier, L1, L2, L3, alpha1, alpha2, alpha3), triangularElement, false);
                }
            }

            return component;
        }

        private double GetIntegral(BasisFunctionTerm basisFunctionTerm, TriangularElement triangularElement, bool isX)
        {
            double x1 = ((Point2D)nodes[triangularElement.Node0]).X;
            double y1 = ((Point2D)nodes[triangularElement.Node0]).Y;
            double x2 = ((Point2D)nodes[triangularElement.Node1]).X;
            double y2 = ((Point2D)nodes[triangularElement.Node1]).Y;
            double x3 = ((Point2D)nodes[triangularElement.Node2]).X;
            double y3 = ((Point2D)nodes[triangularElement.Node2]).Y;
            double detD = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
            double alpha1 = (y2 - y3) / detD;
            double alpha2 = (y3 - y1) / detD;
            double alpha3 = (y1 - y2) / detD;
            if (!isX)
            {
                alpha1 = (x3 - x2) / detD;
                alpha2 = (x1 - x3) / detD;
                alpha3 = (x2 - x1) / detD;
            }

            return basisFunctionTerm.NumericalMultiplier * 
                Math.Pow(alpha1, basisFunctionTerm.Alpha1Degree) *
                Math.Pow(alpha2, basisFunctionTerm.Alpha2Degree) *
                Math.Pow(alpha3, basisFunctionTerm.Alpha3Degree) *
                GetFactorial(basisFunctionTerm.L1Degree) *
                GetFactorial(basisFunctionTerm.L2Degree) *
                GetFactorial(basisFunctionTerm.L3Degree) *
                Math.Abs(detD) / 
                GetFactorial(basisFunctionTerm.L1Degree + basisFunctionTerm.L2Degree + basisFunctionTerm.L3Degree + 2);
        }

        private double GetIntegral(BasisFunctionTerm basisFunctionTerm, TriangularElement triangularElement)
        {
            double x1 = ((Point2D)nodes[triangularElement.Node0]).X;
            double y1 = ((Point2D)nodes[triangularElement.Node0]).Y;
            double x2 = ((Point2D)nodes[triangularElement.Node1]).X;
            double y2 = ((Point2D)nodes[triangularElement.Node1]).Y;
            double x3 = ((Point2D)nodes[triangularElement.Node2]).X;
            double y3 = ((Point2D)nodes[triangularElement.Node2]).Y;
            double detD = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);

            return basisFunctionTerm.NumericalMultiplier *
                GetFactorial(basisFunctionTerm.L1Degree) *
                GetFactorial(basisFunctionTerm.L2Degree) *
                GetFactorial(basisFunctionTerm.L3Degree) *
                Math.Abs(detD) /
                GetFactorial(basisFunctionTerm.L1Degree + basisFunctionTerm.L2Degree + basisFunctionTerm.L3Degree + 2);
        }

        private int GetFactorial(int n)
        {
            if (n <= 0) return 1;
            else
            { 
                return n * GetFactorial(n - 1);
            }
        }
        protected override void AddToVector(Element element)
        {
            TriangularElement triangularElement = (TriangularElement)element;
            for (int i = 0; i < triangularElement.GetNumberOfUnknowns(); i++)
            {
                for (int j = 0; j < triangularElement.GetNumberOfUnknowns(); j++)
                {
                    b[triangularElement.GetIndexOfUnknown(i)] += problem.RightPart(GetPoint(triangularElement, j),element.Subarea) * CMatrixComponent(i, j, triangularElement);
                }
                //Console.WriteLine($"Y({GetPoint(triangularElement, i)}) = {problem.RightPart(GetPoint(triangularElement, i), element.Subarea)}");
            }

        }

        private Point2D GetPoint(TriangularElement element, int iLocal)
        {
            if (iLocal == 0 || iLocal == 1 || iLocal == 2)
                return nodes[element.GetIndexOfVertex(iLocal)];
            else if (iLocal != 9)
            {
                return GetEdgePointByLocal(element, iLocal);
            }
            else
            { 
                Point2D p = new Point2D(0, 0);
                p.X = (nodes[element.Node0].X + nodes[element.Node1].X + nodes[element.Node2].X) / 3;
                p.Y = (nodes[element.Node0].Y + nodes[element.Node1].Y + nodes[element.Node2].Y) / 3;
                return p;
            }
        }

        private Point2D GetEdgePointByLocal(TriangularElement element, int iLocal)
        {

            if (iLocal == 3)
                return GetPointsOfEdge(element.Node0, element.Node1)[0];
            else if (iLocal == 4)
                return GetPointsOfEdge(element.Node0, element.Node1)[1];
            else if (iLocal == 5)
                return GetPointsOfEdge(element.Node1, element.Node2)[0];
            else if (iLocal == 6)
                return GetPointsOfEdge(element.Node1, element.Node2)[1];
            else if (iLocal == 7)
                return GetPointsOfEdge(element.Node0, element.Node2)[1];
            else if (iLocal == 8)
                return GetPointsOfEdge(element.Node0, element.Node2)[0];
            else
                return null;
        }

        protected override void GetGrid() // Это если нужно прочитать классом DataReader уже готовую сетку
        {
            throw new NotImplementedException();
        }

        protected override int GetNumberOfGlobalBasisFunctions()
        {
            throw new NotImplementedException();
        }

        protected override void AddBc2()
        {
            List<NeumannBoundaryCondition> bc2 = DataReader.ReadBC2FromFile();
            for (int i = 0; i < bc2.Count; i++)
            {
                double mes = GetLenght(bc2[i].NodeNumber1, bc2[i].NodeNumber2);
                double[,] C = {
                    {128.0 * mes / 1680.0,99.0 * mes / 1680.0, -36.0 * mes / 1680.0 , 19.0 * mes / 1680.0},
                    {99.0 * mes / 1680.0, 648.0 * mes / 1680.0, -81.0 * mes / 1680.0, -36.0 * mes / 1680.0},
                    {-36.0 * mes / 1680.0, -81.0 * mes / 1680.0, 648.0 * mes / 1680.0, 99.0 * mes / 1680.0},
                    {19.0 * mes / 1680.0, -36.0 * mes / 1680.0, 99.0 * mes / 1680.0, 128.0 * mes / 1680.0}
                };

                List<Point2D> edgePoints = GetPointsOfEdge(bc2[i].NodeNumber1, bc2[i].NodeNumber2);
                List<int> edgeUnknowns = GetUnknownsOfEdge(bc2[i].NodeNumber1, bc2[i].NodeNumber2);
                
                b[unknownsOfNodes[bc2[i].NodeNumber1][0]] += C[0, 0] * problem.Theta(nodes[bc2[i].NodeNumber1], bc2[i].FormulaNumber) +
                    C[0, 1] * problem.Theta(edgePoints[0], bc2[i].FormulaNumber) +
                    C[0, 2] * problem.Theta(edgePoints[1], bc2[i].FormulaNumber) +
                    C[0, 3] * problem.Theta(nodes[bc2[i].NodeNumber2], bc2[i].FormulaNumber);
                
                b[edgeUnknowns[0]] += C[1, 0] * problem.Theta(nodes[bc2[i].NodeNumber1], bc2[i].FormulaNumber) +
                    C[1, 1] * problem.Theta(edgePoints[0], bc2[i].FormulaNumber) +
                    C[1, 2] * problem.Theta(edgePoints[1], bc2[i].FormulaNumber) +
                    C[1, 3] * problem.Theta(nodes[bc2[i].NodeNumber2], bc2[i].FormulaNumber);

                b[edgeUnknowns[1]] += C[2, 0] * problem.Theta(nodes[bc2[i].NodeNumber1], bc2[i].FormulaNumber) +
                    C[2, 1] * problem.Theta(edgePoints[0], bc2[i].FormulaNumber) +
                    C[2, 2] * problem.Theta(edgePoints[1], bc2[i].FormulaNumber) +
                    C[2, 3] * problem.Theta(nodes[bc2[i].NodeNumber2], bc2[i].FormulaNumber);

                b[unknownsOfNodes[bc2[i].NodeNumber2][0]] += C[3, 0] * problem.Theta(nodes[bc2[i].NodeNumber1], bc2[i].FormulaNumber) +
                    C[3, 1] * problem.Theta(edgePoints[0], bc2[i].FormulaNumber) +
                    C[3, 2] * problem.Theta(edgePoints[1], bc2[i].FormulaNumber) +
                    C[3, 3] * problem.Theta(nodes[bc2[i].NodeNumber2], bc2[i].FormulaNumber);

            }

        }

        

        protected override void AddConcentratedSources()
        {
            throw new NotImplementedException();
        }

        protected override void AddBc3()
        {
            throw new NotImplementedException();
        }

        public override double GetSolution(Point2D point, TriangularElement element)
        {
            double result = 0;
            for (int i = 0; i < psi.Count; i++)
            {
                result += Q[element.GetIndexOfUnknown(i)] * GetPsi(i, element, point);
            }

            return result;
        }
        public double GetPsi(int i, TriangularElement element, Point2D point)
        {
            double result = 0;
            for (int j = 0; j < psi[i].Count; j++)
            {
                double multiplier = psi[i][j].NumericalMultiplier;

                double L1 = GetL(1, point, element);
                double L2 = GetL(2, point, element);
                double L3 = GetL(3, point, element);

                result += multiplier * Math.Pow(L1, psi[i][j].L1Degree)
                    * Math.Pow(L2, psi[i][j].L2Degree)
                    * Math.Pow(L3, psi[i][j].L3Degree);
            }
            return result;
        }

        public double GetL(int lNumber, Point2D point, TriangularElement element)
        {
            double x1 = ((Point2D)nodes[element.Node0]).X;
            double y1 = ((Point2D)nodes[element.Node0]).Y;
            double x2 = ((Point2D)nodes[element.Node1]).X;
            double y2 = ((Point2D)nodes[element.Node1]).Y;
            double x3 = ((Point2D)nodes[element.Node2]).X;
            double y3 = ((Point2D)nodes[element.Node2]).Y;
            double absDetD = Math.Abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

            if (lNumber == 1)
                return getS(1, 2, element, point) / absDetD;
            if (lNumber == 2)
                return getS(2, 0, element, point) / absDetD;
            if (lNumber == 3)
                return getS(0, 1, element, point) / absDetD;
            return 0;

        }

     


        public double GetSolution(Point2D point)
        {
            for (int i = 0; i < elems.Count; i++)
            {
                if (IsPointInTriangle(point, elems[i]))
                {
                    Console.WriteLine( $"Triangle: vertex0 = {elems[i].Node0} " +
                        $"vertex1 = {elems[i].Node1} " +
                        $"vertex2 = {elems[i].Node2} ");
                    return GetSolution(point, elems[i]);
                }
            }
            return 0;
        }

        bool IsPointInTriangle(Point2D point, TriangularElement element)
        {
            double x1 = ((Point2D)nodes[element.Node0]).X;
            double y1 = ((Point2D)nodes[element.Node0]).Y;
            double x2 = ((Point2D)nodes[element.Node1]).X;
            double y2 = ((Point2D)nodes[element.Node1]).Y;
            double x3 = ((Point2D)nodes[element.Node2]).X;
            double y3 = ((Point2D)nodes[element.Node2]).Y;
            double absDetD = Math.Abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

            return Math.Abs(absDetD - (getS(1, 2, element, point) + getS(2, 0, element, point) + getS(0, 1, element, point))) < 1e-7;
        }

        double getS(int i, int j, TriangularElement element, Point2D point)
        { 
            List<double> x = new List<double>();
            x.Add((nodes[element.Node0]).X);
            x.Add((nodes[element.Node1]).X);
            x.Add((nodes[element.Node2]).X);
            List<double> y = new List<double>();
            y.Add((nodes[element.Node0]).Y);
            y.Add((nodes[element.Node1]).Y);
            y.Add((nodes[element.Node2]).Y);

            return Math.Abs((x[j] - x[i]) * (point.Y - y[i]) - (point.X - x[i]) * (y[j] - y[i]));
        }
    }
}
