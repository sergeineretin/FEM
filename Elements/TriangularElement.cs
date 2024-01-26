namespace FEM.Elements
{
    public class TriangularElement : Element, IComparable<TriangularElement>
    {
        protected int node2;
        public int Node2 { get { return node2; } set { node2 = value; } }


        protected int basisFunctionsType;

        public int CompareTo(TriangularElement other)
        {
            int minCompare = GetMinValue().CompareTo(other.GetMinValue());
            if (minCompare == 0)
            {
                int nextMinCompare = GetNextMinValue().CompareTo(other.GetNextMinValue());
                if (nextMinCompare == 0)
                {
                    return CompareAllProperties(other);
                }
                return nextMinCompare;
            }
            return minCompare;
        }

        private int GetMinValue()
        {
            return Math.Min(Node0, Math.Min(Node1, Node2));
        }

        private int GetNextMinValue()
        {
            List<int> values = new List<int> { Node0, Node1, Node2 };
            values.Remove(GetMinValue());
            return Math.Min(values[0], values[1]);
        }

        private int CompareAllProperties(TriangularElement other)
        {
            if (Node0 != other.Node0)
                return Node0.CompareTo(other.Node0);

            if (Node1 != other.Node1)
                return Node1.CompareTo(other.Node1);

            return Node2.CompareTo(other.Node2);
        }

        public override int BasisFunctionsType
        {
            get { return basisFunctionsType; }
            set 
            {
                basisFunctionsType = value;
                if ((BasisFunctions)value == BasisFunctions.Linear)
                {
                    globalIndexesOfUnknowns = new List<int>(new int[3]);
                }
                else if ((BasisFunctions)value == BasisFunctions.Quadratic)
                {
                    globalIndexesOfUnknowns = new List<int>(new int[6]);
                }
                else if ((BasisFunctions)value == BasisFunctions.CubicLagrangian)
                {
                    globalIndexesOfUnknowns = new List<int>(new int[10]);
                }
                else if ((BasisFunctions)value == BasisFunctions.CubicHermite)
                {
                    globalIndexesOfUnknowns = new List<int>(new int[10]);
                }
                else throw new NotImplementedException();
            }
        }

        public TriangularElement(int node0, int node1, int node2, int subarea) : base(subarea, node0, node1)
        {
            this.node2 = node2;
            BasisFunctionsType = (int)BasisFunctions.Linear;
        }
        public override int GetNumberOfUnknowns() {
            if ((BasisFunctions)BasisFunctionsType == BasisFunctions.Linear)
                return 3;
            else if ((BasisFunctions)BasisFunctionsType == BasisFunctions.Quadratic)
                return 6;
            else if (((BasisFunctions)BasisFunctionsType == BasisFunctions.CubicHermite || (BasisFunctions)BasisFunctionsType == BasisFunctions.CubicLagrangian))
                return 10;
            else throw new NotImplementedException();
        }
        public override int GetNumberOfVertexes()
        {
            return 3;
        }
        public override int GetIndexOfUnknown(int localNumber)
        {
            return globalIndexesOfUnknowns[localNumber];
        }

        public override int GetIndexOfVertex(int localNumber)
        {
            if (localNumber == 0)
                return node0;
            if (localNumber == 1)
                return node1;
            if (localNumber == 2)
                return node2;
            else throw new ArgumentException("The local vertex number for a triangular element can be 0, 1, 2.");
        }

        public override string ToString()
        {
            return $"{Node0} {Node1} {Node2} {Subarea}";
        }
    }
}
