namespace FEM.Elements
{
    public class Element
    {
        protected int subarea = 0;
        protected int node0 = 0;
        protected int node1 = 0;

        protected int edgeNodesNumber;
        protected int innerNodesNumber;
        protected List<int> globalIndexesOfUnknowns;
        bool isGlobalIndexesOfUnknownsSet = false;
        public int VertexNodesNumber {
            get; init;
        }

        public enum BasisFunctions
        {
            Linear,
            Quadratic,
            CubicLagrangian,
            CubicHermite
        }
        public virtual int BasisFunctionsType
        {
            get;
            set;
        }
        public Element(int subarea, int node0, int node1)
        {
            this.subarea = subarea;
            this.node0 = node0;
            this.node1 = node1;
            BasisFunctionsType = (int)BasisFunctions.Linear;
        }
        public Element(BasisFunctions basisFunctions, int subarea, int node0, int node1)
        {
            this.subarea = subarea;
            this.node0 = node0;
            this.node1 = node1;
            BasisFunctionsType = (int)BasisFunctions.Linear;
        }
        public int Node0 { get { return node0; } set { node0 = value; } }
        public int Node1 { get { return node1; } set { node1 = value; } }
        public int Subarea { get { return subarea; } }
        public virtual int GetNumberOfUnknowns() { return 2; }
        public virtual int GetNumberOfVertexes() { return 2; }
        public virtual int GetIndexOfUnknown(int localNumber)
        {
            if (localNumber == 0)
                return node0;
            else if (localNumber == 1)
                return node1;
            else throw new ArgumentException("Local number for a one-dimensional element can be either 0 or 1");
        }

        public void SetGlobalIndicesOfUnknowns(List<int> globalIndicesOfUnknown)
        {
            try
            {
                if (!isGlobalIndexesOfUnknownsSet)
                {
                    if (globalIndexesOfUnknowns.Count != 0)
                    {
                        globalIndexesOfUnknowns = globalIndicesOfUnknown;
                        isGlobalIndexesOfUnknownsSet = true;
                    }
                    else
                        throw new ArgumentException();
                }
                else
                { 
                    throw new Exception("Attempt to reset the given indices of unknowns on the element");
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
        }

        public virtual int GetIndexOfVertex(int localNumber)
        {
            return -1;
        }

    }
}
