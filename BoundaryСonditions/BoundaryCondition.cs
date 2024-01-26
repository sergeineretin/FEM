namespace FEM.BoundaryСonditions
{
    public abstract class BoundaryCondition 
    {
        public int NodeNumber1 { get; set; }
        public int NodeNumber2 { get; set; }
        public int FormulaNumber { get; set; }

        public BoundaryCondition(int nodeNumber1, int nodeNumber2, int formulaNumber)
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
}
