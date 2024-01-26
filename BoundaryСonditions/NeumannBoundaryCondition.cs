namespace FEM.BoundaryСonditions
{
    public class NeumannBoundaryCondition : BoundaryCondition
    {
        public NeumannBoundaryCondition(int nodeNumber1, int nodeNumber2, int formulaNumber) : base(nodeNumber1, nodeNumber2, formulaNumber)
        {
        }
    }
}
