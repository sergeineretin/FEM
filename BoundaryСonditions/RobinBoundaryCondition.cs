namespace FEM.BoundaryСonditions
{
    public class RobinBoundaryCondition : BoundaryCondition
    {
        public RobinBoundaryCondition(int nodeNumber1, int nodeNumber2, int formulaNumber) : base(nodeNumber1, nodeNumber2, formulaNumber)
        {
        }
    }
}
