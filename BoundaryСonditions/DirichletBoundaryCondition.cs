namespace FEM.BoundaryСonditions
{
    public class DirichletBoundaryCondition : BoundaryCondition
    {
        public DirichletBoundaryCondition(int nodeNumber1, int nodeNumber2, int formulaNumber) : base(nodeNumber1, nodeNumber2, formulaNumber)
        {
        }
    }
}
