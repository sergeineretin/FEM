using FEM.FiniteElementsMethods;
using FEM.Problems;
using FEM.Solvers;

namespace FEM
{
    internal class Program
    {
        static void Main(string[] args)
        {
            

            CubicLagrangianFiniteElementsMethod fem = new CubicLagrangianFiniteElementsMethod();
            fem.GenerateGrid("sreda.txt", "bound.txt");

            Problem problem = new TestBoundaryValueProblem();
            fem.GenSLAE(problem);
            fem.WriteKuslau(1e-18, 1000);

            CGMIncompleteLLTSolver solver = new CGMIncompleteLLTSolver();
            solver.ShowProgress = true;
            solver.Solve();

            fem.Q = solver.GetSolution();
            fem.PrintSolution();
            Console.WriteLine(  $"(4.5,2.5556) = {fem.GetSolution(new(4.5,2.5556))}");
        }
    }
}
