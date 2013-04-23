using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Gurobi;


namespace Metabolism.Analysis.Calculators
{

    [Serializable()]
    public class CalculationException : System.Exception
    {
        public CalculationException() : base() { }
        public CalculationException(string message) : base(message) { }
        public CalculationException(string message, System.Exception inner) : base(message, inner) { }

        // A constructor is needed for serialization when an
        // exception propagates from a remoting server to the client. 
        protected CalculationException(System.Runtime.Serialization.SerializationInfo info,
            System.Runtime.Serialization.StreamingContext context) { }
    }




    public interface IIntFCACalculator
    {
        bool[] CalculateMax(List<int> R, int n, ICollection<bool[]> witnesses, out LinkedList<bool[]> new_witnesses, out bool[] calcs, bool[] MAX, IIntCoupling C = null);
        bool[] CalculateMax(bool[] lower_bound, bool[] upper_bound, IIntCoupling coupling, out LinkedList<bool[]> witnesses, out bool[] calcs);

        int SolverCalls { get; }
    }

    public abstract class FCACalculator : IIntFCACalculator
    {

        public readonly int m, n;
        protected readonly bool[] rev;


        protected GRBModel model;
        protected GRBVar[] x;


        public readonly double TOLERANCE, MAXVALUE;

        protected int solver_calls;

        protected FCACalculator(double[,] matrix, bool[] rev, double max_value, double tolerance)
        {
            explain_bounds = false;
            this.solver_calls = 0;
            this.rev = rev;

            this.n = rev.Length;
            this.m = matrix.Length / n;

            this.TOLERANCE = tolerance;
            this.MAXVALUE = max_value;

            // Create variables
            GRBEnv env = new GRBEnv();
            model = new GRBModel(env);
            model.GetEnv().Set(GRB.IntParam.OutputFlag, 0);
            model.GetEnv().Set(GRB.DoubleParam.FeasibilityTol, 1e-8);

            double[] lb, ub, obj;
            char[] type;

            lb = new double[n];
            ub = new double[n];
            obj = new double[n];
            type = new char[n];

            for (int i = 0; i < n; ++i)
            {
                lb[i] = rev[i] ? -GRB.INFINITY : 0;
                ub[i] = GRB.INFINITY;
                obj[i] = 0;
                type[i] = GRB.CONTINUOUS;
            }

            x = model.AddVars(lb, ub, obj, type, null);
            model.Update();


            // flux conservation
            GRBLinExpr expr;
            for (int i = 0; i < m; ++i)
            {
                expr = new GRBLinExpr();


                for (int j = 0; j < n; ++j)
                    if (Math.Abs(matrix[i, j]) > 1e-6)
                        expr.AddTerm(matrix[i, j], x[j]);

                model.AddConstr(expr, GRB.EQUAL, 0, null);
            }
            model.Update();
        }

        public bool[] CalculateMax(List<int> R, int n, ICollection<bool[]> witnesses, out LinkedList<bool[]> new_witnesses, out bool[] calcs, bool[] MAX = null, IIntCoupling C = null)
        {

            bool[] lb = new bool[2 * n];
            bool keep;
            foreach (bool[] w in witnesses)
            {
                keep = true;
                for (int i = 0; keep && i < R.Count; ++i)
                    keep = !w[R[i]];
                if (keep)
                    AddWitnesses(lb, w);
            }

            bool[] ub = new bool[n];
            if (MAX == null)
                for (int r = 0; r < n; ++r)
                    ub[r] = true;
            else
                for (int r = 0; r < n; ++r)
                    ub[r] = MAX[r];

            ICollection<int> c;
            if (C == null)
                foreach (int r in R)
                    ub[r] = false;
            else
                foreach (int r in R)
                {
                    c = C[r];
                    foreach (int i in c)
                        ub[i] = false;
                }


            return CalculateMax(lb, ub, C, out new_witnesses, out calcs);
        }

        public abstract bool[] CalculateMax(bool[] lb, bool[] ub, IIntCoupling C, out LinkedList<bool[]> new_witnesses, out bool[] calcs);

        public void AddWitnesses(bool[] w1, bool[] w2)
        {
            if (w1.Length != 2 * n || w2.Length != 2 * n)
                return;

            for (int r = 0; r < 2*n; ++r)
                w1[r] |= w2[r];

        }


        private bool explain_bounds;
        protected bool ExplainBounds { get { return explain_bounds; } set { explain_bounds = value; } }

        public bool[] ExtractWitness()
        {

            bool[] new_witness = new bool[2 * n];

            for (int j = 0; j < n; ++j)
            {
                new_witness[j + n] = Math.Abs(x[j].Get(GRB.DoubleAttr.X)) > TOLERANCE;
                new_witness[j] = Math.Abs(x[j].Get(GRB.DoubleAttr.X)) > TOLERANCE / MAXVALUE;
            }

            /*int r = 52, s = 66;
            if (!new_witness[r] && new_witness[s + n])
                Console.WriteLine("We found a feasible solution with x_{0} = {2:0.00000} and x_{1} = {3:0.000}.", r, s, x[r].Get(GRB.DoubleAttr.X), x[s].Get(GRB.DoubleAttr.X));
            */

            if (ExplainBounds)
            {
                bool changed = false;
                for (int r = 0; r < n; ++r)
                    if (
                        (r == 213 && new_witness[r]) || 
                        (!new_witness[r] && Math.Abs(x[r].Get(GRB.DoubleAttr.X)) > (TOLERANCE/MAXVALUE)/10)
                        )
                    {
                        Console.WriteLine("x_{0}={1:0.####################} -> ({2}, {3}).", r, x[r].Get(GRB.DoubleAttr.X), new_witness[n + r] ? 1 : 0, Math.Abs(x[r].Get(GRB.DoubleAttr.X)) > TOLERANCE / MAXVALUE ? 1 : 0);
                        changed = true;
                    }
                if (changed)
                {
                    Console.WriteLine("We used tolerance values of {0:0.##e-0###} and {1:0.##e-0####}.", TOLERANCE, TOLERANCE / MAXVALUE);
                    Console.WriteLine();
                }
            }
            return new_witness;
        }

        public int SolverCalls
        {
            get { return this.solver_calls; }
        }
    }



    public class FastFCIntFCACalculator : FCACalculator
    {
        private LinkedList<bool[]> circles;

        protected GRBVar[] supp;
        
        private LPIntFCACalculator relaxation_solver;

        public bool WitnessContainsCircle(bool[] witness)
        {

            foreach (bool[] c in circles)
            {
                bool contained = true;
                for (int i = 0; i < n && contained; ++i)
                    if (c[i] && !witness[i])
                        contained = false;
                if (contained)
                    return true;
            }
            return false;
        }

        private void ForbidCycles()
        {
            int k;
            GRBLinExpr expr;
            foreach (bool[] cycle in this.circles)
            {

                k = 0;
                expr = new GRBLinExpr();
                for (int i = 0; i < n; ++i)
                {
                    if (cycle[i])
                    {
                        expr.AddTerm(1, supp[i]);
                        ++k;
                    }
                }
                model.AddConstr(expr, GRB.LESS_EQUAL, k - 1, null);
            }
            model.Update();
        }

        public FastFCIntFCACalculator(double[,] matrix, bool[] rev, ICollection<bool[]> circles, double max_value, double tolerance) :
            base(matrix, rev, max_value, tolerance)
        {
            ExplainBounds = false;

            // create circle specific support variables
            this.relaxation_solver = new LPIntFCACalculator(matrix, rev, max_value, tolerance);
            this.circles = new LinkedList<bool[]>(circles);

            double[] lb, ub, obj;
            char[] type;

            lb = new double[n];
            ub = new double[n];
            obj = new double[n];
            type = new char[n];

            for (int i = 0; i < n; ++i)
            {
                lb[i] = 0;
                ub[i] = 1;
                obj[i] = 0;
                type[i] = GRB.BINARY;
            }
            supp = model.AddVars(lb, ub, obj, type, null);
            model.Update();


            // supports fuer a_i = 0 --> v_i = 0
            GRBLinExpr expr;
            for (int i = 0; i < n; ++i)
            {
                expr = new GRBLinExpr();
                expr.AddTerm(max_value, supp[i]);
                model.AddConstr(expr, GRB.GREATER_EQUAL, x[i], null);

                if (rev[i])
                {
                    expr = new GRBLinExpr();
                    expr.AddTerm(-max_value, supp[i]);
                    model.AddConstr(expr, GRB.LESS_EQUAL, x[i], null);
                }
            }


            model.Update();

            ForbidCycles();
        }


        public override bool[] CalculateMax(bool[] lower_bound, bool[] upper_bound, IIntCoupling coupling, out LinkedList<bool[]> witnesses, out bool[] calcs)
        {


            LinkedList<bool[]> wits_relaxed;
            bool[] ub = new bool[upper_bound.Length];
            Array.Copy(upper_bound, ub, ub.Length);

            bool[] lb = new bool[lower_bound.Length];
            Array.Copy(lower_bound, lb, lb.Length);

            upper_bound = relaxation_solver.CalculateMax(lb, ub, coupling, out wits_relaxed, out calcs);


            witnesses = new LinkedList<bool[]>();
            foreach (bool[] w in wits_relaxed)
                if (!this.WitnessContainsCircle(w))
                {
                    AddWitnesses(lower_bound, w);
                    witnesses.AddLast(w);
                }
                else
                {
                    //Console.Out.WriteLine("Found a pathway containing a circle.");
                }


            bool allrev = true;
            for (int i = 0; allrev && i < n; ++i)
                allrev = rev[i] || !upper_bound[i];


            for (int i = 0; i < n; ++i)
            {
                if (!upper_bound[i])
                    supp[i].Set(GRB.DoubleAttr.UB, 0);
            }
            model.Update();

            calcs = new bool[n];
            for (int i = 0; i < n; ++i)
                if (calcs[i] = upper_bound[i] && !lower_bound[i+n])
                {
                    //Console.Out.WriteLine("Solving MILP for r = {0}", i);

                    x[i].Set(GRB.DoubleAttr.LB, 1);
                    model.Update();
                    model.Optimize();
                    this.solver_calls++;

                    if (rev[i] && !allrev && model.Get(GRB.IntAttr.Status) != GRB.Status.OPTIMAL && model.Get(GRB.IntAttr.Status) != GRB.Status.SUBOPTIMAL)
                    {
                        x[i].Set(GRB.DoubleAttr.LB, -GRB.INFINITY);
                        x[i].Set(GRB.DoubleAttr.UB, -1);
                        model.Update();
                        model.Optimize();
                        this.solver_calls++;
                    }

                    switch (model.Get(GRB.IntAttr.Status))
                    {
                        case GRB.Status.INFEASIBLE:
                        case GRB.Status.INF_OR_UNBD:
                            {
                                if (coupling == null)
                                {
                                    upper_bound[i] = false;
                                    supp[i].Set(GRB.DoubleAttr.UB, 0);
                                }
                                else
                                {
                                    ICollection<int> c = coupling[i];
                                    foreach (int r in c)
                                    {
                                        upper_bound[r] = false;
                                        supp[r].Set(GRB.DoubleAttr.UB, 0);
                                    }
                                }
                                break;
                            }
                        case GRB.Status.OPTIMAL:
                        case GRB.Status.SUBOPTIMAL:
                            {
                                bool[] new_witness = ExtractWitness();

                                AddWitnesses(lower_bound, new_witness);
                                witnesses.AddLast(new_witness);
                                break;
                            }
                        default:
                            {
                                break;
                            }
                    }
                    x[i].Set(GRB.DoubleAttr.UB, GRB.INFINITY);
                    x[i].Set(GRB.DoubleAttr.LB, rev[i] ? -GRB.INFINITY : 0);
                    model.Update();

                }

            for (int i = 0; i < n; ++i)
            {
                supp[i].Set(GRB.DoubleAttr.UB, 1);
                if (upper_bound[i] && !lower_bound[i + n])
                    Console.WriteLine("Difference for reaction {0}.", i);
            }
            model.Update();


            return upper_bound;
        }
        
    }

    public class LPIntFCACalculator : FCACalculator
    {

        public LPIntFCACalculator(double[,] matrix, bool[] rev, double max_value, double tolerance): base(matrix, rev, max_value, tolerance)
        {
            ExplainBounds = false;
        }



        public override bool[] CalculateMax(bool[] lower_bound, bool[] upper_bound, IIntCoupling coupling, out LinkedList<bool[]> witnesses, out bool[] calcs)
        {
            bool allrev = true;
            for (int i = 0; allrev && i < n; ++i)
                allrev = rev[i] || !upper_bound[i];


            for (int i = 0; i < n; ++i)
            {
                if (!upper_bound[i])
                {
                    x[i].Set(GRB.DoubleAttr.UB, 0);
                    if (rev[i])
                        x[i].Set(GRB.DoubleAttr.LB, 0);
                }
            }
            model.Update();


            calcs = new bool[n];
            witnesses = new LinkedList<bool[]>();
            for (int i = 0; i < n; ++i)
                if (calcs[i] = upper_bound[i] && !lower_bound[i+n])
                {
                    x[i].Set(GRB.DoubleAttr.LB, 1);
                    model.Update();
                    model.Optimize();
                    this.solver_calls++;

                    if (rev[i] && !allrev && model.Get(GRB.IntAttr.Status) != GRB.Status.OPTIMAL && model.Get(GRB.IntAttr.Status) != GRB.Status.SUBOPTIMAL)
                    {
                        x[i].Set(GRB.DoubleAttr.LB, -GRB.INFINITY);
                        x[i].Set(GRB.DoubleAttr.UB, -1);
                        model.Update();
                        model.Optimize();
                        this.solver_calls++;
                    }

                    switch (model.Get(GRB.IntAttr.Status))
                    {
                        case GRB.Status.INFEASIBLE:
                        case GRB.Status.INF_OR_UNBD:
                            {
                                if (coupling == null)
                                {
                                    upper_bound[i] = false;
                                    x[i].Set(GRB.DoubleAttr.UB, 0);
                                    if (rev[i])
                                        x[i].Set(GRB.DoubleAttr.LB, 0);
                                }
                                else
                                {
                                    ICollection<int> c = coupling[i];
                                    foreach (int r in c)
                                    {
                                        upper_bound[r] = false;
                                        x[r].Set(GRB.DoubleAttr.UB, 0);
                                        if (rev[r])
                                            x[r].Set(GRB.DoubleAttr.LB, 0);
                                    }
                                }
                                break;
                            }
                        case GRB.Status.OPTIMAL:
                        case GRB.Status.SUBOPTIMAL:
                            {
                                bool[] new_witness = ExtractWitness();

                                AddWitnesses(lower_bound, new_witness);
                                witnesses.AddLast(new_witness);
                                break;
                            }
                        default: { break; }
                    }

                    x[i].Set(GRB.DoubleAttr.UB, GRB.INFINITY);
                    x[i].Set(GRB.DoubleAttr.LB, rev[i] ? -GRB.INFINITY : 0);
                    model.Update();

                }

            for (int i = 0; i < n; ++i)
            {
                x[i].Set(GRB.DoubleAttr.UB, GRB.INFINITY);
                if (rev[i])
                    x[i].Set(GRB.DoubleAttr.LB, -GRB.INFINITY);

                if (upper_bound[i] && !lower_bound[i + n])
                    Console.WriteLine("Difference for reaction {0}.", i);
            }
            model.Update();

            return upper_bound;
        }
        

    }
}