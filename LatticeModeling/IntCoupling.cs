/*
 *  Copyright (C) 2013 Dr. Yaron Goldstein, Freie Universitaet Berlin
 *  
 *  This work is licensed under the Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0) License. 
 *  To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ 
 *  or send a letter to Creative Commons, 171 Second Street, Suite 300, San Francisco, California, 94105, USA.
 *  
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 *  See the License for more details. 
 *  You should have received a copy of the License along with this program. 
 *  If not, see <http://creativecommons.org/licenses/by-nc-sa/3.0/>.
 */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Metabolism.Analysis.Calculators;

namespace Metabolism.Analysis
{

    public interface IIntCoupling
    {
        ICollection<int> this[int index] { get; }

        int Length { get; }

        LinkedList<bool[]> Witnesses { get; }

        bool[] Max(int index);

        LinkedList<String> ToDot(String[] names);



        DateTime[] TimeReaction { get; }
        int[] LPsReaction { get; }
        int[] WitsReaction { get; }
        int WitsUsed { get; }
    }


    public class IntCoupling : IIntCoupling
    {

        private readonly int n;
        private readonly LinkedList<int>[] coupled;
        private readonly LinkedList<bool[]> witnesses;
        private readonly bool[][] maxima;
        private readonly bool[] max;

        private readonly DateTime[] time_reaction;
        private readonly int[] lps_reaction;
        private readonly int[] wits_reaction; 
        private readonly int wits_used;

        public DateTime[] TimeReaction { get { return time_reaction; } }
        public int[] LPsReaction { get { return lps_reaction; } }
        public int[] WitsReaction { get { return wits_reaction; } }
        public int WitsUsed { get { return wits_used; } }
        public LinkedList<bool[]> Witnesses { get { 
            return new LinkedList<bool[]>(witnesses); } }

        public IntCoupling(IIntFCACalculator calculator, bool show_output, ICollection<bool[]> witnesses,bool[] max, List<int> R = null, IIntCoupling coupling = null, int start = 0)
        {
            this.max = max;
            this.n = max.Length;

            time_reaction = new DateTime[n];
            lps_reaction = new int[n];
            wits_reaction = new int[n];
            int lps_count = calculator.SolverCalls;


            this.maxima = new bool[n][];
            coupled = new LinkedList<int>[n];
            for (int i = 0; i < n; ++i)
            {
                coupled[i] = new LinkedList<int>();
                coupled[i].AddLast(i);
            }

            if (R == null)
                R = new List<int>();


            LinkedList<bool[]> new_witnesses = new LinkedList<bool[]>(witnesses);
            this.witnesses = new LinkedList<bool[]>();
            bool allowed;
            foreach (bool[] w in new_witnesses)
            {
                allowed = true;
                for (int i = 0; allowed && i < R.Count; ++i)
                    allowed = !w[R[i]];
                if (allowed)
                    this.witnesses.AddLast(w);
            }
            wits_used = this.witnesses.Count;


            if (coupling == null)
                coupling = this;

            bool[] opt, calcs;
            R.Add(-1);
            ICollection<int> c;
            for (int i = start; i < n; ++i)
                if (maxima[i] == null)
                {
                    if (max[i])
                    {


                        R[R.Count - 1] = i;

                        opt = calculator.CalculateMax(R, n, this.witnesses, out new_witnesses, out calcs, max, coupling);
                        lps_reaction[i] = calculator.SolverCalls - lps_count;
                        lps_count = calculator.SolverCalls;

                        c = coupling[i];
                        foreach (int r in c)
                            if (r >= i && coupling.Max(r) != null && !coupling.Max(r)[i])
                            {
                                maxima[r] = opt;
                                for (int j = 0; j < n; ++j)
                                    if (max[j] && !maxima[r][j] && r != j)
                                        coupled[r].AddLast(j);
                            }

                        if (maxima[i] == null)
                        {
                            maxima[i] = opt;
                            if (maxima[i] == null)
                                Console.WriteLine("Fehler!");
                            else
                                for (int j = 0; j < n; ++j)
                                    if (!maxima[i][j] && max[j] && i != j)
                                        coupled[i].AddLast(j);
                        }

                        foreach (bool[] a in new_witnesses)
                                {
                                    this.witnesses.AddLast(a);
                                    wits_reaction[i]++;
                                }

                        time_reaction[i] = DateTime.Now;
                        if (show_output)
                            Console.WriteLine("({0})\tCoupling {1} from {2} calculated. ({3} LPs solved, {4} kept.)\n", time_reaction[i], i + 1, n, lps_reaction[i], this.witnesses.Count);

                    }
                    else
                        maxima[i] = max;
                }
            R.RemoveAt(R.Count - 1);
        }

        public ICollection<int> this[int index]
        {
            get { return coupled[index]; }
        }

        public int Length
        {
            get { return n; }
        }


        public bool[] Max(int index)
        {
            return maxima[index];
        }

        public LinkedList<string> ToDot(String[] reaction_names)
        {
            LinkedList<String> res = new LinkedList<string>();

            res.AddFirst("digraph FCA {");
            res.AddLast("\tnodesep=0.7;");
            res.AddLast("\trankdir=LR;");

            if (reaction_names != null)
                for (int i = 0; i < n; ++i)
                    res.AddLast(String.Format("\t{0}[label=\"{1}\"];", i + 1, reaction_names[i]));

            // compress
            int[] component = new int[n];
            List<int>[] components = new List<int>[n];
            for (int i = 0; i < n; ++i)
            {
                component[i] = -1;
                components[i] = new List<int>();
            }

            for (int i = 0; i < n; ++i)
                if (component[i] < 0)
                {
                    component[i] = i;
                    components[i].Add(i);
                    for (int j = i + 1; j < n; ++j)
                        if (!Max(i)[j] && !Max(j)[i])
                        {
                            component[j] = i;
                            components[i].Add(j);
                        }
                }

            // add components
            for (int i = 0; i < n; ++i)
                if (component[i] == i)
                    switch (components[i].Count)
                    {
                        case 0:
                        case 1:
                            break;
                        case 2:
                            res.AddLast(String.Format("\t{0}->{1}[dir=both];", components[i][0] + 1, components[i][1] + 1));
                            break;
                        default:
                            for (int j = 0; j < components[i].Count; ++j)
                                res.AddLast(String.Format("\t{0}->{1}[dir=both];", components[i][j] + 1, components[i][(j + 1) % components[i].Count] + 1));
                            break;
                    }

            // connect components
            for (int i = 0; i < n; ++i)
                if (component[i] == i)
                    for (int j = 0; j < n; ++j)
                        if (!Max(i)[j] && component[j] == j && i != j)
                            res.AddLast(String.Format("\t{0}->{1};", i + 1, j + 1));

            res.AddLast("}");

            return res;
        }
    }


    public interface IEFCACount
    {

        int this[int i, int j] { get; }

        int Length { get; }

        int WitnessCount { get; }

        int LPCount { get; }

        LinkedList<string> ToCSV(string separator, string[] names = null);



        DateTime[] TimeReaction { get; }
        int[] LPsReaction { get; }
        int[] WitsReaction { get; }
        int[] WitsUsable { get; }
    }

    public class EFCACount : IEFCACount
    {

        private readonly int[,] max_size, targeted_size;
        private readonly bool[] max;
        private readonly List<int> T; 
        public readonly int n;
        private int lp_counter = 0, witness_counter = 0;
        private readonly LinkedList<bool[]> witnesses_kept;
        private readonly Random rand = new Random();

        private readonly DateTime[] time_reaction; 
        private readonly int[] lps_reaction;
        private readonly int[] wits_reaction;
        private readonly int[] wits_usable;

        public DateTime[] TimeReaction { get { return time_reaction; } }
        public int[] LPsReaction { get { return lps_reaction; } }
        public int[] WitsReaction { get { return wits_reaction; } }
        public int[] WitsUsable { get { return wits_usable; } }

        public EFCACount(IIntFCACalculator calculator, ICollection<bool[]> witnesses, bool[] max, IIntCoupling coupling = null, ICollection<int> T = null, bool show_output = true)
        {
            this.max = max;
            this.n = max.Length;

            time_reaction = new DateTime[n];
            lps_reaction = new int[n];
            for (int i = 0; i < n; ++i)
                lps_reaction[i] = -1;
            wits_reaction = new int[n];
            wits_usable = new int[n];

            int temp_counter = calculator.SolverCalls;
            if (coupling == null)
            {
                coupling = new IntCoupling(calculator, false, witnesses, max);
                this.witness_counter = coupling.Witnesses.Count;
                this.lp_counter = calculator.SolverCalls - temp_counter;
                temp_counter = calculator.SolverCalls;
            }

            if (T == null)
            {
                this.T = new List<int>(n);
                for (int i = 0; i < n; ++i)
                    this.T.Add(i);
            }
            else
                this.T = new List<int>(T);

            LinkedList<bool[]> new_witnesses = new LinkedList<bool[]>();
            witnesses_kept = new LinkedList<bool[]>(witnesses);
            IIntCoupling row;

            max_size = new int[n, n];
            targeted_size = new int[n, n];
            for (int i = 0; i < n; ++i)
            {
                max_size[i, i] = -1;
                targeted_size[i, i] = -1;
            }

            List<int> R = new List<int>();
            ICollection<int> c;
            R.Add(-1);
            for (int i = 0; i < n; ++i)
                if (max[i] && max_size[i, i] < 0)
                {
                    R[0] = i;


                    row = new IntCoupling(calculator, false, witnesses_kept, max, R, coupling, i);
                    wits_usable[i] = row.WitsUsed;
                    if (show_output)
                        Console.WriteLine("({0})\tRow {1} from {2} calculated of EFCA. ({3} LPs solved, {4} kept.)\n", System.DateTime.Now, i + 1, n, calculator.SolverCalls - temp_counter, witnesses_kept.Count);

                    this.lp_counter += calculator.SolverCalls - temp_counter;
                    temp_counter = calculator.SolverCalls;
                    new_witnesses = row.Witnesses;
                    this.witness_counter += row.Witnesses.Count;


                    c = coupling[i];
                    foreach (int r in c)
                        if (r >= i && !coupling.Max(r)[i])
                            for (int j = r; j < n; ++j)
                                if (max[j])
                                {
                                    foreach(bool unblocked in row.Max(j))
                                        if(unblocked)
                                            max_size[r, j]++;
                                    max_size[j, r] = max_size[r, j];

                                    foreach(int t in this.T)
                                        if(row.Max(j)[t])
                                            targeted_size[r, j] ++;
                                    targeted_size[j, r] = targeted_size[r, j];
                                }

                    lps_reaction[i] = lp_counter;
                    wits_reaction[i] = witness_counter;
                    time_reaction[i] = DateTime.Now;
                }

            if (show_output)
                Console.WriteLine("({0})\tEFCA uncompressed.\n", System.DateTime.Now);
        }

        public int this[int i, int j]
        {
            get { return max_size[i, j]; }
        }

        public int Length
        {
            get { return n; }
        }

        public int WitnessCount
        {
            get { return witness_counter; }
        }

        public int LPCount
        {
            get { return lp_counter; }
        }


        public LinkedList<string> ToCSV(string separator, string[] names = null)
        {

            if (names == null)
            {
                names = new string[n];
                for (int i = 0; i < n; ++i)
                    names[i] = (i + 1).ToString();
            }

            LinkedList<string> res = new LinkedList<string>();
            res.AddFirst(String.Format("i{0}j{0}Number.of.maximum{0}Number.of.targeted", separator));

            for (int i = 0; i < n; ++i)
                if (this.max[i])
                    for (int j = 0; j < n; ++j)
                        if (this.max[j])
                            res.AddLast(String.Format("{1}{0}{2}{0}{3}{0}{4}", separator, names[i], names[j], this[i, j], this.targeted_size[i, j]));

            return res;
        }
    }
}