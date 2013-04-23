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
using Metabolism.Analysis;
using Metabolism.Analysis.Calculators;
using Tools.IO;
using Tools.Converter.Biology;

namespace FCA
{
    public class FCA
    {


        public static void Main(string[] args)
        {

            //init
            string path2models = args[0], name, separator = args[1];

            bool use_names = true;
            bool do_efca = args.Length > 5 && int.Parse(args[5]) == 1, no_output = args.Length > 4 && int.Parse(args[4]) == 0;

            int n_simulations = args.Length > 6 ? int.Parse(args[6]) : 1;
            double max_value = double.Parse(args[2]), tolerance = double.Parse(args[3]);
            LinkedList<bool[]> forbidden = new LinkedList<bool[]>();


            double[,] matrix; bool[] rev;
            for (int i = 0; i < 1; ++i)
            {
                bool fc = i  == 1;


                LinkedList<string> paths = Reader.ReadFile(path2models), res = new LinkedList<string>(), details = new LinkedList<string>();
                res.AddFirst("\\begin{table}");
                res.AddLast("\\centering");
                res.AddLast("\\begin{tabular}{llrrrr}");
                res.AddLast("\\toprule");
                res.AddLast("\\rowcolor{gray!10}\\textbf{Model}&\\textbf{Value}&\\textbf{Number}&\\textbf{LPs}&\\textbf{Witnesses}&\\textbf{Time}\\\\");
                //res.AddLast("Model Method Metabolites Reactions Reactions.unblocked Reactions.frev Couples Time Time.unblocked Time.frev Time.couples Time.efca LPs LPs.unblocked LPs.frev LPs.couples LPs.efca Witnesses Witnesses.unblocked Witnesses.frev Witnesses.couples Witnesses.efca");
                //details.AddLast("Model Method Metabolites Reactions Reactions.unblocked Reactions.frev Couples Reaction Time LPs LPs.total Witnesses Witnesses.total Witnesses.distribution Witnesses.usable");
                String model, method;
                foreach (string path in paths)
                {

                    // read network
                    forbidden.Clear();

                    String[] metabolite_names = null, reaction_names = null;
                    if (path.Contains("matrev"))
                    {
                        List<double>[] network = Reader.ReadDoubleFileToList(path, separator);
                        NetworkReader.ExtractNetwork(network, out matrix, out rev);
                        name = path.Substring(0, path.IndexOf(".matrev"));
                    }
                    else if (path.Contains("meta"))
                    {
                        NetworkReader.ReadMetatool(Reader.ReadFile(path), out matrix, out rev, out metabolite_names, out reaction_names);
                        name = path.Substring(0, path.IndexOf(".meta"));
                    }
                    else
                    {
                        name = path;
                        NetworkReader.ReadDir(name, separator, out matrix, out rev, out metabolite_names, out reaction_names, out forbidden);
                    }
                    Writer.WriteFile(NetworkReader.NetworkToDot(matrix, rev, metabolite_names, reaction_names), name + "_network.dot");
                    Writer.WriteFile(NetworkReader.NetworkToConnections(matrix, rev, metabolite_names, reaction_names), name + "_connections.dot");
                    Writer.WriteFile(NetworkReader.NetworkToMatlab(matrix, rev), name + ".m");

                    model = name.Substring(name.LastIndexOf('\\') + 1);

                    //bool[] trans = FindTrans(matrix, rev);

                    // Calculations.

                    LinkedList<String> lines = new LinkedList<string>();
                    {
                        FCA fca = null;
                        method = fc ? "fc" : "all";
                        for (int t = 0; t < n_simulations; ++t)
                        {


                            try
                            {
                                fca = new FCA(matrix, rev, forbidden, max_value, tolerance, fc, do_efca, !no_output);


                                fca.Results(model, method, res, details);
                                //Writer.WriteFile(details, String.Format("{0}_details_{1}_{2}-{3}-{4}-{5}-{6}.log", path2models.Substring(0, path2models.Length - 4), fc ? "fc" : "all", DateTime.Now.Year, DateTime.Now.Month, DateTime.Now.Day, DateTime.Now.Hour, DateTime.Now.Minute));
                                Writer.WriteFile(res, String.Format("{0}_{1}_{2}-{3}-{4}-{5}-{6}.log", path2models.Substring(0, path2models.Length - 4), fc ? "fc" : "all", DateTime.Now.Year, DateTime.Now.Month, DateTime.Now.Day, DateTime.Now.Hour, DateTime.Now.Minute));
                            }
                            catch (CalculationException e)
                            {
                                Console.WriteLine("There's happened an error during the calculation. Repeating the simulation.");
                                --t;
                            }
                        }

                        // save results of last simulation

                        String suffix = (fc ? "_fc" : "_all");

                        lines.Clear();
                        bool[] max = fca.max;
                        int n = fca.n;
                        for (int r = 0; r < n; ++r)
                            if (!max[r])
                                if (use_names)
                                    lines.AddLast(reaction_names[r]);
                                else
                                    lines.AddLast(r.ToString());
                        Writer.WriteFile(lines, name + "\\blocked" + suffix + ".txt");

                        lines.Clear();
                        for (int r = 0; r < n; ++r)
                            if (max[r])
                                for (int s = 0; s < n; ++s)
                                    if (max[s] && !fca.coupling.Max(r)[s] && r != s)
                                        if (use_names)
                                            lines.AddLast(String.Format("{0} -> {1}", reaction_names[r], reaction_names[s]));
                                        else
                                            lines.AddLast(String.Format("{0} -> {1}", r, s));
                        Writer.WriteFile(lines, name + "\\coupled" + suffix + ".txt");

                        Writer.WriteFile(fca.coupling.ToDot(reaction_names), name + "\\coupled" + suffix + ".dot");

                        if (do_efca)
                            Writer.WriteFile(fca.efca.ToCSV(separator, reaction_names), name + "\\efca" + suffix + ".txt");

                        Writer.WriteFile(res, String.Format("{0}_{1}_{2}-{3}-{4}-{5}-{6}.log", path2models.Substring(0, path2models.Length - 4), fc ? "fc" : "all", DateTime.Now.Year, DateTime.Now.Month, DateTime.Now.Day, DateTime.Now.Hour, DateTime.Now.Minute));
                    }
                }

                // Finish stats.
                res.AddLast("\\bottomrule");
                res.AddLast("\\end{tabular}");
                res.AddLast("\\end{table}");
                Writer.WriteFile(res, String.Format("{0}_{1}_{2}-{3}-{4}-{5}-{6}.log", path2models.Substring(0, path2models.Length - 4), fc ? "fc" : "all", DateTime.Now.Year, DateTime.Now.Month, DateTime.Now.Day, DateTime.Now.Hour, DateTime.Now.Minute));

            }
            Console.WriteLine("Done. [press Enter]");
            Console.ReadLine();
        }
        
            private bool[] max;
            private IIntCoupling coupling; 
            private EFCACount efca; 
            private int reactions_unblocked, couples; 
            private DateTime time_start, time_unblocked, time_couples, time_efca; 
            private DateTime[] time_reaction;
            private int n;
            private int lps_unblocked, lps_couples, lps_efca; 
            private int[] lps_reaction;
            private int wits_unblocked,  wits_couples, wits_efca; 
            private int[] wits_reaction, wits_distribution, wits_usable; 


        public FCA(double[,] matrix, bool[] rev, 
            LinkedList<bool[]> forbidden_combinations, 
            double max_value, double tolerance, 
            bool fc, bool do_efca, bool show_output, 
            String name = null, LinkedList<bool[]> witnesses = null)
        {

            //initialisation
            time_start = DateTime.Now;
            efca = null;
            wits_usable = null;
            n = rev.Length;

            if (witnesses == null)
                witnesses = new LinkedList<bool[]>();
            else
                witnesses = new LinkedList<bool[]>(witnesses);
            LinkedList<bool[]> blocking_witnesses;
            IIntFCACalculator fcacalculator = !fc ? (IIntFCACalculator)new LPIntFCACalculator(matrix, rev, max_value, tolerance) :
                (IIntFCACalculator)new FastFCIntFCACalculator(matrix, rev, forbidden_combinations, max_value, tolerance);
            
            // calculate blocked reactions
            bool[] calcs;
            int wits_used;

            max = fcacalculator.CalculateMax(new List<int>(), n, witnesses, out blocking_witnesses, out calcs, null, null);
            reactions_unblocked = 0;
            foreach(bool r in max)
                if(r)
                    ++reactions_unblocked;



            foreach (bool[] a in blocking_witnesses)
                        witnesses.AddLast(a);
            if (name != null)
                SaveStatistics(witnesses, n, name + "_max.stats");
            time_unblocked = DateTime.Now;
            lps_unblocked = fcacalculator.SolverCalls;
            wits_unblocked = witnesses.Count;
            if (show_output)
                Console.WriteLine("({0})\tBlocked reactions calculated: {1} of {2} blocked.\n", time_unblocked, n - reactions_unblocked, n);


            // FCA
            Console.WriteLine("\nStarting FCA.\n");

            coupling = new IntCoupling(fcacalculator, show_output, witnesses, max);
            time_reaction = coupling.TimeReaction;
            lps_reaction = coupling.LPsReaction;
            wits_reaction = coupling.WitsReaction;
            wits_used = coupling.WitsUsed;

            if (name != null)
                SaveStatistics(coupling.Witnesses, n, name + "_fca.stats");

            time_couples = DateTime.Now;
            lps_couples = 0; wits_couples = 0;
            for (int i = 0; i < n; ++i)
            {
                lps_couples += lps_reaction[i];
                wits_couples += wits_reaction[i];
            }

            if (show_output)
                Console.WriteLine("({0})\tCouples calculated.\n", time_couples);

            wits_distribution = new int[n + 1];
            foreach (bool[] a in coupling.Witnesses)
            {
                int cardinality = 0;
                for (int i = 0; i < n; ++i)
                    if (a[i])
                        ++cardinality;
                wits_distribution[cardinality]++;
            }


            // EFCA
            lps_efca = 0;
            wits_efca = 0;
            if (do_efca)
            {
                efca = new EFCACount(fcacalculator, coupling.Witnesses, max, coupling, null, show_output);
                time_reaction = efca.TimeReaction;
                lps_reaction = efca.LPsReaction;
                wits_reaction = efca.LPsReaction;
                wits_usable = efca.WitsUsable;

                lps_efca = efca.LPCount;
                wits_efca = efca.WitnessCount;
            }
            time_efca = DateTime.Now;

            couples = -reactions_unblocked; // r <-> r
            for (int r = 0; r < n; ++r)
                if (max[r])
                    for (int s = 0; s < n; ++s)
                        if (max[s] && !coupling.Max(r)[s])
                            ++couples;
        }


        public LinkedList<bool[]> Witnesses { get { return coupling.Witnesses; } } 

        public static void SaveStatistics(ICollection<bool[]> witnesses, int n, String path)
        {
            int[,] coupling = new int[n, n];
            int[] reactions = new int[n];

            foreach (bool[] a in witnesses)
                for (int i = 0; i < n; ++i)
                    if (a[i])
                    {
                        ++reactions[i];
                        for (int j = 0; j < n; ++j)
                            if (!a[j])
                                coupling[i, j]++;
                    }

            int k = witnesses.Count;

            LinkedList<String> file = new LinkedList<string>();
            file.AddLast(String.Format("{0} {1} {2}", -1, -1, k));

            for (int i = 0; i < n; ++i)
            {
                file.AddLast(String.Format("{0} {1} {2}", i + 1, -1, k - reactions[i]));
                file.AddLast(String.Format("{0} {1} {2}", -1, i + 1, reactions[i]));

                for (int j = 0; j < n; ++j)
                {
                    file.AddLast(String.Format("{0} {1} {2}", j + 1, i + 1, coupling[i, j]));
                }
            }

            Writer.WriteFile(file, path);
        }

        public void Results(
            String model, String method, 
            LinkedList<string> res, LinkedList<string> details)
        {
            
            int lps, wits;

            lps = lps_unblocked + lps_couples + lps_efca;
            wits = wits_unblocked + wits_couples + wits_efca;

            String cellcolor = "\\cellcolor{gray!10}";
            res.AddLast("\\midrule\n\\multirowbt{4}{*}{\\Var{" + model.Replace("_", "\\_") + "}}" + String.Format("& {4}Total &{0}&{2}&{3}&{1:0.0}\\\\", n, (time_efca - time_start).TotalSeconds, lps, wits, "", cellcolor));
                        res.AddLast(String.Format("\\cmidrule{4} & {5}$1_L$ &{5}{0}&{5}{2}&{5}{3}&{5}{1:0.0}\\\\", reactions_unblocked, (time_unblocked - time_start).TotalSeconds, lps_unblocked, wits_unblocked, "{2-6}", cellcolor));
                        res.AddLast(String.Format("\\cmidrule{4} & $\\Coupling$ &{0}&{2}&{3}&{1:0.0}\\\\", couples, (time_couples - time_unblocked).TotalSeconds, lps_couples, wits_couples, "{2-6}", cellcolor));
            if(efca!=null)
                res.AddLast(String.Format("\\cmidrule{4} & {5}EFCA &{5}{1}&{5}{2}&{5}{3}&{5}{1:0.0}\\\\", efca.Length, (time_efca - time_couples).TotalSeconds, lps_efca, wits_efca, "{2-6}", cellcolor));
            
        }

        public static bool[] FindTrans(double[,] matrix, bool[] rev)
        {
            bool[] trans = new bool[rev.Length];
            int pos, neg;
            for (int i = 0; i < rev.Length; ++i)
            {
                pos = 0; neg = 0;
                for (int j = 0; j < matrix.Length / rev.Length; ++j)
                {
                    if (matrix[j, i] > 0)
                        ++pos;
                    if (matrix[j, i] < 0)
                        ++neg;
                }
                trans[i] = (pos * neg == 0) && (pos + neg > 0);
                //trans[i] = false;
            }
            return trans;
        }

    }

}