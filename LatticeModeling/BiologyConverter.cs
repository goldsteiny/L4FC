using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Tools.IO;


namespace Tools.Converter.Biology
{


    public class NetworkReader
    {
        public static void ReadDir(string path, string separator, out double[,] matrix, out bool[] rev, out String[] metabolite_names, out String[] reaction_names, out LinkedList<bool[]> forbidden_combinations)
        {
            LinkedList<String> mets_list = Reader.ReadFile(path + "\\mets.txt"),
                rxns_list = Reader.ReadFile(path + "\\rxns.txt");

            List<double>[] S_list = Reader.ReadDoubleFileToList(path + "\\S.txt", separator), 
                lb_list = Reader.ReadDoubleFileToList(path + "\\lb.txt", separator);



            int m = S_list.Length, n = lb_list.Length;

            metabolite_names = new String[m];
            int k = 0;
            foreach (string met in mets_list)
                metabolite_names[k++] = met;


            reaction_names = new String[n];
            k = 0;
            foreach (string r in rxns_list)
                reaction_names[k++] = r;

            rev = new bool[n];
            for (int i = 0; i < n; ++i)
                rev[i] = lb_list[i][0]<0;

            matrix = new double[m, n];
            for (int i = 0; i < m; ++i)
            {
                if ((10 * i) % m == 0)
                    Console.WriteLine("Read {0} of {1} metabolites.", i, m);

                for (int j = 0; j < n; ++j)
                    matrix[i, j] = S_list[i][j];
            }

            forbidden_combinations = new LinkedList<bool[]>();
            try
            {
                LinkedList<bool[]> cycles = Reader.ReadBoolFile(path + "\\cycles.txt", " ");
                foreach (bool[] cycle in cycles)
                    forbidden_combinations.AddLast(cycle);
            }
            catch (Exception) {}
        }

        public static void ReadMetatool(LinkedList<String> file, out double[,] matrix, out bool[] rev, out String[] metabolite_names, out String[] reaction_names)
        {

            // read reactions
            string line;
            do
            {
                line = file.First.Value;
                file.RemoveFirst();
            } while (!line.Equals("-ENZREV"));
            line = file.First.Value;
            file.RemoveFirst();
            string[] ENZREV = line.Length > 0 ? line.Split(new char[] { ' ' }) : new String[0];

            do
            {
                line = file.First.Value;
                file.RemoveFirst();
            } while (!line.Equals("-ENZIRREV"));
            line = file.First.Value;
            file.RemoveFirst();
            string[] ENZIRREV = line.Length > 0 ? line.Split(new char[] { ' ' }) : new String[0];

            Dictionary<string, int> rev_reactions = new Dictionary<string, int>();
            for (int i = 0; i < ENZREV.Length; ++i)
                rev_reactions.Add(ENZREV[i], i);

            Dictionary<string, int> irr_reactions = new Dictionary<string, int>();
            for (int i = 0; i < ENZIRREV.Length; ++i)
                irr_reactions.Add(ENZIRREV[i], i);

            // read metabolites
            do
            {
                line = file.First.Value;
                file.RemoveFirst();
            } while (!line.Equals("-METINT"));
            string[] METINT = file.First.Value.Split(new char[] { ' ' });
            file.RemoveFirst();

            do
            {
                line = file.First.Value;
                file.RemoveFirst();
            } while (!line.Equals("-METEXT"));
            string[] METEXT = file.First.Value.Split(new char[] { ' ' });
            file.RemoveFirst();

            Dictionary<string, int> metabolites = new Dictionary<string, int>();
            for (int i = 0; i < METINT.Length; ++i)
                metabolites.Add(METINT[i], i);
            for (int i = 0; i < METEXT.Length; ++i)
                metabolites.Add(METEXT[i], -1);

            double[,] A_rev = new double[METINT.Length, ENZREV.Length], A_irr = new double[METINT.Length, ENZIRREV.Length];

            bool[] transit_rev = new bool[ENZREV.Length], transit_irr = new bool[ENZIRREV.Length];

            // read stochiometrics
            do
            {
                line = file.First.Value;
                file.RemoveFirst();
            } while (!line.Equals("-CAT"));

            string[] lr, comps, nac;
            int k = -1, l = -1;
            string met, reaction;
            bool transit;
            int signum;
            double coeff;
            while (file.Count > 0)
            {
                line = file.First.Value;
                file.RemoveFirst();
                k = line.IndexOf(" : ");
                l = line.IndexOf(" = ");
                reaction = line.Substring(0, k);
                lr = line.Substring(k + 3, line.Length - k - 5).Split(new char[] { '=' });
                for (int right = 0; right < lr.Length; ++right)
                {
                    signum = (lr.Length > 1 ? (2 * right - 1) : line.Contains("= .") ? 1 : -1);
                    try
                    {
                        comps = lr[right].Trim().Split(new char[] { '+' });
                    }
                    catch (Exception)
                    {
                        comps = new string[0];
                    }
                    foreach (string comp in comps)
                    {
                        if (comp.Length < 1)
                            continue;
                        transit = true;
                        nac = comp.Trim().Split(new char[] { ' ' });
                        met = nac[nac.Length - 1];
                        try
                        {
                            coeff = nac.Length > 1 ? double.Parse(nac[0]) : 1;
                        }
                        catch (Exception)
                        {
                            coeff = 1;
                        }
                        metabolites.TryGetValue(met, out k);
                        if (rev_reactions.ContainsKey(reaction))
                            rev_reactions.TryGetValue(reaction, out l);
                        else
                            irr_reactions.TryGetValue(reaction, out l);
                        if (k >= 0)
                        {
                            if (rev_reactions.ContainsKey(reaction))
                                A_rev[k, l] = signum * coeff;
                            else
                                A_irr[k, l] = signum * coeff;
                            transit = false;
                        }
                        if (transit || lr.Length < 2) if (rev_reactions.ContainsKey(reaction))
                                transit_rev[l] = true;
                            else
                                transit_irr[l] = true;
                    }
                }
            }

            int m = METINT.Length, n_rev = ENZREV.Length, n_irr = ENZIRREV.Length, n = n_rev + n_irr;
            matrix = new double[m, n];
            rev = new bool[n];
            for (int j = 0; j < n_rev; ++j)
            {
                for (int i = 0; i < m; ++i)
                    matrix[i, j] = A_rev[i, j];
                rev[j] = true;
            }
            for (int j = 0; j < n_irr; ++j)
                for (int i = 0; i < m; ++i)
                    matrix[i, j + n_rev] = A_irr[i, j];


            metabolite_names = METINT;
            reaction_names = new String[n];
            for (int i = 0; i < n_rev; ++i)
                reaction_names[i] = ENZREV[i];
            for (int i = 0; i < n_irr; ++i)
                reaction_names[i + n_rev] = ENZIRREV[i];


        }

        public static LinkedList<String> NetworkToConnections(double[,] matrix, bool[] rev, String[] metabolite_names = null, String[] reaction_names = null)
        {

            if (metabolite_names == null)
            {
                metabolite_names = new String[matrix.Length / rev.Length];
                for (int i = 0; i < matrix.Length / rev.Length; ++i)
                    metabolite_names[i] = "M_" + (i + 1);
            }

            if (reaction_names == null)
            {
                reaction_names = new String[rev.Length];
                for (int i = 0; i < rev.Length; ++i)
                    reaction_names[i] = "R_" + (i + 1);
            }

            LinkedList<String> res = new LinkedList<string>();

            res.AddFirst("digraph Network {");
            res.AddLast("\tnodesep=0.7;");
            res.AddLast("\trankdir=LR;");

            int n = rev.Length, m = matrix.Length / n;

            for (int i = 0; i < m; ++i)
                res.AddLast(String.Format("\t{0}[label=\"{1}\"];", i + 1, metabolite_names[i]));


            LinkedList<int> produced = new LinkedList<int>(), consumed = new LinkedList<int>();
            for (int j = 0; j < n; ++j)
            {
                produced.Clear();
                consumed.Clear();

                for (int i = 0; i < m; ++i)
                {
                    if (matrix[i, j] < 0)
                        consumed.AddLast(i);
                    if (matrix[i, j] > 0)
                        produced.AddLast(i);
                }

                foreach (int m1 in produced)
                    foreach (int m2 in consumed)
                        if (rev[j])
                            res.AddLast(String.Format("{0}->{1}[dir=both]", m1 + 1, m2 + 1));
                        else
                            res.AddLast(String.Format("{0}->{1}", m1 + 1, m2 + 1));
            }


            res.AddLast("}");

            return res;
        }

        public static LinkedList<String> NetworkToDot(double[,] matrix, bool[] rev, String[] metabolite_names = null, String[] reaction_names = null)
        {
            if (metabolite_names == null)
            {
                metabolite_names = new String[matrix.Length/rev.Length];
                for (int i = 0; i < matrix.Length / rev.Length; ++i)
                    metabolite_names[i] = "M_" + (i + 1);
            }

            if (reaction_names == null)
            {
                reaction_names = new String[rev.Length];
                for (int i = 0; i < rev.Length; ++i)
                    reaction_names[i] = "R_" + (i + 1);
            }

            LinkedList<String> res = new LinkedList<string>();

            res.AddFirst("digraph Network {");
            res.AddLast("\tnodesep=0.5;");
            res.AddLast("\trankdir=LR;");

            int n = rev.Length, m = matrix.Length / n;

            for (int i = 0; i < m; ++i)
                res.AddLast(String.Format("\t{0}[label=\"{1}\"];", i + 1, metabolite_names[i]));

            int is_consumer, is_producer;
            bool[] detailed = new bool[n];
            for (int i = 0; i < n; ++i)
            {
                is_consumer = 0;
                is_producer = 0;
                for (int k = 0; k < m; ++k)
                {
                    if (matrix[k, i] > 0)
                        is_producer++;
                    if (matrix[k, i] < 0)
                        is_consumer++;
                }
                detailed[i] = is_producer > 1 || is_consumer > 1;
                if (detailed[i])
                {
                    res.AddLast(String.Format("\t{0}[label=\"\", height=\"0.01\", width=\"0.01\"];", m + i + 1));
                    res.AddLast(String.Format("\t{0}[label=\"\", height=\"0.01\", width=\"0.01\"];", m + n + i + 1));
                    if (is_producer > 0 && (is_consumer > 0 || !rev[i]))
                        res.AddLast(String.Format("\t{0}->{1}[label=\"{2}\", arrowhead=\"none\"];", m + i + 1, m + n + i + 1, reaction_names[i]));
                    else if (is_producer == 0)
                        res.AddLast(String.Format("\t{0}->{1}[label=\"{2}\"];", m + i + 1, m + n + i + 1, reaction_names[i]));
                    else
                        res.AddLast(String.Format("\t{0}->{1}[label=\"{2}\", dir=back];", m + i + 1, m + n + i + 1, reaction_names[i]));
                }
            }

            bool fwd_direction;
            for (int j = 0; j < n; ++j)
                if (detailed[j])
                {
                    for (int i = 0; i < m; ++i)
                        if (matrix[i, j] != 0)
                        {
                            fwd_direction = matrix[i, j] > 0;
                            if (rev[j])
                                if (fwd_direction)
                                    res.AddLast(String.Format("\t{0}->{1}[headlabel=\"{2:0.#}\"]", m + j + 1 + (fwd_direction ? n : 0), i + 1, Math.Abs(matrix[i, j])));
                                else
                                    res.AddLast(String.Format("\t{1}->{0}[taillabel=\"{2:0.#}\", dir=back]", m + j + 1 + (fwd_direction ? n : 0), i + 1, Math.Abs(matrix[i, j])));
                            else
                                if (fwd_direction)
                                    res.AddLast(String.Format("\t{0}->{1}[headlabel=\"{2:0.#}\"]", m + n + j + 1, i + 1, Math.Abs(matrix[i, j])));
                                else
                                    res.AddLast(String.Format("\t{1}->{0}[arrowhead=\"none\", taillabel=\"{2:0.#}\"]", m + j + 1, i + 1, Math.Abs(matrix[i, j])));
                        }
                }
                else
                {
                    int p = -1; int c = -1;
                    for (int i = 0; i < m; ++i)
                        if (matrix[i, j] > 0)
                            p = i;
                        else if (matrix[i, j] < 0)
                            c = i;
                    if (c > -1)
                    {
                        if (p > -1)
                            res.AddLast(String.Format("\t{0}->{1}[headlabel=\"{4:0.#}\", taillabel=\"{3:0.#}\", label={2}{5}]", c + 1, p + 1, reaction_names[j], -matrix[c, j], matrix[p, j], rev[j] ? ", dir=both" : ""));
                        else
                        {
                            res.AddLast(String.Format("\t{0}[label=\"\", height=\"0.01\", width=\"0.01\"];", m + j + 1));
                            res.AddLast(String.Format("\t{0}->{1}[taillabel=\"{3:0.#}\", label={2}{4}]", c + 1, m + j + 1, reaction_names[j], -matrix[c, j], rev[j] ? ", dir=both" : ""));
                        }
                    }
                    else if (p > -1)
                    {
                        res.AddLast(String.Format("\t{0}[label=\"\", height=\"0.01\", width=\"0.01\"];", m + j + 1));
                        res.AddLast(String.Format("\t{0}->{1}[headlabel=\"{3:0.#}\", label={2}{4}]", m + j + 1, p + 1, reaction_names[j], matrix[p, j], rev[j] ? ", dir=both" : ""));
                    }

                }

            res.AddLast("}");

            return res;
        }

        public static LinkedList<String> NetworkToMatlab(double[,] matrix, bool[] rev)
        {
            LinkedList<string> res = new LinkedList<string>(), temp_list;
            res.AddFirst("function [S, rev] = Network()");

            int n = rev.Length, m = matrix.Length / n;

            List<int> rev_list = new List<int>(n);
            for (int i = 0; i < n; ++i)
                rev_list.Add(rev[i] ? 1 : 0);

            LinkedList<List<int>> rev_listlist = new LinkedList<List<int>>();
            rev_listlist.AddFirst(rev_list);
            temp_list = ConvertingRoutines.List2Strings<List<int>, int>(rev_listlist, ", ", "");
            res.AddLast(String.Format("rev = [{0}];", temp_list.First.Value));


            List<List<double>> network = new List<List<double>>(m);
            for (int i = 0; i < m; ++i)
            {
                network.Add(new List<double>(n));
                for (int j = 0; j < n; ++j)
                    network[i].Add(matrix[i, j]);
            }
            temp_list = ConvertingRoutines.List2Strings<List<double>, double>(network, ", ", ";");

            string temp;
            temp = temp_list.First.Value;
            temp_list.RemoveFirst();
            res.AddLast("S = [" + temp);

            temp = temp_list.Last.Value;
            temp_list.RemoveLast();

            foreach (string row in temp_list)
                res.AddLast(row);
            res.AddLast(temp.Replace(";", "];"));

            return res;
        }



        public static void ExtractNetwork(List<double>[] network, out double[,] matrix, out bool[] rev)
        {
            int m = network.Length - 1, n = network[0].Count;

            rev = new bool[n];
            for (int i = 0; i < n; ++i)
                rev[i] = network[0].ElementAt(i) > 0;

            matrix = new double[m, n];
            for (int i = 0; i < m; ++i)
            {
                if ((10 * i) % m == 0)
                    Console.WriteLine("Read {0} of {1} metabolites.", i, m);

                for (int j = 0; j < n; ++j)
                    matrix[i, j] = network[i + 1].ElementAt(j);
            }
        }
    }
}
