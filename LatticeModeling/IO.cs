using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using Tools.Converter;

namespace Tools.IO
{

    public class Writer
    {
        public static void WriteFile(ICollection<String> rows, String path)
        {

            using (StreamWriter sw = File.CreateText(path))
            {
                foreach (string line in rows)
                    sw.WriteLine(line);
            }
        }
    }

    public class Reader
    {

        public static LinkedList<string> ReadFile(string path)
        {
            LinkedList<string> file = new LinkedList<string>();

            using (StreamReader sr = new StreamReader(path, Encoding.Default))
            {
                while (sr.Peek() != -1)
                    file.AddLast(sr.ReadLine().Trim());
            }

            return file;
        }

        public static LinkedList<bool[]> ReadBoolFile(string path, string separator)
        {

            LinkedList<string> file = ReadFile(path);
            LinkedList<bool[]> res = new LinkedList<bool[]>();

            bool[] now;
            string[] splitted_line;
            foreach (string line in file)
            {
                splitted_line = line.Split(separator.ToCharArray());
                now = new bool[splitted_line.Length];

                for (int i = 0; i < now.Length; ++i)
                    now[i] = int.Parse(splitted_line[i]) != 0;

                res.AddLast(now);
            }

            return res;
        }

        public static List<int>[] ReadIntFile(string path, string separator)
        {

            LinkedList<string> file = ReadFile(path);
            List<int>[] res = new List<int>[file.Count];

            List<int> now;
            string[] splitted_line;
            int k = 0;
            foreach (string line in file)
            {
                now = new List<int>();
                splitted_line = line.Split(separator.ToCharArray());

                foreach (string number in splitted_line)
                    now.Add(int.Parse(number) - 1);
                now.Sort();

                res[k++] = now;
            }

            return res;
        }

        public static List<double>[] ReadDoubleFileToList(string path, string separator)
        {

            LinkedList<string> file = ReadFile(path);
            List<double>[] res = new List<double>[file.Count];

            List<double> now;
            string[] splitted_line;
            int k = 0;
            foreach (string line in file)
            {
                now = new List<double>();
                splitted_line = line.Trim().Split(separator.ToCharArray());

                foreach (string number in splitted_line)
                    now.Add(double.Parse(number));

                res[k++] = now;
            }

            return res;
        }
    }
}
