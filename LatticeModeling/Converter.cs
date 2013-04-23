using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Tools.IO;


namespace Tools.Converter
{
    public class ConvertingRoutines
    {
        public static LinkedList<string> List2Strings<S, T>(ICollection<S> values, string separator, string separator2) where S : ICollection<T>
        {
            LinkedList<string> res = new LinkedList<string>();
            StringBuilder str = new StringBuilder();
            foreach (S row in values)
            {
                str.Clear();
                foreach (T value in row)
                    str.Append(value).Append(separator);
                res.AddLast(str.ToString(0, str.Length - separator.Length) + separator2);
            }
            return res;
        }

        public static LinkedList<string> List2Strings<S, T>(ICollection<S> values, string separator) where S : ICollection<T>
        {
            return List2Strings<S, T>(values, separator, "");
        }
    }
}
