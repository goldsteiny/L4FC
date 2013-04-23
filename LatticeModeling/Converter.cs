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
