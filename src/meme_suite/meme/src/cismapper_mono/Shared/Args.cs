//--------------------------------------------------------------------------------
// File: TRFScorer.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Shared
{
    using System;
    using System.Collections.Generic;
    using System.Linq;

    /// <summary>
    /// Command line argument object representation
    /// </summary>
    public class Args
    {
        private EnumDictionary<string> stringEnumArgs;

        /// <summary>
        /// Initializes a new instance of the <see cref="Shared.Args"/> class.
        /// </summary>
        /// <param name="args">Arguments.</param>
        public Args(string[] args)
        {
            int temp;
            StringArgs = new CheckLookupDictionary<string, string>(
                args.Select((x, i) => new {Arg = x, ValueIndex = i + 1})
                    .Where(x => 
                        this.IsArgName(x.Arg) &&
                         x.ValueIndex < args.Length &&
                        !this.IsArgName(args[x.ValueIndex]))
                    .Select(x => new { Arg = x.Arg, Value = args[x.ValueIndex] })
                    .ToDictionary(x => x.Arg.Substring(1), x => x.Value));
            
            IntArgs = new CheckLookupDictionary<string, int>(
                args.Select((x, i) => new {Arg = x, ValueIndex = i + 1})
                    .Where(x => 
                        this.IsArgName(x.Arg) &&
                        x.ValueIndex < args.Length &&
                        !this.IsArgName(args[x.ValueIndex]) && 
                        int.TryParse(args[x.ValueIndex], out temp))
                    .Select(x => new { Arg = x.Arg, Value = int.Parse(args[x.ValueIndex]) })
                    .ToDictionary(x => x.Arg.Substring(1), x => x.Value));
            
            Flags = new HashSet<string>(
                args.Select((x, i) => new {Arg = x, ValueIndex = i + 1})
                .Where(x => 
                    this.IsArgName(x.Arg) &&
                   (x.ValueIndex >= args.Length ||
                    this.IsArgName(args[x.ValueIndex])))
                .Select(x => x.Arg.Substring(1)));
        }

        /// <summary>
        /// Validate the specified argumentEnum.
        /// </summary>
        /// <param name="argumentEnum">Argument enum.</param>
        public void Validate(Type argumentEnum)
        {
            HashSet<string> names = new HashSet<string>(Enum.GetNames(argumentEnum).Concat(new string[] { "Mode", "Help" }));

            var invalidArgs = this.StringArgs.d.Keys.Where(arg => !names.Contains(arg))
                .Concat(this.IntArgs.d.Keys.Where(arg => !names.Contains(arg)))
                .Concat(this.Flags.Where(arg => !names.Contains(arg)))
                .ToList();

            if (invalidArgs.Count > 0)
            {
                throw new Exception("Invalid arguments found: " + string.Join(", ", invalidArgs));
            }
        }

        /// <summary>
        /// Gets the string arguments.
        /// </summary>
        /// <value>The string arguments.</value>
        public CheckLookupDictionary<string, string> StringArgs { get; private set; }

        /// <summary>
        /// Gets the string arguments.
        /// </summary>
        /// <value>The string arguments.</value>
        public EnumDictionary<string> StringEnumArgs 
        { 
            get
            {
                return Helpers.CheckInit(
                    ref this.stringEnumArgs,
                    () => new EnumDictionary<string>(this.StringArgs.d));
            }
        }

        /// <summary>
        /// Gets the int arguments.
        /// </summary>
        /// <value>The int arguments.</value>
        public CheckLookupDictionary<string, int> IntArgs { get; private set; }

        /// <summary>
        /// Gets the flags.
        /// </summary>
        /// <value>The flags.</value>
        public HashSet<string> Flags { get; private set; }

        /// <summary>
        /// Determines whether c is a digit.
        /// </summary>
        /// <returns><c>true</c> if this instance is digit the specified c; otherwise, <c>false</c>.</returns>
        /// <param name="c">C.</param>
        private bool IsDigit(char c)
        {
            int val;
            return int.TryParse(string.Format("{0}", c), out val);
        }

        /// <summary>
        /// Determines whether string is an argument name.
        /// </summary>
        /// <returns><c>true</c> if this instance is argument name the specified arg; otherwise, <c>false</c>.</returns>
        /// <param name="arg">Argument.</param>
        private bool IsArgName(string arg)
        {
            return arg.Length > 0 &&
                arg[0] == '-' &&
                arg.Length > 1 &&
                !this.IsDigit(arg[1]);
        }
    }
}

