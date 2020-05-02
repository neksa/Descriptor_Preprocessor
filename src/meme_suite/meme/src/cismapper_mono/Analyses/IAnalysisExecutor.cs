//--------------------------------------------------------------------------------
// <copyright file="IAnalysisExecutor.cs" 
//            company="The University of Queensland"
//            author="Timothy O'Connor">
//     Copyright © The University of Queensland, 2012-2014. All rights reserved.
// </copyright>
// License: 
//--------------------------------------------------------------------------------

namespace Analyses
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Genomics;
    using Shared;

    /// <summary>
    /// Analysis executor interface
    /// </summary>
    /// <typeparam name="TAnalysis">Type of the analysis to be executed.</typeparam>
    /// <typeparam name="TArgEnum">Type of the command line arguments to be used.</typeparam>
    public abstract class IAnalysisExecutor<TAnalysis, TArgEnum> : ITaskExecutor
    {
        /// <summary>
        /// Gets the command line options
        /// </summary>
        override public Dictionary<string, string> Options 
        {
            get
            {
                return this.OptionsData.ToDictionary(x => x.Key.ToString(), x => x.Value);
            }
        }

        /// <summary>
        /// Gets the name of the analysis.
        /// </summary>
        override public string Name
        {
            get
            {
                return typeof(TAnalysis).Name;
            }
        }

        /// <summary>
        /// Gets the command line option descriptions.
        /// </summary>
        protected abstract Dictionary<TArgEnum, string> OptionsData { get; }

        /// <summary>
        /// Gets the optional string argument.
        /// </summary>
        /// <returns>The optional string argument.</returns>
        /// <param name="args">Arguments.</param>
        /// <param name="arg">Argument.</param>
        protected string GetOptionalStringArg(Args args, TArgEnum arg)
        {
            return args.StringEnumArgs.ContainsKey(arg.ToString()) ? 
                args.StringEnumArgs[arg.ToString()] : 
                null;
        }

        /// <summary>
        /// Gets the optional int argument.
        /// </summary>
        /// <returns><c>true</c>, if optional int argument was gotten, <c>false</c> otherwise.</returns>
        /// <param name="args">Arguments.</param>
        /// <param name="arg">Argument.</param>
        /// <param name="i">The index.</param>
        protected bool GetOptionalIntArg(Args args, TArgEnum arg, ref int i)
        {
            if (args.StringEnumArgs.ContainsKey(arg.ToString()))
            {
                i = int.Parse(args.StringEnumArgs[arg.ToString()]);
                return true;
            }

            return false;
        }
        /// <summary>
        /// Gets the enum argument.
        /// </summary>
        /// <returns>The enum argument.</returns>
        /// <param name="args">Arguments.</param>
        /// <param name="arg">Argument.</param>
        /// <typeparam name="TEnum">The 1st type parameter.</typeparam>
        protected TEnum GetEnumArg<TEnum>(Args args, TArgEnum arg)
        {
            return (TEnum)Enum.Parse(typeof(TEnum), args.StringEnumArgs[arg.ToString()]);
        }

        /// <summary>
        /// Reflects the string.
        /// </summary>
        /// <param name="o">O.</param>
        /// <param name="arg">Argument.</param>
        protected void ReflectString(object o, TArgEnum arg)
        {
            this.GetProperty(o, arg).SetValue(o, CommandArgs.StringEnumArgs[arg.ToString()]);
        }

        /// <summary>
        /// Reflects the string arguments.
        /// </summary>
        /// <param name="analysis">Analysis.</param>
        /// <param name="args">Arguments.</param>
        protected void ReflectStringArgs(BaseAnalysis analysis, TArgEnum[] args)
        {
            var argEnumAsString = args.Select(x => x.ToString()).ToList();
            var registeredArgs = analysis.ElementArgRegistry
                .Where(x => argEnumAsString.Contains(x.Key)).ToList();

            foreach (var arg in registeredArgs)
            {
                arg.Value(this.CommandArgs.StringArgs[arg.Key]);
            }

            var unregisteredArgs = args.Where(x => !analysis.ElementArgRegistry.ContainsKey(x.ToString())).ToArray();

            this.ReflectStringArgs((object)analysis, unregisteredArgs);
        }

        /// <summary>
        /// Reflects the string arguments.
        /// </summary>
        /// <param name="analysis">Analysis.</param>
        /// <param name="args">Arguments.</param>
        protected void ReflectOptionalStringArgs(BaseAnalysis analysis, TArgEnum[] args)
        {
            var argEnumAsString = args.Select(x => x.ToString()).ToList();
            var registeredArgs = analysis.ElementArgRegistry
                .Where(x => argEnumAsString.Contains(x.Key)).ToList();

            foreach (var arg in registeredArgs)
            {
                if (this.CommandArgs.StringArgs.ContainsKey(arg.Key))
                {
                    arg.Value(this.CommandArgs.StringArgs[arg.Key]);
                }
            }

            var unregisteredArgs = args.Where(x => !analysis.ElementArgRegistry.ContainsKey(x.ToString())).ToArray();
                
            this.ReflectOptionalStringArgs((object)analysis, unregisteredArgs);
        }

        /// <summary>
        /// Reflects the string arguments.
        /// </summary>
        /// <param name="o">O.</param>
        /// <param name="args">Arguments.</param>
        protected void ReflectStringArgs(object o, TArgEnum[] args)
        {
            foreach (var arg in args)
            {
                this.ReflectString(o, arg);
            }
        }

        /// <summary>
        /// Reflects the optional string.
        /// </summary>
        /// <param name="o">O.</param>
        /// <param name="arg">Argument.</param>
        /// <param name="defaultValue">optional default value.</param>
        protected void ReflectOptionalString(object o, TArgEnum arg, string defaultValue)
        {
            string value = GetOptionalStringArg(this.CommandArgs, arg);

            if (value != null)
            {
                this.GetProperty(o, arg).SetValue(o, value);
            } 
            else if (defaultValue != null) 
            {
                this.GetProperty(o, arg).SetValue(o, defaultValue);
            }
        }

        /// <summary>
        /// Reflects the optional int.
        /// </summary>
        /// <param name="o">O.</param>
        /// <param name="arg">Argument.</param>
        /// <param name="defaultValue">Default value.</param>
        protected void ReflectOptionalInt(object o, TArgEnum arg, int defaultValue)
        {
            int value = this.CommandArgs.StringEnumArgs.ContainsKey(arg.ToString()) ?
                int.Parse(this.CommandArgs.StringEnumArgs[arg.ToString()]) :
                defaultValue;
            
            this.GetProperty(o, arg).SetValue(o, value);
        }

        /// <summary>
        /// Reflects the flag.
        /// </summary>
        /// <param name="o">O.</param>
        /// <param name="arg">Argument.</param>
        protected void ReflectFlag(object o, TArgEnum arg)
        {
            if (this.CommandArgs.Flags.Contains(arg.ToString()))
            {
                this.GetProperty(o, arg).SetValue(o, true);
            }
        }

        /// <summary>
        /// Reflects the optional string arguments.
        /// </summary>
        /// <param name="o">O.</param>
        /// <param name="args">Arguments.</param>
        protected void ReflectOptionalStringArgs(object o, TArgEnum[] args)
        {
            foreach (var arg in args)
            {
                this.ReflectOptionalString(o, arg, null);
            }
        }

        /// <summary>
        /// Gets the property.
        /// </summary>
        /// <returns>The property.</returns>
        /// <param name="o">O.</param>
        /// <param name="arg">Argument.</param>
        protected System.Reflection.PropertyInfo GetProperty(object o, TArgEnum arg)
        {
            System.Reflection.PropertyInfo prop = o.GetType().GetProperty(arg.ToString());
            if (prop == null)
            {
                throw new Exception(string.Format("Missing property for argument {0}", arg.ToString()));
            }

            return prop;
        }
    }
}
