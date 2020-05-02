//--------------------------------------------------------------------------------
// <copyright file="ITaskExecutor.cs" 
//            company="The University of Queensland"
//            author="Timothy O'Connor">
//     Copyright © The University of Queensland, 2012-2014. All rights reserved.
// </copyright>
// License: 
//--------------------------------------------------------------------------------

namespace Genomics
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Shared;

    /// <summary>
    /// Interface for items that can be executed as a task given a set of command line arguments
    /// </summary>
    public abstract class ITaskExecutor
    {
        /// <summary>
        /// Gets or sets the command arguments.
        /// </summary>
        /// <value>The command arguments.</value>
        public Args CommandArgs { get; set; }

        /// <summary>
        /// Gets the name of the class implementing the executed task.
        /// </summary>
        /// <value>The name.</value>
        public abstract string Name { get; }

        /// <summary>
        /// Gets the description of the task.
        /// </summary>
        /// <value>The description.</value>
        public abstract string Description { get; }

        /// <summary>
        /// Gets the command line options of the task.
        /// </summary>
        /// <value>The options.</value>
        public abstract Dictionary<string, string> Options { get; }

        /// <summary>
        /// Gets the help string.
        /// </summary>
        /// <value>The help string.</value>
        public string HelpString
        {
            get
            {
                return string.Format(
                    "\n\n{0}\n\n{1}\n\n\t{2}\n\n", 
                    this.Name,
                    this.Description,
                    string.Join("\n\t", this.Options.Select(x => string.Format("{0,-20}\t{1}", x.Key, x.Value))));
            }
        }

        /// <summary>
        /// Registers the task.
        /// </summary>
        /// <returns>The task.</returns>
        /// <typeparam name="TExecutor">Type of the task to be registered.</typeparam>
        public static KeyValuePair<string, Type> RegisterTask<TExecutor>() where TExecutor : ITaskExecutor
        {
            var e = (ITaskExecutor)System.Activator.CreateInstance<TExecutor>();

            return new KeyValuePair<string, Type>(e.Name, typeof(TExecutor));
        }

        /// <summary>
        /// Execute the analysis for which this class is the factor with the given command line arguments.
        /// </summary>
        /// <param name="commandArgs">Command arguments.</param>
        public abstract void Execute(Args commandArgs);
    }
}